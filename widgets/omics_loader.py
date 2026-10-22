"""Omics constellation loader: the animated view of a model run in progress.

One card with four parts: the run's input files, a small illustrative graph
(your sample in the centre, a ring of genes, an outer ring of drugs), a
console showing the backend's real log lines, and the current stage with a
five-part progress bar.

Everything that moves is painted by hand and advanced from a single ~60 fps
clock (`OmicsLoader._tick`). Each animated value eases toward a target that
depends only on the current step, so a step change never jumps: whatever is
on screen glides from where it is to where it should be.

The graph is an illustration, not the data: 12 well-known genes and 9 GDSC2
drugs with their real targets, standing in for the 16,214 proteins and ~500
drugs of `hetero_graph.pt`. The stage copy describes what `HeteroIC50GNN`
actually does (PCA omics features joined per cell line, Morgan fingerprints,
two message-passing layers).

Driving it:
    loader.start(files, cell_line)   # a new run
    loader.log(line)                 # each real backend log line
    loader.set_progress(fraction)    # the backend's own progress, 0..1
    loader.complete() / loader.fail(message)
    loader.finished                  # emitted once the "Done" stage has shown

The steps follow the run but never outrun it: a step is only entered once the
backend's log shows the run has got that far, and each step is held for a
short minimum so its animation can play.
"""

from __future__ import annotations

import html
import math
import re
import time

from PySide6.QtCore import QEasingCurve, QPointF, QPropertyAnimation, QRectF, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QFontMetricsF, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import (
    QFrame,
    QGraphicsOpacityEffect,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QScrollArea,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from styles.theme import (
    BORDER,
    CARD_RADIUS,
    CONSOLE_BACKGROUND,
    CONSOLE_BORDER,
    CONSOLE_HEADER,
    CONSOLE_TEXT,
    CONSOLE_TEXT_MUTED,
    DIVIDER,
    LABEL_CAPS_STYLE,
    OMICS_EXPRESSION,
    OMICS_MUTATION,
    OMICS_PROTEOMICS,
    SELECTED_DRUG_SURFACE,
    SUCCESS_SOFT,
    SURFACE,
    SURFACE_LOW,
    TEXT,
    TEXT_MUTED,
    TEXT_SOFT,
    label_style,
)
from widgets.help import HELP

_MONO = "Consolas"
_SANS = "Segoe UI"

_ACCENT = QColor(SELECTED_DRUG_SURFACE)
_WHITE = QColor("#ffffff")
_INK = QColor(TEXT)
_MUTED = QColor(TEXT_MUTED)
_NODE_BORDER = QColor("#8a8d8d")
_QUIET = QColor("#5f6262")

# Omics channels, in the order their squares sit on a gene: letter, name,
# colour, and the step at which that layer is added.
_CHANNELS = [
    ("E", "expression", QColor(OMICS_EXPRESSION), 2),
    ("M", "mutation / CNV", QColor(OMICS_MUTATION), 3),
    ("P", "proteomics", QColor(OMICS_PROTEOMICS), 4),
]

# Illustrative slice of the graph: STRING proteins, interactions between
# them, and GDSC2 drugs with their targets. Drug angles are degrees on the
# outer ring (0 = right, clockwise).
_GENES = ["EGFR", "ERBB2", "PIK3CA", "PTEN", "AKT1", "MTOR", "MDM2", "TP53", "MAPK1", "MAP2K1", "BRAF", "KRAS"]
_PPI = [
    ("EGFR", "ERBB2"), ("EGFR", "KRAS"), ("EGFR", "PIK3CA"), ("ERBB2", "PIK3CA"), ("KRAS", "BRAF"),
    ("KRAS", "PIK3CA"), ("BRAF", "MAP2K1"), ("MAP2K1", "MAPK1"), ("PIK3CA", "AKT1"), ("PIK3CA", "PTEN"),
    ("PTEN", "AKT1"), ("AKT1", "MTOR"), ("AKT1", "MDM2"), ("MDM2", "TP53"), ("TP53", "PTEN"),
]
_DRUGS = [
    ("Erlotinib", ["EGFR"], -103), ("Lapatinib", ["EGFR", "ERBB2"], -66), ("Alpelisib", ["PIK3CA"], -28),
    ("MK-2206", ["AKT1"], 28), ("AZD8055", ["MTOR"], 62), ("Nutlin-3a", ["MDM2"], 92),
    ("Ulixertinib", ["MAPK1"], 150), ("Trametinib", ["MAP2K1"], 180), ("Dabrafenib", ["BRAF"], 210),
]

# (heading, description) per step.
_STAGES = [
    ("Starting up", "Picking the compute device."),
    ("Loading the knowledge graph", "Cell lines, drugs and proteins, with the STRING interactions between them."),
    ("Adding gene expression", "E fills in: your sample's expression profile, compressed with PCA."),
    ("Adding mutations and copy number", "M fills in: the mutation and copy-number profile."),
    ("Adding proteomics", "P fills in: the protein-abundance profile."),
    ("Adding the drugs", "Each drug joins the proteins it targets; its structure is a Morgan fingerprint."),
    ("Combining the omics for your sample", "The three profiles are joined into one feature vector for your sample."),
    ("Passing messages through the graph", "Two graph layers pass signal along the protein interactions."),
    ("Predicting drug response", "Each drug gets a predicted ln(IC50) for your sample."),
    ("Done", "Opening your results."),
]
_LAST_STEP = len(_STAGES) - 1
# (label, first step, last step) of the five progress segments.
_SEGMENTS = [("Graph", 0, 1), ("Omics", 2, 4), ("Drugs", 5, 5), ("Model", 6, 7), ("Predict", 8, 9)]

# Shortest time a step stays on screen, so its animation can play even when
# the backend races through it. Scaled down once the run has finished and the
# view is only catching up.
_MIN_DWELL = [0.5, 1.2, 0.9, 0.9, 0.9, 1.2, 1.3, 1.5, 1.3]
_CATCH_UP = 0.45
# How long "Done" shows before `finished` is emitted.
_DONE_HOLD = 0.9

# (role, key into `files`, step at which it is read, glyph, HELP key)
_FILES = [
    ("Graph", "graph", 1, "graph", "file_graph"),
    ("Expression", "expression", 2, "E", "omics_transcriptomics"),
    ("Mutation / CNV", "mutation", 3, "M", "omics_genomics"),
    ("Proteomics", "proteomics", 4, "P", "omics_proteomics"),
    ("Drug response", "response", 5, "drug", "file_dataset"),
    ("Trained model", "model", 6, "model", "file_checkpoint"),
]

_TAG_COLOURS = {
    "setup": "#a1a1aa", "graph": "#a1a1aa", "data": "#93c5fd", "cell-line": "#93c5fd", "predict": "#93c5fd",
    "model": "#f4f4f5", "eval": "#f4f4f5", "epoch": "#f4f4f5", "train": "#f4f4f5",
    "done": "#4ade80", "error": "#f87171",
}
_LOG_TAG = re.compile(r"^\[([a-z\-]+)[^\]]*\]\s*(.*)$")


def _font(family: str, px: float, weight: QFont.Weight = QFont.Weight.Normal, spacing: float = 0.0) -> QFont:
    font = QFont(family)
    font.setPixelSize(round(px))
    font.setWeight(weight)
    if spacing:
        font.setLetterSpacing(QFont.SpacingType.AbsoluteSpacing, spacing)
    return font


def _mix(a: QColor, b: QColor, t: float) -> QColor:
    t = max(0.0, min(1.0, t))
    return QColor.fromRgbF(
        a.redF() + (b.redF() - a.redF()) * t,
        a.greenF() + (b.greenF() - a.greenF()) * t,
        a.blueF() + (b.blueF() - a.blueF()) * t,
    )


def _alpha(colour: QColor, alpha: float) -> QColor:
    out = QColor(colour)
    out.setAlphaF(max(0.0, min(1.0, alpha)))
    return out


def _ease(current: float, target: float, dt: float, tau: float) -> float:
    """Move `current` toward `target`, covering ~63% of the gap every `tau` seconds.

    Lands exactly on the target once close: a node at scale 0.9997 is drawn
    through the slow scaled-text path, one at exactly 1.0 is not.
    """
    value = target + (current - target) * math.exp(-dt / tau)
    return target if abs(value - target) < 0.002 else value


def _draw_check(painter: QPainter, box: QRectF, colour: QColor) -> None:
    s = box.width() / 24.0
    pen = QPen(colour, 3 * s, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin)
    painter.setPen(pen)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawPolyline(QPolygonF([
        QPointF(box.x() + 5 * s, box.y() + 12.5 * s),
        QPointF(box.x() + 9.5 * s, box.y() + 17 * s),
        QPointF(box.x() + 19 * s, box.y() + 7.5 * s),
    ]))


def _draw_pip(painter: QPainter, box: QRectF, letter: str, colour: QColor, fill: float, font: QFont) -> None:
    """An omics square: outlined when empty, solid once that layer is loaded."""
    painter.setPen(QPen(colour, 1.5))
    painter.setBrush(_mix(_WHITE, colour, fill))
    painter.drawRoundedRect(box, 3, 3)
    painter.setPen(_mix(colour, _WHITE, fill))
    painter.setFont(font)
    painter.drawText(box, Qt.AlignmentFlag.AlignCenter, letter)


class _Constellation(QWidget):
    """The graph: your sample, a ring of genes, an outer ring of drugs."""

    selected = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumSize(420, 390)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMouseTracking(True)

        self._step = -1
        self._clock = 0.0
        self._cell_line = ""
        self._drug_total: int | None = None
        self._selection: str | None = None
        # Seconds since each staggered group was switched on (-1 while off).
        self._gates: dict[str, float] = {}
        # (id, bounding box, is_round) of every clickable node, from the last paint.
        self._hits: list[tuple[str, QRectF, bool]] = []

        index = {name: i for i, name in enumerate(_GENES)}
        # Positions in percent of the widget, as in a 100 x 100 box.
        self._gene_at = [
            (50 + 25 * math.cos(math.radians(-90 + 30 * i)), 50 + 30 * math.sin(math.radians(-90 + 30 * i)))
            for i in range(len(_GENES))
        ]
        self._drug_at = [
            (50 + 40 * math.cos(math.radians(angle)), 50 + 43 * math.sin(math.radians(angle)))
            for _name, _targets, angle in _DRUGS
        ]

        # Interactions appear in clockwise order, so the ring draws itself.
        def clock_angle(pair: tuple[int, int]) -> float:
            (ax, ay), (bx, by) = self._gene_at[pair[0]], self._gene_at[pair[1]]
            return (math.degrees(math.atan2((ay + by) / 2 - 50, (ax + bx) / 2 - 50)) + 90) % 360

        self._ppi = sorted(((index[a], index[b]) for a, b in _PPI), key=clock_angle)
        self._targets = [(d, index[t]) for d, (_n, targets, _a) in enumerate(_DRUGS) for t in targets]

        self._font_gene = _font(_MONO, 11, QFont.Weight.Bold)
        self._font_pip = _font(_MONO, 9, QFont.Weight.Bold)
        self._font_drug = _font(_SANS, 12, QFont.Weight.DemiBold)
        self._font_sample = _font(_MONO, 13, QFont.Weight.Bold)
        self._font_cell = _font(_SANS, 11)
        self._font_note = _font(_MONO, 12)
        metrics = QFontMetricsF(self._font_drug)
        self._drug_text_width = [metrics.horizontalAdvance(name) for name, _t, _a in _DRUGS]

        self.reset("")

    # -- state ----------------------------------------------------------

    def reset(self, cell_line: str) -> None:
        """Snap back to the empty graph for a new run."""
        genes, drugs = len(_GENES), len(_DRUGS)
        self._step = -1
        self._cell_line = cell_line
        self._cell_name = ""
        self._drug_total = None
        self._selection = None
        self._gates.clear()
        self._gene_on = [0.0] * genes
        self._pip_fill = [[0.0] * len(_CHANNELS) for _ in range(genes)]
        self._pip_pop = [[0.0] * len(_CHANNELS) for _ in range(genes)]
        self._drug_on = [0.0] * drugs
        self._drug_scored = [0.0] * drugs
        self._ppi_on = [0.0] * len(self._ppi)
        self._spoke_on = [0.0] * genes
        self._target_on = [0.0] * len(self._targets)
        self._ppi_flow = self._spoke_flow = self._target_flow = 0.0
        self._ppi_phase = 0.0
        self._ppi_speed = 1 / 2.4
        self._spoke_rgb = [_ACCENT.redF(), _ACCENT.greenF(), _ACCENT.blueF()]
        self._ping = self._bars = self._pulse = self._note = 0.0
        self.update()

    def set_step(self, step: int) -> None:
        self._step = step

    def set_drug_total(self, total: int) -> None:
        self._drug_total = total

    def set_cell_name(self, name: str) -> None:
        self._cell_name = name

    def set_selection(self, node_id: str | None) -> None:
        self._selection = node_id
        self.update()

    def _since(self, key: str, on: bool, dt: float) -> float:
        """Seconds since `on` became true, or -1 while it is false."""
        if not on:
            self._gates[key] = -1.0
            return -1.0
        previous = self._gates.get(key, -1.0)
        self._gates[key] = 0.0 if previous < 0 else previous + dt
        return self._gates[key]

    def advance(self, dt: float) -> None:
        """Ease every animated value toward its target for the current step."""
        step = self._step
        self._clock += dt

        def stagger(values: list[float], since: float, gap: float, tau: float = 0.2) -> None:
            for i, value in enumerate(values):
                values[i] = _ease(value, 1.0 if since >= i * gap else 0.0, dt, tau)

        stagger(self._gene_on, self._since("genes", step >= 1, dt), 0.07)
        stagger(self._ppi_on, self._since("ppi", step >= 1, dt), 0.07, 0.22)
        stagger(self._spoke_on, self._since("spokes", step >= 2, dt), 0.05, 0.22)
        stagger(self._drug_on, self._since("drugs", step >= 5, dt), 0.08)
        stagger(self._target_on, self._since("targets", step >= 5, dt), 0.08, 0.22)
        stagger(self._drug_scored, self._since("scored", step >= 8, dt), 0.09, 0.16)

        for c, (_letter, _name, _colour, at_step) in enumerate(_CHANNELS):
            since = self._since(f"pip{c}", step >= at_step, dt)
            for g in range(len(_GENES)):
                filled = since >= g * 0.06
                self._pip_fill[g][c] = _ease(self._pip_fill[g][c], 1.0 if filled else 0.0, dt, 0.16)
                popping = filled and step == at_step
                self._pip_pop[g][c] = _ease(self._pip_pop[g][c], 1.0 if popping else 0.0, dt, 0.14)

        self._ppi_flow = _ease(self._ppi_flow, 0.9 if step == 7 else 0.18 if step >= 1 else 0.0, dt, 0.3)
        self._ppi_speed = _ease(self._ppi_speed, 1 / 0.8 if step == 7 else 1 / 2.4, dt, 0.4)
        self._ppi_phase += self._ppi_speed * dt
        self._spoke_flow = _ease(self._spoke_flow, 0.85 if 2 <= step <= 4 or step == 6 else 0.0, dt, 0.25)
        self._target_flow = _ease(self._target_flow, 0.9 if step == 8 else 0.0, dt, 0.28)

        spoke_colour = next((colour for _l, _n, colour, at_step in _CHANNELS if at_step == step), _ACCENT)
        for k, channel in enumerate((spoke_colour.redF(), spoke_colour.greenF(), spoke_colour.blueF())):
            self._spoke_rgb[k] = _ease(self._spoke_rgb[k], channel, dt, 0.25)

        self._ping = _ease(self._ping, 1.0 if 6 <= step <= 8 else 0.0, dt, 0.3)
        self._bars = _ease(self._bars, 1.0 if step >= 6 else 0.0, dt, 0.25)
        self._pulse = _ease(self._pulse, 1.0 if step == 7 else 0.0, dt, 0.3)
        self._note = _ease(self._note, 1.0 if step >= 5 else 0.0, dt, 0.25)
        self.update()

    # -- painting -------------------------------------------------------

    def paintEvent(self, _event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        width, height, clock = self.width(), self.height(), self._clock
        centre = QPointF(width / 2, height / 2)
        tau = 2 * math.pi

        # Nodes drift a few pixels on slow, out-of-step loops; edges follow.
        genes: list[QPointF] = []
        for i, (x, y) in enumerate(self._gene_at):
            period, on = 4 + (i % 3) * 0.8, self._gene_on[i]
            genes.append(QPointF(
                x * width / 100 + 2 * on * math.sin(tau * (clock / period + i * 0.17)),
                y * height / 100 + 3 * on * math.sin(tau * (clock / (period * 1.3) + i * 0.29)),
            ))
        drugs: list[QPointF] = []
        drug_widths: list[float] = []
        for i, (x, y) in enumerate(self._drug_at):
            period, on = 4.6 + (i % 3) * 0.7, self._drug_on[i]
            pill = 24 + self._drug_text_width[i] + 18 * self._drug_scored[i]
            drug_widths.append(pill)
            # Keep the whole pill inside the widget on a narrow window.
            px = min(max(x * width / 100, pill / 2 + 6), width - pill / 2 - 6)
            drugs.append(QPointF(
                px + 2 * on * math.sin(tau * (clock / period + i * 0.23)),
                y * height / 100 + 3 * on * math.sin(tau * (clock / (period * 1.3) + i * 0.31)),
            ))

        def lines(pairs, opacities, colour: QColor, pen_width: float, dashes=None) -> None:
            pen = QPen(colour, pen_width, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap)
            if dashes:
                pen.setCapStyle(Qt.PenCapStyle.FlatCap)
                pen.setDashPattern([d / pen_width for d in dashes])
            for (a, b), opacity in zip(pairs, opacities):
                if opacity > 0.01:
                    pen.setColor(_alpha(colour, opacity))
                    painter.setPen(pen)
                    painter.drawLine(a, b)

        def flow(pairs, opacity: float, colour: QColor, pen_width: float, dashes, phase: float) -> None:
            """Dashes travelling from each pair's first point to its second."""
            if opacity <= 0.01:
                return
            pen = QPen(_alpha(colour, opacity), pen_width, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap)
            pen.setDashPattern([d / pen_width for d in dashes])
            pen.setDashOffset((1 - phase % 1) * 2 * sum(dashes) / pen_width)
            painter.setPen(pen)
            for a, b in pairs:
                painter.drawLine(a, b)

        spokes = [(centre, gene) for gene in genes]
        ppi = [(genes[a], genes[b]) for a, b in self._ppi]
        # Drawn gene -> drug, so the prediction flows outward to each drug.
        targets = [(genes[g], drugs[d]) for d, g in self._targets]

        lines(spokes, self._spoke_on, QColor("#d4d6d6"), 1, [3, 4])
        flow(spokes, self._spoke_flow, QColor.fromRgbF(*self._spoke_rgb), 2, [4, 8], clock / 0.9)
        lines(ppi, self._ppi_on, QColor("#b9bcbc"), 1.5)
        flow(ppi, self._ppi_flow, _ACCENT, 2.5, [5, 7], self._ppi_phase)
        lines(targets, self._target_on, _NODE_BORDER, 1.25, [2, 4])
        flow(targets, self._target_flow, _ACCENT, 2, [4, 8], clock / 0.9)

        self._hits = []
        self._paint_sample(painter, centre)
        for i, point in enumerate(genes):
            self._paint_gene(painter, i, point)
        for i, point in enumerate(drugs):
            self._paint_drug(painter, i, point, drug_widths[i])

        if self._note > 0.01:
            painter.setOpacity(self._note)
            painter.setPen(QColor(TEXT_SOFT))
            painter.setFont(self._font_note)
            shown = f"{len(_DRUGS)} of {self._drug_total} drugs shown" if self._drug_total else f"{len(_DRUGS)} drugs shown"
            painter.drawText(QRectF(0, height - 18, width, 18), Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, shown)
        painter.end()

    def _paint_sample(self, painter: QPainter, centre: QPointF) -> None:
        clock, tau = self._clock, 2 * math.pi
        if self._ping > 0.01:
            progress = (clock / 1.8) % 1
            progress = 1 - (1 - progress) ** 2
            radius = 52 * (1 + 0.7 * progress)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(_alpha(_ACCENT, 0.35 * (1 - progress) * self._ping))
            painter.drawEllipse(centre, radius, radius)

        radius = 52 * (1 + 0.0225 * (1 - math.cos(tau * clock / 3.2)))
        selected = self._selection == "sample"
        painter.setPen(QPen(QColor("#8fd3f0") if selected else _ACCENT, 3))
        painter.setBrush(_ACCENT)
        painter.drawEllipse(centre, radius - 1.5, radius - 1.5)
        self._hits.append(("sample", QRectF(centre.x() - radius, centre.y() - radius, 2 * radius, 2 * radius), True))

        # The three bars slide in under the labels once the omics are combined.
        bars = self._bars
        top = centre.y() - (31 + 22 * bars) / 2
        painter.setPen(_WHITE)
        painter.setFont(self._font_sample)
        painter.drawText(QRectF(centre.x() - 52, top, 104, 14), Qt.AlignmentFlag.AlignCenter, self._cell_line)
        # The cell line's name, once the backend's log has reported it.
        painter.setPen(_alpha(_WHITE, 0.85))
        painter.setFont(self._font_cell)
        name = QFontMetricsF(self._font_cell).elidedText(self._cell_name, Qt.TextElideMode.ElideRight, 92)
        painter.drawText(QRectF(centre.x() - 50, top + 18, 100, 13), Qt.AlignmentFlag.AlignCenter, name)
        if bars > 0.01:
            wave = 0.5 - 0.5 * math.cos(tau * clock / 1.6)
            levels = (0.30 + 0.55 * wave, 0.70 - 0.30 * wave, 0.45 + 0.30 * wave)
            painter.setPen(Qt.PenStyle.NoPen)
            for k, (colour, level) in enumerate(zip(("#8fd3f0", "#f6b183", "#c9b3f2"), levels)):
                bar_height = 18 * level
                painter.setBrush(_alpha(QColor(colour), bars))
                painter.drawRoundedRect(QRectF(centre.x() - 13 + k * 10, top + 53 - bar_height, 6, bar_height), 1, 1)

    def _paint_gene(self, painter: QPainter, i: int, at: QPointF) -> None:
        on = self._gene_on[i]
        scale = 0.8 + 0.2 * on
        box_width, box_height = 64.0, 44.0
        box = QRectF(-box_width / 2, -box_height / 2, box_width, box_height)
        painter.save()
        painter.translate(at)
        if scale != 1.0:
            painter.scale(scale, scale)
        painter.setOpacity(0.12 + 0.88 * on)

        if self._pulse > 0.01:
            beat = math.sin(math.pi * ((self._clock / 1.2 + (i % 4) * 0.21) % 1)) ** 2
            spread = 6 * beat
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(_alpha(_ACCENT, 0.22 * beat * self._pulse))
            painter.drawRoundedRect(box.adjusted(-spread, -spread, spread, spread), 10 + spread, 10 + spread)

        selected = self._selection == _GENES[i]
        painter.setPen(QPen(_ACCENT, 2) if selected else QPen(_NODE_BORDER, 1.5))
        painter.setBrush(_WHITE)
        painter.drawRoundedRect(box, 10, 10)
        painter.setPen(_INK)
        painter.setFont(self._font_gene)
        painter.drawText(QRectF(box.x(), box.y() + 6, box_width, 13), Qt.AlignmentFlag.AlignCenter, _GENES[i])

        for c, (letter, _name, colour, _at_step) in enumerate(_CHANNELS):
            pop = 1 + 0.15 * self._pip_pop[i][c]
            painter.save()
            painter.translate(-17 + c * 17, box.y() + 30)
            if pop != 1.0:
                painter.scale(pop, pop)
            _draw_pip(painter, QRectF(-7, -7, 14, 14), letter, colour, self._pip_fill[i][c], self._font_pip)
            painter.restore()
        painter.restore()

        if on > 0.5:
            self._hits.append((_GENES[i], QRectF(at.x() - 32 * scale, at.y() - 22 * scale, 64 * scale, 44 * scale), False))

    def _paint_drug(self, painter: QPainter, i: int, at: QPointF, pill_width: float) -> None:
        on, scored = self._drug_on[i], self._drug_scored[i]
        if on <= 0.01:
            return
        name = _DRUGS[i][0]
        scale = 0.6 + 0.4 * on
        box = QRectF(-pill_width / 2, -16, pill_width, 32)
        painter.save()
        painter.translate(at)
        if scale != 1.0:
            painter.scale(scale, scale)
        painter.setOpacity(on)
        if self._selection == name:
            painter.setPen(QPen(_ACCENT, 2))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRoundedRect(box.adjusted(-3.5, -3.5, 3.5, 3.5), 19.5, 19.5)
        painter.setPen(QPen(_ACCENT, 1.5))
        painter.setBrush(_mix(_WHITE, _ACCENT, scored))
        painter.drawRoundedRect(box, 16, 16)
        if scored > 0.01:
            painter.setOpacity(on * scored)
            _draw_check(painter, QRectF(box.x() + 12, -6, 12, 12), _WHITE)
            painter.setOpacity(on)
        painter.setPen(_mix(_INK, _WHITE, scored))
        painter.setFont(self._font_drug)
        painter.drawText(
            QRectF(box.x() + 12 + 18 * scored, -16, self._drug_text_width[i] + 2, 32),
            Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
            name,
        )
        painter.restore()
        if on > 0.5:
            self._hits.append((name, QRectF(at.x() - pill_width / 2, at.y() - 16, pill_width, 32), False))

    # -- interaction ----------------------------------------------------

    def _node_at(self, point: QPointF) -> str | None:
        for node_id, box, is_round in reversed(self._hits):
            if not box.contains(point):
                continue
            if is_round:
                centre = box.center()
                if math.hypot(point.x() - centre.x(), point.y() - centre.y()) > box.width() / 2:
                    continue
            return node_id
        return None

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        over = self._node_at(event.position()) is not None
        self.setCursor(Qt.CursorShape.PointingHandCursor if over else Qt.CursorShape.ArrowCursor)

    def mousePressEvent(self, event) -> None:  # noqa: N802
        node_id = self._node_at(event.position())
        if node_id is not None:
            self.selected.emit(node_id)


class _FileCard(QWidget):
    """One input of the run: its role, where it comes from, and its load state."""

    def __init__(self, role: str, glyph: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._role, self._glyph = role, glyph
        self._name = ""
        self._state = "waiting"  # "waiting" | "active" | "done"
        self._active = 0.0
        self._done = 0.0
        self.setFixedHeight(78)
        self.setMinimumWidth(118)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._font_role = _font(_SANS, 12, QFont.Weight.Bold)
        self._font_name = _font(_MONO, 12)
        self._font_status = _font(_SANS, 12, QFont.Weight.DemiBold)
        self._font_pip = _font(_MONO, 10, QFont.Weight.Bold)

    def set_name(self, name: str) -> None:
        self._name = name
        self.update()

    def set_state(self, state: str, snap: bool = False) -> None:
        self._state = state
        if snap:
            self._active = 1.0 if state == "active" else 0.0
            self._done = 1.0 if state == "done" else 0.0
            self.update()

    def advance(self, dt: float) -> None:
        active = _ease(self._active, 1.0 if self._state == "active" else 0.0, dt, 0.14)
        done = _ease(self._done, 1.0 if self._state == "done" else 0.0, dt, 0.14)
        if abs(active - self._active) > 0.001 or abs(done - self._done) > 0.001:
            self._active, self._done = active, done
            self.update()

    def paintEvent(self, _event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        frame = QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5)
        painter.setPen(QPen(_mix(QColor(DIVIDER), _ACCENT, self._active), 1))
        painter.setBrush(_mix(_WHITE, QColor("#f2f6fa"), self._active))
        painter.drawRoundedRect(frame, 8, 8)

        self._paint_glyph(painter, QRectF(12, 11, 16, 16))
        text_width = self.width() - 24
        painter.setPen(_INK)
        painter.setFont(self._font_role)
        painter.drawText(QRectF(36, 10, text_width - 24, 18), Qt.AlignmentFlag.AlignVCenter, self._role)

        painter.setPen(_MUTED)
        painter.setFont(self._font_name)
        name = QFontMetricsF(self._font_name).elidedText(self._name, Qt.TextElideMode.ElideMiddle, text_width)
        painter.drawText(QRectF(12, 32, text_width, 16), Qt.AlignmentFlag.AlignVCenter, name)

        waiting = _mix(_QUIET, _ACCENT, self._active)
        painter.setFont(self._font_status)
        if self._done < 0.99:
            painter.setOpacity(1 - self._done)
            painter.setPen(waiting)
            label = "Reading…" if self._state == "active" else "Waiting"
            painter.drawText(QRectF(12, 53, text_width, 16), Qt.AlignmentFlag.AlignVCenter, label)
        if self._done > 0.01:
            painter.setOpacity(self._done)
            green = QColor(SUCCESS_SOFT)
            _draw_check(painter, QRectF(12, 55, 12, 12), green)
            painter.setPen(green)
            painter.drawText(QRectF(30, 53, text_width - 18, 16), Qt.AlignmentFlag.AlignVCenter, "Loaded")
        painter.end()

    def _paint_glyph(self, painter: QPainter, box: QRectF) -> None:
        channel = next((c for c in _CHANNELS if c[0] == self._glyph), None)
        if channel is not None:
            _draw_pip(painter, box, channel[0], channel[2], 1.0, self._font_pip)
        elif self._glyph == "graph":
            a, b, c = (QPointF(box.x() + x, box.y() + y) for x, y in ((3, 4), (12, 6), (6, 13)))
            painter.setPen(QPen(_NODE_BORDER, 1.5))
            painter.drawPolyline(QPolygonF([a, b, c, a]))
            painter.setPen(QPen(_QUIET, 1.4))
            painter.setBrush(_WHITE)
            for point in (a, b, c):
                painter.drawEllipse(point, 2.2, 2.2)
        elif self._glyph == "drug":
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(_ACCENT)
            painter.drawRoundedRect(QRectF(box.x() - 1, box.y() + 3, 18, 10), 5, 5)
        else:
            painter.setPen(QPen(_INK, 1.4))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRoundedRect(box.adjusted(2, 2, -2, -2), 2.5, 2.5)
            painter.drawRect(box.adjusted(5.5, 5.5, -5.5, -5.5))


class _Legend(QWidget):
    """Key to the graph's shapes; wraps onto more rows on a narrow card."""

    _ITEMS = [
        ("gene", "Gene in the STRING network"),
        ("pips", "Omics layer loaded"),
        ("drug", "GDSC2 drug and its targets"),
        ("ppi", "Protein interaction"),
        ("flow", "Signal in the model"),
    ]
    _ROW, _GAP, _SWATCH = 24, 20, {"gene": 22, "pips": 48, "drug": 22, "ppi": 22, "flow": 22}

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._font = _font(_SANS, 12)
        self._font_pip = _font(_MONO, 9, QFont.Weight.Bold)
        metrics = QFontMetricsF(self._font)
        self._widths = [self._SWATCH[kind] + 6 + metrics.horizontalAdvance(text) for kind, text in self._ITEMS]
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setFixedHeight(self._ROW)

    def _layout(self, width: int) -> list[tuple[float, float]]:
        x, y, places = 0.0, 0.0, []
        for item_width in self._widths:
            if x > 0 and x + item_width > width:
                x, y = 0.0, y + self._ROW
            places.append((x, y))
            x += item_width + self._GAP
        return places

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self.setFixedHeight(int(self._layout(self.width())[-1][1]) + self._ROW)

    def paintEvent(self, _event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        for (kind, text), (x, y) in zip(self._ITEMS, self._layout(self.width())):
            mid = y + self._ROW / 2
            painter.setBrush(Qt.BrushStyle.NoBrush)
            if kind == "gene":
                painter.setPen(QPen(_NODE_BORDER, 1.5))
                painter.drawRoundedRect(QRectF(x + 1, mid - 7, 20, 14), 4, 4)
            elif kind == "pips":
                for c, (letter, _name, colour, _step) in enumerate(_CHANNELS):
                    _draw_pip(painter, QRectF(x + c * 16 + 1, mid - 7, 14, 14), letter, colour, 1.0, self._font_pip)
            elif kind == "drug":
                painter.setPen(QPen(_ACCENT, 1.5))
                painter.drawRoundedRect(QRectF(x + 1, mid - 6, 20, 12), 6, 6)
            elif kind == "ppi":
                painter.setPen(QPen(QColor("#b9bcbc"), 1.5))
                painter.drawLine(QPointF(x + 1, mid), QPointF(x + 21, mid))
            else:
                pen = QPen(_ACCENT, 2.5)
                pen.setDashPattern([2, 1.6])
                painter.setPen(pen)
                painter.drawLine(QPointF(x + 1, mid), QPointF(x + 21, mid))
            painter.setPen(_MUTED)
            painter.setFont(self._font)
            painter.drawText(
                QRectF(x + self._SWATCH[kind] + 6, y, self.width(), self._ROW), Qt.AlignmentFlag.AlignVCenter, text
            )
        painter.end()


class _Console(QFrame):
    """Dark console of the backend's log lines; follows the newest line."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("OlConsole")
        self.setStyleSheet(f"#OlConsole {{ background: {CONSOLE_BACKGROUND}; border: none; border-radius: 0px; }}")
        self.setMinimumWidth(300)
        self._rows: list[QWidget] = []
        self._fades: list[QPropertyAnimation] = []
        # Whether the view is pinned to the newest line (off once the user scrolls up).
        self._follow = True

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        head = QFrame()
        head.setObjectName("OlConsoleHead")
        head.setFixedHeight(44)
        head.setStyleSheet(
            f"#OlConsoleHead {{ background: {CONSOLE_HEADER}; border: none; border-radius: 0px;"
            f" border-bottom: 1px solid {CONSOLE_BORDER}; }}"
        )
        head_layout = QHBoxLayout(head)
        head_layout.setContentsMargins(16, 0, 16, 0)
        head_layout.setSpacing(8)
        mono = f"font-family: {_MONO}, monospace; font-size: 12px; color: {CONSOLE_TEXT};"
        self._dot = QLabel()
        self._dot.setFixedSize(8, 8)
        self._status = QLabel()
        self._status.setStyleSheet(label_style(mono + " font-weight: 700;"))
        title = QLabel("console")
        title.setStyleSheet(label_style(mono))
        self._count = QLabel()
        self._count.setStyleSheet(label_style(mono))
        head_layout.addWidget(self._dot)
        head_layout.addWidget(self._status)
        head_layout.addSpacing(4)
        head_layout.addWidget(title)
        head_layout.addStretch(1)
        head_layout.addWidget(self._count)
        layout.addWidget(head)

        body = QWidget()
        body.setObjectName("OlConsoleBody")
        body.setStyleSheet(f"#OlConsoleBody {{ background: {CONSOLE_BACKGROUND}; }}")
        self._lines = QVBoxLayout(body)
        self._lines.setContentsMargins(16, 14, 12, 14)
        self._lines.setSpacing(6)
        self._caret = QFrame()
        self._caret.setFixedSize(8, 15)
        self._caret.setStyleSheet(f"background: {CONSOLE_TEXT}; border: none; border-radius: 0px;")
        self._lines.addWidget(self._caret)
        self._lines.addStretch(1)

        scroll = QScrollArea()
        scroll.setObjectName("OlConsoleScroll")
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setStyleSheet(
            f"#OlConsoleScroll {{ background: {CONSOLE_BACKGROUND}; border: none; border-radius: 0px; }}"
            "QScrollBar::handle:vertical { background: rgba(255, 255, 255, 0.22); }"
        )
        scroll.setWidget(body)
        layout.addWidget(scroll, 1)

        self._bar = scroll.verticalScrollBar()
        self._glide = QPropertyAnimation(self._bar, b"value", self)
        self._glide.setDuration(260)
        self._glide.setEasingCurve(QEasingCurve.Type.OutCubic)
        self._bar.rangeChanged.connect(self._on_range_changed)
        self._bar.valueChanged.connect(self._on_scrolled)

        self._blink = QTimer(self)
        self._blink.setInterval(500)
        self._blink.timeout.connect(self._toggle_caret)
        self.set_status("idle")

    def _toggle_caret(self) -> None:
        on = not self._caret.property("on")
        self._caret.setProperty("on", on)
        self._caret.setStyleSheet(
            f"background: {CONSOLE_TEXT if on else 'transparent'}; border: none; border-radius: 0px;"
        )

    def _on_range_changed(self, _minimum: int, maximum: int) -> None:
        if self._follow and maximum > self._bar.value():
            self._glide.stop()
            self._glide.setStartValue(self._bar.value())
            self._glide.setEndValue(maximum)
            self._glide.start()

    def _on_scrolled(self, value: int) -> None:
        if self._glide.state() != QPropertyAnimation.State.Running:
            self._follow = value >= self._bar.maximum() - 2

    def set_status(self, status: str) -> None:
        """`status` is "idle", "live", "done" or "failed"."""
        colour = {"live": "#4ade80", "failed": "#f87171"}.get(status, "#a1a1aa")
        self._dot.setStyleSheet(f"background: {colour}; border: none; border-radius: 4px;")
        self._status.setText(status.upper())
        live = status == "live"
        self._caret.setVisible(live)
        if live:
            self._blink.start()
        else:
            self._blink.stop()

    def clear(self) -> None:
        for row in self._rows:
            self._lines.removeWidget(row)
            row.deleteLater()
        self._rows.clear()
        self._fades.clear()
        self._follow = True
        self._count.setText("0 lines")

    def add_line(self, stamp: str, tag: str, text: str) -> None:
        row = QWidget()
        row_layout = QHBoxLayout(row)
        row_layout.setContentsMargins(0, 0, 0, 0)
        row_layout.setSpacing(10)
        mono = f"font-family: {_MONO}, monospace; font-size: 12px;"
        top = Qt.AlignmentFlag.AlignTop
        stamp_label = QLabel(stamp)
        stamp_label.setStyleSheet(label_style(f"{mono} color: {CONSOLE_TEXT_MUTED};"))
        tag_label = QLabel(tag)
        tag_label.setFixedWidth(72)
        tag_label.setStyleSheet(label_style(f"{mono} font-weight: 700; color: {_TAG_COLOURS.get(tag, '#e4e4e7')};"))
        text_label = QLabel(text)
        text_label.setWordWrap(True)
        text_label.setTextFormat(Qt.TextFormat.PlainText)
        text_label.setStyleSheet(label_style(f"{mono} color: {CONSOLE_TEXT};"))
        row_layout.addWidget(stamp_label, 0, top)
        row_layout.addWidget(tag_label, 0, top)
        row_layout.addWidget(text_label, 1, top)

        self._lines.insertWidget(len(self._rows), row)
        self._rows.append(row)
        self._count.setText(f"{len(self._rows)} line{'' if len(self._rows) == 1 else 's'}")

        # Fade the new line in, then drop the effect so old lines cost nothing.
        effect = QGraphicsOpacityEffect(row)
        row.setGraphicsEffect(effect)
        fade = QPropertyAnimation(effect, b"opacity", row)
        fade.setDuration(280)
        fade.setStartValue(0.0)
        fade.setEndValue(1.0)
        fade.setEasingCurve(QEasingCurve.Type.OutCubic)

        def settle() -> None:
            row.setGraphicsEffect(None)
            if fade in self._fades:
                self._fades.remove(fade)

        fade.finished.connect(settle)
        self._fades.append(fade)
        fade.start()


class _StageText(QWidget):
    """The current stage: step label, heading and one line of description.

    A stage change cross-fades: the old text lifts away as the new rises in.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedHeight(76)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._current: tuple[str, str, str] = ("", "", "")
        self._previous: tuple[str, str, str] | None = None
        self._blend = 1.0
        self._font_label = _font(_SANS, 12, QFont.Weight.Bold, 0.6)
        self._font_heading = _font(_SANS, 24, QFont.Weight.DemiBold, -0.4)
        self._font_text = _font(_SANS, 14)

    def set_stage(self, label: str, heading: str, text: str, snap: bool = False) -> None:
        stage = (label, heading, text)
        if stage == self._current:
            return
        self._previous = None if snap else self._current
        self._current = stage
        self._blend = 1.0 if snap else 0.0
        self.update()

    def advance(self, dt: float) -> None:
        if self._blend < 0.999:
            self._blend = min(1.0, _ease(self._blend, 1.0, dt, 0.16) + 0.001)
            self.update()

    def paintEvent(self, _event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        if self._previous is not None and self._blend < 0.999:
            self._paint_stage(painter, self._previous, 1 - self._blend, -10 * self._blend)
        self._paint_stage(painter, self._current, self._blend, 10 * (1 - self._blend))
        painter.end()

    def _paint_stage(self, painter: QPainter, stage: tuple[str, str, str], opacity: float, shift: float) -> None:
        label, heading, text = stage
        width = self.width()
        left = Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        painter.setOpacity(max(0.0, min(1.0, opacity)))
        painter.setPen(_ACCENT)
        painter.setFont(self._font_label)
        painter.drawText(QRectF(0, shift, width, 16), left, label.upper())
        painter.setPen(_INK)
        painter.setFont(self._font_heading)
        painter.drawText(QRectF(0, shift + 20, width, 32), left, heading)
        painter.setPen(_MUTED)
        painter.setFont(self._font_text)
        elided = QFontMetricsF(self._font_text).elidedText(text, Qt.TextElideMode.ElideRight, width)
        painter.drawText(QRectF(0, shift + 55, width, 20), left, elided)


class _Segments(QWidget):
    """Five-part progress bar: Graph, Omics, Drugs, Model, Predict."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedHeight(30)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._step = -1
        self._fill = [0.0] * len(_SEGMENTS)
        self._font = _font(_SANS, 12, QFont.Weight.Bold, 0.5)

    def set_step(self, step: int, snap: bool = False) -> None:
        self._step = step
        if snap:
            self._fill = self._targets()
            self.update()

    def _targets(self) -> list[float]:
        step = self._step
        return [
            1.0 if step > last else 0.0 if step < first else (step - first + 1) / (last - first + 1)
            for _label, first, last in _SEGMENTS
        ]

    def advance(self, dt: float) -> None:
        moved = False
        for i, target in enumerate(self._targets()):
            value = _ease(self._fill[i], target, dt, 0.3)
            if abs(value - self._fill[i]) > 0.0005:
                self._fill[i] = value
                moved = True
        if moved:
            self.update()

    def paintEvent(self, _event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        gap = 6
        cell = (self.width() - gap * (len(_SEGMENTS) - 1)) / len(_SEGMENTS)
        painter.setFont(self._font)
        for i, (label, first, last) in enumerate(_SEGMENTS):
            x = i * (cell + gap)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor(DIVIDER))
            painter.drawRoundedRect(QRectF(x, 0, cell, 6), 3, 3)
            if self._fill[i] > 0.002:
                painter.setBrush(_ACCENT)
                painter.drawRoundedRect(QRectF(x, 0, max(6.0, cell * self._fill[i]), 6), 3, 3)
            painter.setPen(_INK if first <= self._step <= last else _QUIET)
            painter.drawText(QRectF(x, 12, cell, 18), Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, label.upper())
        painter.end()


class OmicsLoader(QFrame):
    """The whole loader card. See the module docstring for how to drive it."""

    # Emitted once per run, after the "Done" stage has been on screen briefly.
    finished = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("OlCard")
        self.setStyleSheet(
            f"#OlCard {{ background: {SURFACE}; border: 1px solid {BORDER}; border-radius: {CARD_RADIUS}px; }}"
        )

        self._step = -1
        self._reached = 0       # furthest step the backend's log has confirmed
        self._in_step = 0.0     # seconds the current step has been shown
        self._running = False
        self._started_at = time.perf_counter()
        self._last_tick = self._started_at
        self._pending: list[tuple[str, str, str]] = []
        self._line_wait = 0.0

        root = QVBoxLayout(self)
        root.setContentsMargins(1, 1, 1, 1)  # inside the card's 1px border
        root.setSpacing(0)

        # -- inputs --
        files = QWidget()
        files_layout = QVBoxLayout(files)
        files_layout.setContentsMargins(24, 16, 24, 16)
        files_layout.setSpacing(12)
        head = QHBoxLayout()
        eyebrow = QLabel("INPUTS FOR THIS RUN")
        eyebrow.setStyleSheet(LABEL_CAPS_STYLE)
        self._loaded = QLabel()
        self._loaded.setStyleSheet(label_style(f"font-family: {_MONO}, monospace; font-size: 12px; color: {TEXT_SOFT};"))
        head.addWidget(eyebrow)
        head.addStretch(1)
        head.addWidget(self._loaded)
        files_layout.addLayout(head)
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(8)
        self._cards: list[_FileCard] = []
        for column, (role, _key, _stage, glyph, help_key) in enumerate(_FILES):
            card = _FileCard(role, glyph)
            card.setToolTip(HELP[help_key])
            grid.addWidget(card, 0, column)
            self._cards.append(card)
        files_layout.addLayout(grid)
        root.addWidget(files)
        root.addWidget(self._rule())

        # -- graph + console --
        body = QHBoxLayout()
        body.setContentsMargins(0, 0, 0, 0)
        body.setSpacing(0)
        graph_column = QVBoxLayout()
        graph_column.setContentsMargins(20, 12, 20, 14)
        graph_column.setSpacing(10)
        self._graph = _Constellation()
        self._graph.selected.connect(self._select)
        graph_column.addWidget(self._graph, 1)
        graph_column.addWidget(self._rule("#eeeeee"))
        graph_column.addWidget(_Legend())
        detail = QFrame()
        detail.setObjectName("OlDetail")
        detail.setStyleSheet(f"#OlDetail {{ background: {SURFACE_LOW}; border: none; border-radius: 8px; }}")
        detail.setFixedHeight(52)
        detail_layout = QHBoxLayout(detail)
        detail_layout.setContentsMargins(12, 6, 12, 6)
        self._detail = QLabel()
        self._detail.setWordWrap(True)
        self._detail.setTextFormat(Qt.TextFormat.RichText)
        self._detail.setStyleSheet(label_style("font-size: 13px;"))
        detail_layout.addWidget(self._detail)
        graph_column.addWidget(detail)
        body.addLayout(graph_column, 3)
        self._console = _Console()
        body.addWidget(self._console, 2)
        root.addLayout(body, 1)
        root.addWidget(self._rule())

        # -- stage + progress --
        foot = QVBoxLayout()
        foot.setContentsMargins(24, 16, 24, 18)
        foot.setSpacing(14)
        self._stage = _StageText()
        self._segments = _Segments()
        foot.addWidget(self._stage)
        foot.addWidget(self._segments)
        root.addLayout(foot)

        self._timer = QTimer(self)
        self._timer.setTimerType(Qt.TimerType.PreciseTimer)
        self._timer.setInterval(16)
        self._timer.timeout.connect(self._tick)

        self._selection: str | None = None
        self._cell_line = self._cell_name = ""
        self._stage.set_stage("Idle", "No run yet", "Start a run from Upload.", snap=True)
        self._render_files(snap=True)
        self._render_detail()

    @staticmethod
    def _rule(colour: str = DIVIDER) -> QFrame:
        rule = QFrame()
        rule.setFixedHeight(1)
        rule.setStyleSheet(f"background: {colour}; border: none; border-radius: 0px;")
        return rule

    # -- driving --------------------------------------------------------

    def start(self, files: dict[str, str], cell_line: str) -> None:
        """Reset for a new run. `files` maps each `_FILES` key to the name shown on its card."""
        for card, (_role, key, _stage, _glyph, _help) in zip(self._cards, _FILES):
            card.set_name(files.get(key, ""))
        self._graph.reset(cell_line)
        self._console.clear()
        self._console.set_status("live")
        self._pending.clear()
        self._selection = None
        self._cell_line, self._cell_name = cell_line, ""
        self._reached = 0
        self._running = True
        self._started_at = self._last_tick = time.perf_counter()
        self._step = -1
        self._enter(0, snap=True)
        self._timer.start()

    def log(self, line: str) -> None:
        """Show one real backend log line, and note how far the run has got."""
        match = _LOG_TAG.match(line)
        tag, text = match.groups() if match else ("predict" if line.startswith("Running inference") else "run", line)
        self._pending.append((self._stamp(), tag, text))

        if tag == "setup":
            self._reach(1)
        elif tag == "graph" and "loading" not in text:
            # The graph carries the omics features, so loading it covers the omics steps too.
            self._reach(5)
            drugs = re.search(r"(\d+) drugs", text)
            if drugs:
                self._graph.set_drug_total(int(drugs.group(1)))
        elif tag == "data":
            self._reach(6)
        elif tag in ("model", "epoch", "train", "eval"):
            self._reach(7)
        elif tag in ("cell-line", "predict"):
            self._reach(8)
            # "[cell-line] SIDM00136 -> A549 (Lung, NSCLC)"
            named = re.match(r"\S+ -> (.+?) \(", text) if tag == "cell-line" else None
            if named and named.group(1) != "unknown":
                self._cell_name = named.group(1)
                self._graph.set_cell_name(self._cell_name)

    def set_progress(self, fraction: float) -> None:
        """Backend-reported progress, 0..1: keeps the steps moving for a backend whose log lines aren't recognised."""
        self._reach(min(8, int(max(0.0, fraction) * 8)))

    def complete(self) -> None:
        """The run finished: play out the remaining steps, then emit `finished`."""
        if self._running:
            self._pending.append((self._stamp(), "done", "results ready -> opening Results"))
            self._reach(_LAST_STEP)

    def fail(self, message: str) -> None:
        if not self._running:
            return
        self._running = False
        self._pending.append((self._stamp(), "error", message))
        self._console.set_status("failed")
        self._stage.set_stage("Stopped", "Run failed", message)
        for card, (_role, _key, stage, _glyph, _help) in zip(self._cards, _FILES):
            if stage == self._step:
                card.set_state("waiting")

    def _reach(self, step: int) -> None:
        self._reached = max(self._reached, step)

    def _stamp(self) -> str:
        minutes, seconds = divmod(time.perf_counter() - self._started_at, 60)
        return f"{int(minutes)}:{seconds:04.1f}"

    def _enter(self, step: int, snap: bool = False) -> None:
        self._step = step
        self._in_step = 0.0
        self._graph.set_step(step)
        self._segments.set_step(step, snap)
        segment = next(i for i, (_label, first, last) in enumerate(_SEGMENTS) if first <= step <= last)
        label = "Complete" if step == _LAST_STEP else f"Step {segment + 1} of {len(_SEGMENTS)} · {_SEGMENTS[segment][0]}"
        self._stage.set_stage(label, *_STAGES[step], snap=snap)
        self._render_files(snap)
        self._render_detail()
        if step == _LAST_STEP:
            self._console.set_status("done")

    def _tick(self) -> None:
        now = time.perf_counter()
        dt = min(0.1, now - self._last_tick)
        self._last_tick = now

        if self._running:
            self._in_step += dt
            if self._step < self._reached:
                hurry = _CATCH_UP if self._reached == _LAST_STEP else 1.0
                if self._in_step >= _MIN_DWELL[self._step] * hurry:
                    self._enter(self._step + 1)
            elif self._step == _LAST_STEP and self._in_step >= _DONE_HOLD:
                self._running = False
                if not self.isVisible():
                    self._timer.stop()
                self.finished.emit()

        # A poll delivers its lines in one burst; let them in one at a time.
        self._line_wait -= dt
        if self._pending and self._line_wait <= 0:
            self._console.add_line(*self._pending.pop(0))
            self._line_wait = 0.07

        self._graph.advance(dt)
        self._stage.advance(dt)
        self._segments.advance(dt)
        for card in self._cards:
            card.advance(dt)

    # The clock runs while the card is on screen, and also off screen for as
    # long as a run is being followed, so `finished` still fires if the user
    # has clicked away to another page.
    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._last_tick = time.perf_counter()
        self._timer.start()

    def hideEvent(self, event) -> None:  # noqa: N802
        super().hideEvent(event)
        if not self._running:
            self._timer.stop()

    # -- rendering ------------------------------------------------------

    def _render_files(self, snap: bool = False) -> None:
        loaded = 0
        for card, (_role, _key, stage, _glyph, _help) in zip(self._cards, _FILES):
            state = "done" if self._step > stage else "active" if self._step == stage else "waiting"
            loaded += state == "done"
            card.set_state(state, snap)
        self._loaded.setText(f"{loaded} of {len(_FILES)} loaded")

    def _select(self, node_id: str) -> None:
        self._selection = node_id
        self._graph.set_selection(node_id)
        self._render_detail()

    def _render_detail(self) -> None:
        title, text = self._describe(self._selection)
        self._detail.setText(
            f'<span style="font-family: {_MONO}; font-weight: 700; color: {TEXT};">{html.escape(title)}</span>'
            f'&nbsp;&nbsp;<span style="color: {TEXT_MUTED};">{html.escape(text)}</span>'
        )

    def _describe(self, node_id: str | None) -> tuple[str, str]:
        if node_id is None:
            return "Click a shape", "Each one stands for a part of this run."
        if node_id == "sample":
            title = f"{self._cell_line} · {self._cell_name}" if self._cell_name else self._cell_line
            return (
                title or "Sample",
                "The cell line being scored. Its three omics profiles are joined into one feature vector on its graph node.",
            )
        if node_id in _GENES:
            neighbours = [b if a == node_id else a for a, b in _PPI if node_id in (a, b)]
            drugs = [name for name, targets, _angle in _DRUGS if node_id in targets]
            targeted = f"Targeted by {', '.join(drugs)}." if drugs else "No targeting drug shown."
            return f"{node_id} · gene", f"Interacts with {', '.join(neighbours)}. {targeted}"
        targets = next(targets for name, targets, _angle in _DRUGS if name == node_id)
        state = "Response predicted." if self._step >= 8 else "Waiting to be scored."
        return f"{node_id} · drug", f"In GDSC2. Targets {', '.join(targets)}. Structure as a Morgan fingerprint. {state}"
