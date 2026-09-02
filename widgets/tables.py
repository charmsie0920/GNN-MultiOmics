"""Shared data-table styling and cell widgets.

The pipeline-run table on the model execution log page and the predicted
drug results table on the final results page both present a list of records
with a status/ranking column and a progress-style measurement column. This
module gives both tables one shared visual vocabulary:

- `style_data_table` applies the same header, border, and row QSS to any
  `QTableWidget`.
- `build_status_badge` renders a state as a filled dot, a bordered pill, or
  plain muted text — used for pipeline stage ("QC PASSED" / "NORMALIZING" /
  "PENDING") and drug sensitivity ranking ("HIGH SENSITIVITY" / "MEDIUM" /
  "LOW") alike.
- `build_mini_progress_bar` renders a thin two-tone progress track, reused
  for pipeline completion, prediction confidence, and expression-level bars.

Qt/Fusion quirk: a plain `QWidget` placed via `QTableWidget.setCellWidget`
renders with an unwanted embossed border under the Fusion style unless it
explicitly opts out with its own transparent, borderless stylesheet. Rather
than repeat that workaround everywhere, `transparent_cell_widget()` builds
the base container every custom cell widget in this module starts from.
"""

from __future__ import annotations

from PySide6.QtCore import QTimer, Qt
from PySide6.QtGui import QColor, QPainter, QPen
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)

from styles.theme import BORDER, PRIMARY, SURFACE, SURFACE_CONTAINER, TEXT_FAINT, TEXT_MUTED, WINDOW_BACKGROUND, table_stylesheet
from widgets.icons import icon_text

# Visual language for the three states a badge can represent. "positive"
# reads as complete/favorable, "neutral" as in-progress/moderate, and
# "muted" as not-yet-relevant/unfavorable.
_BADGE_TONES = ("positive", "neutral", "muted")


class SpinningIconLabel(QLabel):
    """A small rotating ring, standing in for the static "sync" icon.

    Rotating the tiny "sync" text glyph itself (an earlier version of this
    widget) read as choppy and clipped at small font sizes — glyph hinting
    at 12px doesn't rotate cleanly. Drawing a plain arc directly avoids
    that entirely and reads as a standard loading spinner.
    """

    def __init__(self, text: str, color: str, *, interval_ms: int = 16, degrees_per_tick: float = 6.0) -> None:
        del text  # kept for call-site compatibility; nothing is drawn from it
        super().__init__()
        self._angle = 0.0
        self._degrees_per_tick = degrees_per_tick
        self._color = QColor(color)
        self.setFixedSize(14, 14)
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._advance)
        self._timer.start(interval_ms)

    def _advance(self) -> None:
        self._angle = (self._angle + self._degrees_per_tick) % 360
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802 (Qt override)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        pen = QPen(self._color)
        pen.setWidth(2)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        rect = self.rect().adjusted(1, 1, -1, -1)
        # Qt angles are in 1/16ths of a degree, measured counter-clockwise
        # from the 3 o'clock position.
        start_angle = int(-self._angle * 16)
        span_angle = int(270 * 16)
        painter.drawArc(rect, start_angle, span_angle)
        painter.end()


def transparent_cell_widget() -> QWidget:
    """Build a borderless, transparent container for use as a table cell widget.

    Returns:
        A `QWidget` pre-styled to avoid the Fusion style's default embossed
        panel look, so its own layout/children fully control appearance.
    """
    wrapper = QWidget()
    wrapper.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
    wrapper.setStyleSheet("background: transparent; border: none;")
    wrapper.setAutoFillBackground(False)
    return wrapper


def style_data_table(table: QTableWidget, *, header_background: str = WINDOW_BACKGROUND) -> None:
    """Apply the shared read-only, borderless data-table appearance and behavior.

    Configures the table as a non-editable, non-selectable list view (rows
    are informational, not interactive) and installs the shared QSS from
    `styles.theme.table_stylesheet`.

    Args:
        table: The table to configure, already sized (rows/columns) and
            populated by the caller.
        header_background: Background color for the header row; pass the
            hosting card's header strip color so the two blend together.
    """
    table.verticalHeader().setVisible(False)
    table.setShowGrid(False)
    table.setAlternatingRowColors(False)
    table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
    table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    table.setFocusPolicy(Qt.FocusPolicy.NoFocus)
    table.setWordWrap(False)
    table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
    table.horizontalHeader().setStretchLastSection(True)
    table.horizontalHeader().setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
    table.horizontalHeader().setHighlightSections(False)
    table.verticalHeader().setDefaultSectionSize(52)
    table.setStyleSheet(table_stylesheet(header_background))


def build_status_badge(text: str, tone: str, *, icon_name: str | None = None) -> QWidget:
    """Render a row's state as a dot, bordered pill, or plain muted label.

    Args:
        text: The status text to display (e.g. "QC PASSED", "HIGH SENSITIVITY").
        tone: One of `"positive"`, `"neutral"`, or `"muted"` — see module
            docstring for the meaning of each.
        icon_name: Optional leading icon (key into `widgets.icons.ICON_MAP`),
            only rendered for the `"neutral"` pill tone (e.g. a spinning
            "sync" glyph for an in-progress pipeline stage).

    Returns:
        A `QWidget` suitable for `QTableWidget.setCellWidget`.

    Raises:
        ValueError: If `tone` is not one of the supported values.
    """
    if tone not in _BADGE_TONES:
        raise ValueError(f"tone must be one of {_BADGE_TONES}, got {tone!r}")

    wrapper = transparent_cell_widget()
    row = QHBoxLayout(wrapper)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(0)

    # Every tone renders inside the same bordered pill container so all
    # three states read as one consistent "badge" shape; only the inner
    # content (dot/icon/text color) changes per tone.
    pill = QFrame()
    pill.setStyleSheet(f"background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 4px;")
    pill_layout = QHBoxLayout(pill)
    pill_layout.setContentsMargins(8, 3, 8, 3)
    pill_layout.setSpacing(6)

    if tone == "positive":
        dot = QFrame()
        dot.setFixedSize(8, 8)
        dot.setStyleSheet(f"background: {PRIMARY}; border-radius: 4px;")
        pill_layout.addWidget(dot)
        text_color = TEXT_MUTED
    elif tone == "neutral":
        if icon_name is not None:
            if icon_name == "sync":
                icon = SpinningIconLabel(icon_text(icon_name), PRIMARY)
            else:
                icon = QLabel(icon_text(icon_name))
            icon.setStyleSheet(f"background: transparent; border: none; color: {PRIMARY}; font-size: 12px;")
            pill_layout.addWidget(icon)
        text_color = PRIMARY
    else:  # muted
        text_color = TEXT_FAINT

    label = QLabel(text)
    label.setStyleSheet(
        f"background: transparent; border: none; font-size: 11px; font-weight: 700;"
        f" letter-spacing: 0.05em; color: {text_color};"
    )
    pill_layout.addWidget(label)

    row.addWidget(pill)
    row.addStretch(1)
    return wrapper


def build_mini_progress_bar(percent: int, *, height: int = 4) -> QWidget:
    """Render a thin two-tone progress track.

    Args:
        percent: Fill amount from 0 to 100.
        height: Track height in pixels.

    Returns:
        A `QWidget` suitable for `QTableWidget.setCellWidget`, or for
        embedding directly in any other layout that wants a compact
        progress indicator.
    """
    wrapper = transparent_cell_widget()
    outer_layout = QVBoxLayout(wrapper)
    outer_layout.setContentsMargins(16, 0, 16, 0)

    track = QFrame()
    track.setFixedHeight(height)
    track.setStyleSheet(f"background: {SURFACE_CONTAINER}; border: none;")
    track_layout = QHBoxLayout(track)
    track_layout.setContentsMargins(0, 0, 0, 0)
    track_layout.setSpacing(0)

    fill = QFrame()
    fill.setFixedHeight(height)
    fill.setStyleSheet(f"background: {PRIMARY}; border: none;")
    track_layout.addWidget(fill, percent)
    if percent <= 0:
        fill.setFixedWidth(0)
    track_layout.addStretch(max(0, 100 - percent))

    outer_layout.addWidget(track)
    return wrapper
