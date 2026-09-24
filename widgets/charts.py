"""Custom-painted, hoverable charts for the Model Analytics page.

Both charts are fed real per-run data, never mocked:

- `TrainingCurveWidget`: how well the model fits, from `on_training_history`
  (backend/model_backends/hetero_gnn.py).
- `PredictionScatterWidget`: predicted IC50 vs. confidence, one dot per
  ranked drug, from the same rows as the results table.

They share `_HoverChart`, which gives every chart the same hover behaviour:
the pointed-at mark is highlighted (the rest dim, dashed guides run to the
axes) and a tooltip names it, without waiting for Qt's tooltip delay.
Pointing at a mark is deliberate, so the chart responds at once.

Visual language: plain plot area on the card's own surface, recessive solid
gridlines, tick labels in muted text ink, and series colour only on marks,
never on text.
"""

from __future__ import annotations

import math

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import QWidget

from styles.theme import SURFACE, SURFACE_CONTAINER, TEXT_MUTED, TEXT_SOFT
from widgets.formatting import format_ic50
from widgets.tooltip import hide_tooltip, show_tooltip

# Series colours: a warm/cool pair, so the two training metrics never read
# as the same series.
TRAIN_LOSS_COLOR = "#8a3419"
VAL_RMSE_COLOR = "#1c4a7a"
VALIDATION_POINT_COLOR = "#1c4a7a"
# Prediction confidence runs amber (uncertain) -> green (trustworthy).
LOW_CONFIDENCE_COLOR = QColor(133, 77, 24)
HIGH_CONFIDENCE_COLOR = QColor(27, 94, 54)

_AXIS_FONT_PX = 11
_WAITING_TEXT_COLOR = TEXT_MUTED


def _nice_ticks(lo: float, hi: float, target: int = 4) -> list[float]:
    """Round-numbered ticks covering [lo, hi] -- steps of 1, 2 or 5 x 10^k."""
    span = hi - lo
    if span <= 0:
        return [lo]
    raw_step = span / target
    magnitude = 10 ** math.floor(math.log10(raw_step))
    step = next(m * magnitude for m in (1, 2, 5, 10) if m * magnitude >= raw_step)
    first = math.ceil(lo / step) * step
    ticks = []
    value = first
    while value <= hi + step * 1e-9:
        ticks.append(round(value, 10))
        value += step
    return ticks


def _tick_text(value: float) -> str:
    return f"{value:.0f}" if float(value).is_integer() else f"{value:g}"


class _HoverChart(QWidget):
    """Base for charts whose marks can be pointed at.

    Subclasses record where they painted each hoverable mark in
    `self._hit_points` during `paintEvent`, and implement `_tooltip_for`.
    `_hit_test` defaults to "nearest mark within `_HOVER_RADIUS`", which is
    what a scatter needs; line charts override it to snap to the nearest x.
    """

    # Generous: readers aim near a dot, not dead centre on it.
    _HOVER_RADIUS = 12.0
    _HOVER_CURSOR = Qt.CursorShape.PointingHandCursor
    # Alpha of the other marks while one is hovered, so the picked one stands out.
    _DIMMED_ALPHA = 55

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMouseTracking(True)
        self._hit_points: list[QPointF] = []
        self._hovered_index: int | None = None

    # -- hover plumbing ---------------------------------------------------

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        self._set_hovered(self._hit_test(event.position()))
        super().mouseMoveEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802
        self._set_hovered(None)
        super().leaveEvent(event)

    def _hit_test(self, position: QPointF) -> int | None:
        nearest_index, nearest_distance = None, self._HOVER_RADIUS
        for index, centre in enumerate(self._hit_points):
            distance = math.hypot(centre.x() - position.x(), centre.y() - position.y())
            if distance <= nearest_distance:
                nearest_index, nearest_distance = index, distance
        return nearest_index

    def _tooltip_for(self, index: int) -> tuple[str, str, QPointF]:
        """(title, body, anchor point in widget coordinates) for mark `index`."""
        raise NotImplementedError

    def _set_hovered(self, index: int | None) -> None:
        if index == self._hovered_index:
            return
        self._hovered_index = index
        self.update()
        if index is None:
            self.unsetCursor()
            hide_tooltip(self)
            return
        self.setCursor(self._HOVER_CURSOR)
        title, body, anchor = self._tooltip_for(index)
        show_tooltip(body, self.mapToGlobal(anchor.toPoint()), self, title=title, above=True)

    def _reset_hover(self) -> None:
        """Forget the hovered mark when the data underneath it changes."""
        self._hit_points = []
        self._hovered_index = None
        self.unsetCursor()
        hide_tooltip(self)

    # -- shared painting ----------------------------------------------------

    def _begin_paint(self) -> QPainter:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor(SURFACE))
        font = QFont(painter.font())
        font.setPixelSize(_AXIS_FONT_PX)
        painter.setFont(font)
        return painter

    def _paint_waiting(self, painter: QPainter, text: str) -> None:
        painter.setPen(QColor(_WAITING_TEXT_COLOR))
        painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, text)

    @staticmethod
    def _grid_pen() -> QPen:
        return QPen(QColor(SURFACE_CONTAINER), 1)

    @staticmethod
    def _guide_pen() -> QPen:
        pen = QPen(QColor(TEXT_SOFT), 1)
        pen.setStyle(Qt.PenStyle.DashLine)
        return pen

    @staticmethod
    def _paint_highlight(painter: QPainter, centre: QPointF, color: QColor, radius: float = 5.5) -> None:
        """The picked mark: a soft halo in its own colour, then a white-ringed dot."""
        solid = QColor(color)
        solid.setAlpha(255)
        halo = QColor(color)
        halo.setAlpha(55)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(halo)
        painter.drawEllipse(centre, radius + 5, radius + 5)
        painter.setPen(QPen(QColor(SURFACE), 2))
        painter.setBrush(solid)
        painter.drawEllipse(centre, radius, radius)

    @staticmethod
    def _paint_dot(painter: QPainter, centre: QPointF, color: QColor, radius: float = 4.0) -> None:
        """A standard mark: filled dot with a thin surface ring, so overlaps stay legible."""
        painter.setPen(QPen(QColor(SURFACE), 1.2))
        painter.setBrush(color)
        painter.drawEllipse(centre, radius, radius)

    @staticmethod
    def _draw_right_aligned(painter: QPainter, right_x: float, baseline_y: float, text: str) -> None:
        width = painter.fontMetrics().horizontalAdvance(text)
        painter.drawText(QPointF(right_x - width, baseline_y), text)

    @staticmethod
    def _draw_centered(painter: QPainter, centre_x: float, baseline_y: float, text: str) -> None:
        width = painter.fontMetrics().horizontalAdvance(text)
        painter.drawText(QPointF(centre_x - width / 2, baseline_y), text)


class PredictionScatterWidget(_HoverChart):
    """Predicted IC50 (log x) vs. confidence (y), one dot per ranked drug.

    Each dot is coloured along a low-to-high confidence gradient, so the
    colour itself carries information rather than being decorative. Clicking
    a dot emits `drugClicked` with that drug's id.
    """

    drugClicked = Signal(str)

    _MARGINS = (46, 14, 16, 42)  # left, top, right, bottom

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(340)
        # (ic50, confidence, drug name, drug id)
        self._points: list[tuple[float, float, str, str]] = []

    def set_points(self, points: list[tuple[float, float, str, str]]) -> None:
        self._points = points
        self._reset_hover()
        self.update()

    def mousePressEvent(self, event) -> None:  # noqa: N802
        index = self._hit_test(event.position())
        if event.button() == Qt.MouseButton.LeftButton and index is not None:
            self.drugClicked.emit(self._points[index][3])
            return
        super().mousePressEvent(event)

    @staticmethod
    def confidence_color(confidence_percent: float) -> QColor:
        t = max(0.0, min(1.0, confidence_percent / 100.0))
        low, high = LOW_CONFIDENCE_COLOR, HIGH_CONFIDENCE_COLOR
        return QColor(
            int(low.red() + (high.red() - low.red()) * t),
            int(low.green() + (high.green() - low.green()) * t),
            int(low.blue() + (high.blue() - low.blue()) * t),
            210,
        )

    def _tooltip_for(self, index: int) -> tuple[str, str, QPointF]:
        ic50, confidence, name, _drug_id = self._points[index]
        centre = self._hit_points[index]
        return (
            name,
            f"IC50 {format_ic50(ic50)} µM · {confidence:.0f}% confidence",
            QPointF(centre.x(), centre.y() - 8),
        )

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = self._begin_paint()
        if not self._points:
            self._paint_waiting(painter, "Waiting for results...")
            return

        left, top, right, bottom = self._MARGINS
        plot = QRectF(left, top, max(1, self.width() - left - right), max(1, self.height() - top - bottom))

        log_ic50 = [math.log10(max(ic50, 1e-6)) for ic50, _c, _n, _id in self._points]
        lo, hi = min(log_ic50), max(log_ic50)
        if hi - lo < 1e-9:
            lo, hi = lo - 0.5, hi + 0.5
        pad = (hi - lo) * 0.04
        lo, hi = lo - pad, hi + pad

        def map_x(log_value: float) -> float:
            return plot.left() + (log_value - lo) / (hi - lo) * plot.width()

        def map_y(confidence: float) -> float:
            return plot.bottom() - confidence / 100.0 * plot.height()

        # Gridlines + y ticks (confidence %).
        for tick in (0, 25, 50, 75, 100):
            y = map_y(tick)
            painter.setPen(self._grid_pen())
            painter.drawLine(QPointF(plot.left(), y), QPointF(plot.right(), y))
            painter.setPen(QColor(TEXT_SOFT))
            self._draw_right_aligned(painter, plot.left() - 8, y + 4, f"{tick}%")

        # X ticks at whole decades of IC50 (log axis).
        painter.setPen(QColor(TEXT_SOFT))
        decades = range(math.ceil(lo), math.floor(hi) + 1)
        for decade in decades:
            x = map_x(decade)
            self._draw_centered(painter, x, plot.bottom() + 16, format_ic50(10.0**decade))
        self._draw_centered(painter, plot.center().x(), self.height() - 6, "Predicted IC50 (µM, log scale)")

        hovered = self._hovered_index
        self._hit_points = []
        for index, ((_ic50, confidence, _name, _id), log_value) in enumerate(zip(self._points, log_ic50)):
            centre = QPointF(map_x(log_value), map_y(confidence))
            self._hit_points.append(centre)
            color = self.confidence_color(confidence)
            if hovered is not None and index != hovered:
                color.setAlpha(self._DIMMED_ALPHA)
            self._paint_dot(painter, centre, color)

        if hovered is not None and hovered < len(self._hit_points):
            centre = self._hit_points[hovered]
            painter.setPen(self._guide_pen())
            painter.drawLine(QPointF(plot.left(), centre.y()), centre)
            painter.drawLine(centre, QPointF(centre.x(), plot.bottom()))
            self._paint_highlight(painter, centre, self.confidence_color(self._points[hovered][1]))


class TrainingCurveWidget(_HoverChart):
    """Model performance for this run, in one of two real-data modes.

    - **Curve** (the backend actually trained): train loss and validation
      RMSE per epoch, drawn as two stacked panels sharing the epoch axis.
      Each metric keeps its own y-scale rather than being squashed onto one
      shared axis, since loss and RMSE are different units. Hovering snaps a
      crosshair to the nearest epoch across both panels.
    - **Checkpoint** (the backend loaded pretrained weights, so there are no
      epochs): held-out validation pairs, actual vs. predicted ln(IC50), with
      a dashed perfect-prediction diagonal. Hovering picks the nearest point.

    `mode()` says which applies, based on the fields present on the first
    point (see `backend/schemas/results.py::TrainingHistoryPoint`).
    """

    _CURVE_MARGINS = (54, 12, 16, 40)  # left, top, right, bottom
    _PANEL_GAP = 22
    _SCATTER_MARGINS = (46, 14, 16, 42)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(340)
        self._history: list[dict] = []

    def set_history(self, history: list[dict]) -> None:
        self._history = history
        self._reset_hover()
        self.update()

    def mode(self) -> str:
        """"curve", "scatter" (checkpoint validation pairs) or "empty"."""
        if not self._history:
            return "empty"
        first = self._history[0]
        if first.get("epoch") is not None:
            return "curve" if len(self._history) >= 2 else "empty"
        if first.get("actual_ln_ic50") is not None:
            return "scatter"
        return "empty"

    # -- hover --------------------------------------------------------------

    def _hit_test(self, position: QPointF) -> int | None:
        if self.mode() != "curve":
            return super()._hit_test(position)
        # Crosshair: snap to the nearest epoch by x alone, anywhere in the plot.
        # `_hit_points` holds the train-loss point per epoch (x is shared).
        left, top, right, bottom = self._CURVE_MARGINS
        if not self._hit_points or not (left - 10 <= position.x() <= self.width() - right + 10):
            return None
        if not (top <= position.y() <= self.height() - bottom):
            return None
        return min(range(len(self._hit_points)), key=lambda i: abs(self._hit_points[i].x() - position.x()))

    def _tooltip_for(self, index: int) -> tuple[str, str, QPointF]:
        point = self._history[index]
        centre = self._hit_points[index]
        if self.mode() == "curve":
            return (
                f"Epoch {point['epoch']}",
                f"Train loss {point['train_loss']:.3f} · Val RMSE {point['val_rmse']:.3f}",
                QPointF(centre.x(), centre.y() - 8),
            )
        actual, predicted = point["actual_ln_ic50"], point["predicted_ln_ic50"]
        return (
            "Validation sample",
            f"Actual {actual:.2f} · Predicted {predicted:.2f} ln(IC50)",
            QPointF(centre.x(), centre.y() - 8),
        )

    # -- painting -------------------------------------------------------------

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = self._begin_paint()
        mode = self.mode()
        if mode == "curve":
            self._paint_training_curve(painter)
        elif mode == "scatter":
            self._paint_validation_scatter(painter)
        else:
            self._paint_waiting(painter, "Waiting for training history...")

    def _paint_validation_scatter(self, painter: QPainter) -> None:
        actual = [point["actual_ln_ic50"] for point in self._history]
        predicted = [point["predicted_ln_ic50"] for point in self._history]
        ticks = _nice_ticks(min(min(actual), min(predicted)), max(max(actual), max(predicted)), target=6)
        lo = min(ticks[0], min(actual), min(predicted))
        hi = max(ticks[-1], max(actual), max(predicted))
        span = (hi - lo) or 1.0

        left, top, right, bottom = self._SCATTER_MARGINS
        plot = QRectF(left, top, max(1, self.width() - left - right), max(1, self.height() - top - bottom))

        def map_xy(x_value: float, y_value: float) -> QPointF:
            return QPointF(
                plot.left() + (x_value - lo) / span * plot.width(),
                plot.bottom() - (y_value - lo) / span * plot.height(),
            )

        for tick in ticks:
            y = map_xy(lo, tick).y()
            painter.setPen(self._grid_pen())
            painter.drawLine(QPointF(plot.left(), y), QPointF(plot.right(), y))
            painter.setPen(QColor(TEXT_SOFT))
            self._draw_right_aligned(painter, plot.left() - 8, y + 4, _tick_text(tick))
            self._draw_centered(painter, map_xy(tick, lo).x(), plot.bottom() + 16, _tick_text(tick))
        self._draw_centered(painter, plot.center().x(), self.height() - 6, "Actual ln(IC50)")
        painter.save()
        painter.translate(12, plot.center().y())
        painter.rotate(-90)
        self._draw_centered(painter, 0, 0, "Predicted ln(IC50)")
        painter.restore()

        painter.setPen(self._guide_pen())
        painter.drawLine(map_xy(lo, lo), map_xy(hi, hi))

        hovered = self._hovered_index
        self._hit_points = []
        for index, (x_value, y_value) in enumerate(zip(actual, predicted)):
            centre = map_xy(x_value, y_value)
            self._hit_points.append(centre)
            color = QColor(VALIDATION_POINT_COLOR)
            color.setAlpha(self._DIMMED_ALPHA if hovered is not None and index != hovered else 190)
            self._paint_dot(painter, centre, color, radius=3.5)

        if hovered is not None and hovered < len(self._hit_points):
            centre = self._hit_points[hovered]
            painter.setPen(self._guide_pen())
            painter.drawLine(QPointF(plot.left(), centre.y()), centre)
            painter.drawLine(centre, QPointF(centre.x(), plot.bottom()))
            self._paint_highlight(painter, centre, QColor(VALIDATION_POINT_COLOR))

    def _paint_training_curve(self, painter: QPainter) -> None:
        epochs = [point["epoch"] for point in self._history]
        e_lo, e_hi = min(epochs), max(epochs)
        e_span = (e_hi - e_lo) or 1

        left, top, right, bottom = self._CURVE_MARGINS
        width = max(1, self.width() - left - right)
        panel_height = max(1, (self.height() - top - bottom - self._PANEL_GAP) / 2)
        panels = [
            ("Train loss", "train_loss", TRAIN_LOSS_COLOR, QRectF(left, top, width, panel_height)),
            (
                "Validation RMSE",
                "val_rmse",
                VAL_RMSE_COLOR,
                QRectF(left, top + panel_height + self._PANEL_GAP, width, panel_height),
            ),
        ]

        def map_x(epoch: float) -> float:
            return left + (epoch - e_lo) / e_span * width

        hovered = self._hovered_index
        panel_points: list[list[QPointF]] = []
        for label, key, color, rect in panels:
            values = [point[key] for point in self._history]
            ticks = _nice_ticks(min(values), max(values), target=3)
            v_lo, v_hi = min(ticks[0], min(values)), max(ticks[-1], max(values))
            v_span = (v_hi - v_lo) or 1.0

            def map_y(value: float, rect=rect, v_lo=v_lo, v_span=v_span) -> float:
                return rect.bottom() - (value - v_lo) / v_span * rect.height()

            for tick in ticks:
                y = map_y(tick)
                painter.setPen(self._grid_pen())
                painter.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))
                painter.setPen(QColor(TEXT_SOFT))
                self._draw_right_aligned(painter, rect.left() - 8, y + 4, f"{tick:g}")

            points = [QPointF(map_x(epoch), map_y(value)) for epoch, value in zip(epochs, values)]
            path = QPainterPath(points[0])
            for point in points[1:]:
                path.lineTo(point)
            painter.setPen(QPen(QColor(color), 2))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPath(path)

            # Panel label in text ink, top-right: loss curves fall from the
            # top-left, so that corner is where the line starts. The card's
            # legend row carries the colour key.
            painter.setPen(QColor(TEXT_MUTED))
            self._draw_right_aligned(painter, rect.right(), rect.top() + 11, label)
            panel_points.append(points)

        # Shared epoch axis under the bottom panel.
        painter.setPen(QColor(TEXT_SOFT))
        bottom_rect = panels[-1][3]
        for tick in _nice_ticks(e_lo, e_hi, target=6):
            if e_lo <= tick <= e_hi:
                self._draw_centered(painter, map_x(tick), bottom_rect.bottom() + 16, _tick_text(tick))
        self._draw_centered(painter, left + width / 2, self.height() - 6, "Epoch")

        self._hit_points = panel_points[0]
        if hovered is not None and hovered < len(epochs):
            x = map_x(epochs[hovered])
            painter.setPen(self._guide_pen())
            painter.drawLine(QPointF(x, panels[0][3].top()), QPointF(x, bottom_rect.bottom()))
            for (_label, _key, color, _rect), points in zip(panels, panel_points):
                self._paint_highlight(painter, points[hovered], QColor(color), radius=4.5)
