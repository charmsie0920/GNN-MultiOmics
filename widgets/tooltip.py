"""App-wide custom tooltip, replacing Qt's native one.

Qt's native tooltip can't be made to look consistent here: it inherits the
selector-less inline stylesheets of whatever widget it belongs to (most cards
set `background`/`border`/`border-radius` that way), and those outrank any
app-level `QToolTip` rule. This tooltip is a separate, self-painted window, so
it looks identical everywhere: a dark rounded card with a soft shadow and a
short fade-in.

`install_tooltips(app)` routes every tooltip in the app through it with no
per-widget changes -- plain `setToolTip`, table header tooltips and table item
tooltips all keep working as before. Widgets that want a tooltip on their own
terms (the results scatter's per-dot hover) call `show_tooltip` /
`hide_tooltip` directly.

Timing still comes from Qt (wake-up delay, see `styles.theme`); only the
display is replaced. Escape, a click, scrolling, or leaving the owning widget
dismisses it.
"""

from __future__ import annotations

import atexit

import shiboken6
from PySide6.QtCore import QEvent, QObject, QPoint, QPropertyAnimation, QRect, QRectF, Qt
from PySide6.QtGui import QColor, QGuiApplication, QPainter, QPen
from PySide6.QtWidgets import QAbstractItemView, QApplication, QHeaderView, QLabel, QVBoxLayout, QWidget

from styles.theme import TOOLTIP_BACKGROUND, TOOLTIP_TEXT, TOOLTIP_TEXT_MUTED, label_style


class _ToolTipWindow(QWidget):
    """The floating tooltip itself: an optional bold title over a line of text."""

    _RADIUS = 8
    # Transparent margin around the card that the drop shadow is painted into.
    _SHADOW = 12
    _MAX_TEXT_WIDTH = 320
    # Where the card sits relative to the cursor / anchor point.
    _CURSOR_OFFSET = QPoint(12, 18)
    _ANCHOR_GAP = 10
    _FADE_MS = 90

    def __init__(self) -> None:
        super().__init__(
            None,
            Qt.WindowType.ToolTip | Qt.WindowType.FramelessWindowHint | Qt.WindowType.NoDropShadowWindowHint,
        )
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.setAttribute(Qt.WidgetAttribute.WA_ShowWithoutActivating)
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)

        layout = QVBoxLayout(self)
        margin = self._SHADOW
        layout.setContentsMargins(margin + 12, margin + 8, margin + 12, margin + 9)
        layout.setSpacing(2)
        self._title = QLabel()
        self._title.setStyleSheet(label_style(f"color: {TOOLTIP_TEXT}; font-size: 12px; font-weight: 600;"))
        self._body = QLabel()
        layout.addWidget(self._title)
        layout.addWidget(self._body)

        self._fade = QPropertyAnimation(self, b"windowOpacity", self)
        self._fade.setDuration(self._FADE_MS)
        self._fade.setStartValue(0.0)
        self._fade.setEndValue(1.0)

        # The widget the tooltip belongs to, and the area of it (in its own
        # coordinates) the cursor must stay inside; None = the whole widget.
        self.owner: QWidget | None = None
        self.owner_rect: QRect | None = None

    def show_text(
        self,
        text: str,
        global_pos: QPoint,
        owner: QWidget,
        owner_rect: QRect | None = None,
        *,
        title: str | None = None,
        above: bool = False,
    ) -> None:
        self.owner = owner
        self.owner_rect = owner_rect

        self._title.setVisible(bool(title))
        if title:
            self._fit(self._title, title)
        # Under a title the text is secondary detail, so it's dimmed.
        body_color = TOOLTIP_TEXT_MUTED if title else TOOLTIP_TEXT
        self._body.setStyleSheet(label_style(f"color: {body_color}; font-size: 12px;"))
        self._fit(self._body, text)
        self.adjustSize()
        self.move(self._position(global_pos, above))

        if self.isVisible():
            # Already showing (e.g. moving between table headers or scatter
            # dots): swap content in place rather than fading in again.
            self.update()
            return
        self.setWindowOpacity(0.0)
        self.show()
        self._fade.start()

    def hide_tip(self) -> None:
        self._fade.stop()
        self.hide()
        self.owner = None
        self.owner_rect = None

    def _fit(self, label: QLabel, text: str) -> None:
        """Single line when it fits, otherwise wrapped at `_MAX_TEXT_WIDTH`."""
        label.setText(text)
        label.ensurePolished()
        label.setWordWrap(False)
        label.setMinimumWidth(0)
        label.setMaximumWidth(16777215)
        if label.fontMetrics().horizontalAdvance(text) > self._MAX_TEXT_WIDTH:
            label.setWordWrap(True)
            label.setFixedWidth(self._MAX_TEXT_WIDTH)

    def _position(self, global_pos: QPoint, above: bool) -> QPoint:
        """Top-left for the window so the card sits below the cursor (or centred
        above an anchor point), kept on screen."""
        margin = self._SHADOW
        width, height = self.width(), self.height()
        screen = QGuiApplication.screenAt(global_pos) or QGuiApplication.primaryScreen()
        area = screen.availableGeometry()

        if above:
            x = global_pos.x() - width // 2
            y = global_pos.y() - self._ANCHOR_GAP - height + margin
            if y + margin < area.top():
                y = global_pos.y() + self._ANCHOR_GAP - margin
        else:
            x = global_pos.x() + self._CURSOR_OFFSET.x() - margin
            y = global_pos.y() + self._CURSOR_OFFSET.y() - margin
            if y + height - margin > area.bottom():
                y = global_pos.y() - 8 - height + margin

        x = max(area.left() - margin, min(x, area.right() - width + margin))
        return QPoint(x, y)

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        margin = self._SHADOW
        card = QRectF(self.rect()).adjusted(margin, margin, -margin, -margin)

        # Soft shadow: stacked, slightly offset rounded rects, faintest outermost.
        painter.setPen(Qt.PenStyle.NoPen)
        steps = margin - 3
        for step in range(steps, 0, -1):
            alpha = int(26 * (1 - step / steps) ** 2) + 1
            painter.setBrush(QColor(0, 0, 0, alpha))
            painter.drawRoundedRect(
                card.adjusted(-step, -step + 3, step, step + 3), self._RADIUS + step, self._RADIUS + step
            )

        painter.setBrush(QColor(TOOLTIP_BACKGROUND))
        painter.setPen(QPen(QColor(255, 255, 255, 22), 1))
        painter.drawRoundedRect(card, self._RADIUS, self._RADIUS)


def _item_view_tip(view: QAbstractItemView, pos: QPoint) -> tuple[str, QRect | None]:
    """The tooltip under `pos` (viewport coordinates) of a table or its header."""
    if isinstance(view, QHeaderView):
        section = view.logicalIndexAt(pos)
        if section < 0:
            return "", None
        text = view.model().headerData(section, view.orientation(), Qt.ItemDataRole.ToolTipRole)
        start, size = view.sectionViewportPosition(section), view.sectionSize(section)
        if view.orientation() == Qt.Orientation.Horizontal:
            rect = QRect(start, 0, size, view.viewport().height())
        else:
            rect = QRect(0, start, view.viewport().width(), size)
        return (str(text) if text else ""), rect

    index = view.indexAt(pos)
    if not index.isValid():
        return "", None
    text = index.data(Qt.ItemDataRole.ToolTipRole)
    return (str(text) if text else ""), view.visualRect(index)


class _ToolTipController(QObject):
    """App-wide event filter that shows and dismisses the custom tooltip."""

    def __init__(self, tip: _ToolTipWindow, parent: QObject) -> None:
        super().__init__(parent)
        self._tip = tip

    def eventFilter(self, watched, event) -> bool:  # noqa: N802
        event_type = event.type()
        if event_type == QEvent.Type.ToolTip:
            return isinstance(watched, QWidget) and self._show_for(watched, event)

        if not self._tip.isVisible():
            return False
        if event_type == QEvent.Type.MouseMove:
            if not self._in_owner(event.globalPosition().toPoint()):
                self._tip.hide_tip()
        elif event_type in (QEvent.Type.Leave, QEvent.Type.Hide) and watched is self._tip.owner:
            self._tip.hide_tip()
        elif event_type in (
            QEvent.Type.MouseButtonPress,
            QEvent.Type.MouseButtonDblClick,
            QEvent.Type.Wheel,
            QEvent.Type.WindowDeactivate,
        ):
            self._tip.hide_tip()
        elif event_type == QEvent.Type.KeyPress:
            self._tip.hide_tip()
            # Consume Escape so the same press doesn't also, say, cancel a dialog.
            return event.key() == Qt.Key.Key_Escape
        return False

    def _show_for(self, widget: QWidget, event) -> bool:
        """Show `widget`'s tooltip, if it has one. Returning False lets Qt pass
        the event on to the parent widget, which comes back through here."""
        parent = widget.parentWidget()
        if isinstance(parent, QAbstractItemView) and widget is parent.viewport():
            text, rect = _item_view_tip(parent, event.pos())
        else:
            text, rect = widget.toolTip(), None
        if not text:
            return False
        self._tip.show_text(text, event.globalPos(), widget, rect)
        return True

    def _in_owner(self, global_pos: QPoint) -> bool:
        """Whether `global_pos` is still over the tooltip's owner (or its area)."""
        owner = self._tip.owner
        if owner is None or not shiboken6.isValid(owner):
            return False
        area = self._tip.owner_rect or owner.rect()
        return area.contains(owner.mapFromGlobal(global_pos))


_tip: _ToolTipWindow | None = None


def install_tooltips(app: QApplication) -> None:
    """Route every tooltip in the app through the custom tooltip window."""
    global _tip
    _tip = _ToolTipWindow()
    # Parented to `app` so it lives as long as the application does.
    controller = _ToolTipController(_tip, app)
    app.installEventFilter(controller)

    # Detach before teardown: otherwise Qt keeps feeding shutdown events into
    # this Python filter while the Python side is being dismantled, which
    # recurses until it errors on exit. `atexit` covers scripts that never
    # call `app.exec()`; it runs before PySide's own shutdown hook, which was
    # registered earlier.
    def _uninstall() -> None:
        if shiboken6.isValid(app) and shiboken6.isValid(controller):
            app.removeEventFilter(controller)
        if _tip is not None and shiboken6.isValid(_tip):
            _tip.hide_tip()

    app.aboutToQuit.connect(_uninstall)
    atexit.register(_uninstall)


def show_tooltip(
    text: str,
    global_pos: QPoint,
    owner: QWidget,
    *,
    title: str | None = None,
    above: bool = False,
) -> None:
    """Show the tooltip for `owner` directly, bypassing Qt's hover delay.

    Args:
        text: Tooltip text (the secondary line when `title` is given).
        global_pos: Cursor position, or with `above`, the point to centre the
            tooltip over.
        owner: The widget it belongs to; leaving it dismisses the tooltip.
        title: Optional bold first line.
        above: Centre the tooltip above `global_pos` instead of below-right.
    """
    if _tip is not None:
        _tip.show_text(text, global_pos, owner, title=title, above=above)


def hide_tooltip(owner: QWidget | None = None) -> None:
    """Hide the tooltip -- only if it belongs to `owner`, when given."""
    if _tip is not None and (owner is None or _tip.owner is owner):
        _tip.hide_tip()
