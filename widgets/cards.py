from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QPainter
from PySide6.QtWidgets import (
    QFrame,
    QGraphicsDropShadowEffect,
    QHBoxLayout,
    QLabel,
    QRadioButton,
    QVBoxLayout,
    QWidget,
)


class SurfaceCard(QFrame):
    def __init__(self, parent: QWidget | None = None, object_name: str = "SurfaceCard") -> None:
        super().__init__(parent)
        self.setObjectName(object_name)
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        shadow = QGraphicsDropShadowEffect(self)
        shadow.setBlurRadius(24)
        shadow.setOffset(0, 8)
        shadow.setColor(QColor(0, 0, 0, 18))
        self.setGraphicsEffect(shadow)


class UploadCard(SurfaceCard):
    fileDropped = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent, object_name="UploadCard")
        self.setAcceptDrops(True)

    def dragEnterEvent(self, event) -> None:  # noqa: N802
        if self._first_csv_path(event.mimeData()) is not None:
            event.acceptProposedAction()
        else:
            event.ignore()

    def dropEvent(self, event) -> None:  # noqa: N802
        path = self._first_csv_path(event.mimeData())
        if path is not None:
            self.fileDropped.emit(path)
            event.acceptProposedAction()
        else:
            event.ignore()

    @staticmethod
    def _first_csv_path(mime_data) -> str | None:
        for url in mime_data.urls():
            if url.isLocalFile() and url.toLocalFile().lower().endswith(".csv"):
                return url.toLocalFile()
        return None

    def paintEvent(self, event) -> None:  # noqa: N802
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(18, 18, 18, 14))
        width = self.width()
        height = self.height()
        spacing = 20
        radius = 1.5
        for y in range(0, height, spacing):
            for x in range(0, width, spacing):
                painter.drawEllipse(x, y, radius * 2, radius * 2)


class OptionTile(SurfaceCard):
    clicked = Signal()

    def __init__(self, title: str, checked: bool = False, parent: QWidget | None = None) -> None:
        super().__init__(parent, object_name="OptionTile")
        self._title_text = title
        self._selected = checked
        self.setProperty("selected", checked)

        self._indicator = QFrame()
        self._indicator.setObjectName("OptionIndicator")
        self._indicator.setFixedSize(20, 20)
        self._indicator.setProperty("selected", checked)
        self._indicator.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)

        self._title = QLabel(title)
        self._title.setWordWrap(True)
        self._title.setStyleSheet(
            "font-size: 16px; font-weight: 600; color: #1a1c1c;"
        )

        row = QHBoxLayout(self)
        row.setContentsMargins(14, 12, 14, 12)
        row.setSpacing(12)
        row.addWidget(self._indicator, 0, Qt.AlignmentFlag.AlignTop)
        row.addWidget(self._title, 1)

        self._sync_state()

    def mousePressEvent(self, event) -> None:  # noqa: N802
        self.clicked.emit()
        super().mousePressEvent(event)

    def setSelected(self, selected: bool) -> None:
        self._selected = selected
        self.setProperty("selected", selected)
        self._indicator.setProperty("selected", selected)
        self._sync_state()

    def isSelected(self) -> bool:
        return self._selected

    def _sync_state(self) -> None:
        self.style().unpolish(self)
        self.style().polish(self)
        self._indicator.style().unpolish(self._indicator)
        self._indicator.style().polish(self._indicator)
        self._title.setStyleSheet(
            "font-size: 16px; font-weight: 600; color: #1a1c1c;"
            if self._selected
            else "font-size: 16px; font-weight: 400; color: #1a1c1c;"
        )
