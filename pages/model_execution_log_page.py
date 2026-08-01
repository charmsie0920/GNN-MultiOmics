"""Model Execution Log page.

Shown while a model run is in progress. It surfaces two live views side by
side: a per-file pipeline status table (which sequencing files have cleared
QC, which are still normalizing/cleaning, which haven't started) and a
scrolling terminal-style log panel.

Main UI components:
    - Shared sidebar (`widgets.navigation.build_sidebar`) with "Model Logs"
      marked as the active nav item.
    - Shared header bar (`widgets.navigation.build_header_bar`) with
      "Model Running" marked as the active workflow tab.
    - An "Active Pipeline Elements" card containing a status table built
      with the shared table helpers in `widgets.tables`.
    - A "Model logs" card containing a dark terminal-style scroll panel with
      severity-colored lines and a pulsing "live" line.

Interactions with other pages:
    - `on_upload_clicked` navigates back to `DatasetInitializationPage`.
    - `on_model_visualization_clicked` navigates to `ModelVisualizationPage`.
    Both callbacks are supplied and wired by `main.py`.

All pipeline rows and log lines are placeholder data until a backend
execution service is connected.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import QEasingCurve, QPropertyAnimation, Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QFrame,
    QGraphicsOpacityEffect,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpacerItem,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from styles.theme import (
    BORDER,
    CARD_CONTAINER_STYLE,
    CARD_TITLE_STYLE,
    PAGE_MARGIN,
    PAGE_SUBTITLE_STYLE,
    PAGE_TITLE_STYLE,
    PRIMARY,
    PRIMARY_BUTTON_STYLE,
    SECONDARY_BUTTON_STYLE,
    SECTION_SPACING,
    SURFACE_CONTAINER,
    TEXT,
    TEXT_FAINT,
    TEXT_MUTED,
    WINDOW_BACKGROUND,
)
from widgets.icons import icon_text
from widgets.navigation import (
    build_header_bar,
    build_sidebar,
    make_primary_cta_button,
    make_sidebar_nav_button,
    make_top_tab,
)
from widgets.tables import build_mini_progress_bar, build_status_badge, style_data_table, transparent_cell_widget

# Pipeline row state -> status badge tone (see widgets.tables.build_status_badge).
_STATE_TO_BADGE_TONE = {"done": "positive", "active": "neutral", "pending": "muted"}

# (filename, size, status label, completion percent, state)
_PIPELINE_ROWS = [
    # * Replace with records fetched from backend /api/v1/model/logs when pipeline is integrated.
    ("SEQ_001_A.fq.gz", "2.4GB", "QC PASSED", 100, "done"),
    ("SEQ_002_A.fq.gz", "2.1GB", "NORMALIZING", 80, "active"),
    ("SEQ_003_B.fq.gz", "3.0GB", "CLEANING", 35, "active"),
    ("SEQ_004_C.fq.gz", "1.8GB", "PENDING", 0, "pending"),
    ("SEQ_005_C.fq.gz", "2.2GB", "PENDING", 0, "pending"),
]

# (timestamp, message, severity style key)
_LOG_LINES = [
    # * Replace with streaming backend logs once execution service is wired.
    ("[10:04:12]", "INFO: Initializing Preprocessing Engine v2.4.1", "info"),
    ("[10:04:15]", "INFO: Loading parameter configuration 'Standard_Clinical_V2'", "info"),
    ("[10:04:16]", "INFO: Starting batch RUN_A7B9 (24 items)", "info"),
    ("[10:04:18]", "PROCESS: SEQ_001_A.fq.gz -> Cleaning sequence...", "process"),
    ("[10:05:42]", "PROCESS: SEQ_001_A.fq.gz -> Normalizing Read Depth (Target: 30x)", "process"),
    ("[10:08:15]", "WARN: Mild GC bias detected in SEQ_001_A. Applying correction matrix.", "warn"),
    ("[10:12:01]", "PROCESS: SEQ_001_A.fq.gz -> Executing Quality Control bounds check.", "process"),
    ("[10:12:30]", "SUCCESS: SEQ_001_A.fq.gz QC Passed. Phred>30 = 94.2%", "process"),
    ("[10:12:35]", "PROCESS: SEQ_002_A.fq.gz -> Cleaning sequence...", "process"),
    ("[10:14:10]", "PROCESS: SEQ_002_A.fq.gz -> Normalizing Read Depth (Target: 30x)", "process"),
    ("[10:14:15]", "PROCESS: SEQ_003_B.fq.gz -> Cleaning sequence...", "process"),
    ("[10:15:00]", "INFO: Thread pool utilization at 85%. Memory stable at 14.2GB.", "info"),
    ("[10:16:22]", "Running normalization algorithm on block 4/12...", "pulse"),
]

# Severity -> log line text color.
_LOG_MESSAGE_COLOR = {
    "info": "#d4d4d8",
    "warn": "#d4d4d8",
    "process": "#ffffff",
    "pulse": "#a1a1aa",
}


class ModelExecutionLogPage(QWidget):
    """Live view of an in-progress model run: per-file pipeline status and logs."""

    def __init__(
        self,
        parent: QWidget | None = None,
        on_upload_clicked: Callable[[], None] | None = None,
        on_model_visualization_clicked: Callable[[], None] | None = None,
    ) -> None:
        """Build the page.

        Args:
            parent: Optional Qt parent widget.
            on_upload_clicked: Invoked when the header's "Upload" tab is
                clicked; should navigate back to the dataset upload page.
            on_model_visualization_clicked: Invoked when the sidebar's
                "Model Visualization" link is clicked; should navigate to
                the model visualization page.
        """
        super().__init__(parent)
        self.setObjectName("ModelExecutionLogPage")
        # * UI and table rows are placeholders until backend model-run data is connected.
        self._on_upload_clicked = on_upload_clicked
        self._on_model_visualization_clicked = on_model_visualization_clicked
        # Keeps the pulse animation alive; QPropertyAnimation is not retained
        # by Qt once its local variable goes out of scope.
        self._pulse_animation: QPropertyAnimation | None = None
        self._build_ui()

    def _build_ui(self) -> None:
        """Lay out the sidebar, header, and scrollable body content."""
        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        root.addWidget(self._build_sidebar(), 0)

        content_shell = QFrame()
        content_shell.setObjectName("RootShell")
        content_shell_layout = QVBoxLayout(content_shell)
        content_shell_layout.setContentsMargins(0, 0, 0, 0)
        content_shell_layout.setSpacing(0)
        content_shell.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        content_shell_layout.addWidget(self._build_header(), 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(PAGE_MARGIN, PAGE_MARGIN, PAGE_MARGIN, PAGE_MARGIN)
        content_layout.setSpacing(SECTION_SPACING)

        content_layout.addLayout(self._build_page_header())

        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(SECTION_SPACING)
        grid.setVerticalSpacing(SECTION_SPACING)
        grid.setColumnStretch(0, 2)
        grid.setColumnStretch(1, 1)
        grid.addWidget(self._build_pipeline_card(), 0, 0)
        grid.addWidget(self._build_log_card(), 0, 1)
        content_layout.addLayout(grid)

        content_layout.addItem(QSpacerItem(0, 0, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))
        scroll.setWidget(content)
        content_shell_layout.addWidget(scroll, 1)
        root.addWidget(content_shell, 1)

    def _build_sidebar(self) -> QFrame:
        """Build the shared sidebar with "Model Logs" as the active nav item."""
        nav_widgets = [
            make_sidebar_nav_button(
                "Model Visualization",
                "insights",
                callback=self._on_model_visualization_clicked,
            ),
            make_sidebar_nav_button("Model Logs", "terminal", active=True),
        ]
        footer_widgets = [
            make_sidebar_nav_button("History", "history"),
            make_sidebar_nav_button("Settings", "settings"),
            make_sidebar_nav_button("Support", "help_outline"),
        ]
        return build_sidebar(
            nav_widgets,
            footer_widgets,
            cta_widget=make_primary_cta_button(),
        )

    def _build_header(self) -> QFrame:
        """Build the shared header bar with "Model Running" as the active tab."""
        tabs = [
            make_top_tab("Upload", callback=self._on_upload_clicked),
            make_top_tab("Model Running", active=True),
            make_top_tab("Results"),
        ]
        return build_header_bar(tabs)

    def _build_page_header(self) -> QVBoxLayout:
        """Build the title/subtitle row, action buttons, and overall progress bar."""
        page_header = QVBoxLayout()
        page_header.setSpacing(8)

        title_row = QHBoxLayout()
        title_col = QVBoxLayout()
        title_col.setSpacing(4)
        title = QLabel("Model Execution Log")
        title.setStyleSheet(PAGE_TITLE_STYLE)
        subtitle = QLabel("Real-time status of neural network training and validation across distributed compute nodes.")
        subtitle.setStyleSheet(PAGE_SUBTITLE_STYLE)
        title_col.addWidget(title)
        title_col.addWidget(subtitle)
        title_row.addLayout(title_col)
        title_row.addStretch(1)

        action_row = QHBoxLayout()
        action_row.setSpacing(12)
        pause_button = QPushButton("Pause Execution")
        pause_button.setStyleSheet(SECONDARY_BUTTON_STYLE)
        halt_button = QPushButton("Halt Execution")
        halt_button.setStyleSheet(PRIMARY_BUTTON_STYLE)
        action_row.addWidget(pause_button)
        action_row.addWidget(halt_button)
        title_row.addLayout(action_row)
        page_header.addLayout(title_row)

        page_header.addWidget(self._build_overall_progress_bar(percent=45))

        labels_row = QHBoxLayout()
        progress_label = QLabel("Overall Progress: 45%")
        remaining_label = QLabel("Est. Time Remaining: 12m 30s")
        for label in (progress_label, remaining_label):
            label.setStyleSheet(f"font-size: 13px; font-family: Consolas, monospace; color: {TEXT_MUTED};")
        labels_row.addWidget(progress_label)
        labels_row.addStretch(1)
        labels_row.addWidget(remaining_label)
        page_header.addLayout(labels_row)

        return page_header

    @staticmethod
    def _build_overall_progress_bar(percent: int) -> QFrame:
        """Build the thin, page-level run-completion progress bar."""
        track = QFrame()
        track.setFixedHeight(4)
        track.setStyleSheet(f"background: {SURFACE_CONTAINER}; border: none;")
        track_layout = QHBoxLayout(track)
        track_layout.setContentsMargins(0, 0, 0, 0)
        track_layout.setSpacing(0)
        fill = QFrame()
        fill.setFixedHeight(4)
        fill.setStyleSheet(f"background: {PRIMARY}; border: none;")
        track_layout.addWidget(fill, percent)
        track_layout.addStretch(100 - percent)
        return track

    def _build_pipeline_card(self) -> QFrame:
        """Build the "Active Pipeline Elements" card with the per-file status table."""
        card = QFrame()
        card.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        header = QFrame()
        header.setStyleSheet(f"background: {WINDOW_BACKGROUND}; border-bottom: 1px solid {SURFACE_CONTAINER};")
        header_layout = QHBoxLayout(header)
        header_layout.setContentsMargins(16, 12, 16, 12)
        header_title = QLabel("Active Pipeline Elements")
        header_title.setStyleSheet(CARD_TITLE_STYLE)
        header_layout.addWidget(header_title)
        header_layout.addStretch(1)
        badge = QLabel("Queue Number")
        badge.setStyleSheet(
            f"background: {SURFACE_CONTAINER}; border: none; color: {TEXT_MUTED}; border-radius: 4px;"
            " padding: 4px 8px; font-size: 12px; font-weight: 700; letter-spacing: 0.05em;"
        )
        header_layout.addWidget(badge)
        layout.addWidget(header)

        table = QTableWidget(len(_PIPELINE_ROWS), 4)
        table.setHorizontalHeaderLabels(["File Identifier", "Size", "Current State", "Pipeline Progress"])
        style_data_table(table, header_background=WINDOW_BACKGROUND)

        for row_index, (filename, size, status, percent, state) in enumerate(_PIPELINE_ROWS):
            table.setCellWidget(row_index, 0, self._build_file_cell(filename, state))

            size_item = QTableWidgetItem(size)
            size_item.setForeground(QColor(TEXT_MUTED))
            table.setItem(row_index, 1, size_item)

            icon_name = "sync" if state == "active" else None
            table.setCellWidget(row_index, 2, build_status_badge(status, _STATE_TO_BADGE_TONE[state], icon_name=icon_name))
            table.setCellWidget(row_index, 3, build_mini_progress_bar(percent))

        table.setColumnWidth(0, 220)
        table.setColumnWidth(1, 100)
        table.setColumnWidth(2, 190)
        layout.addWidget(table)
        return card

    @staticmethod
    def _build_file_cell(filename: str, state: str) -> QWidget:
        """Build the "File Identifier" cell: a document icon plus filename.

        Args:
            filename: The file's display name.
            state: Pipeline row state ("done", "active", or "pending"),
                used to dim the icon/text for not-yet-started files.

        Returns:
            A `QWidget` suitable for `QTableWidget.setCellWidget`.
        """
        wrapper = transparent_cell_widget()
        row = QHBoxLayout(wrapper)
        row.setContentsMargins(16, 0, 8, 0)
        row.setSpacing(8)

        icon_color = PRIMARY if state == "active" else "#8a8d8d"
        icon = QLabel(icon_text("description"))
        icon.setStyleSheet(f"background: transparent; border: none; color: {icon_color}; font-size: 14px;")

        text_color = TEXT_FAINT if state == "pending" else TEXT
        label = QLabel(filename)
        label.setStyleSheet(
            f"background: transparent; border: none; font-family: Consolas, monospace; font-size: 13px; color: {text_color};"
        )

        row.addWidget(icon)
        row.addWidget(label)
        row.addStretch(1)
        return wrapper

    def _build_log_card(self) -> QFrame:
        """Build the dark terminal-style "Model logs" card."""
        card = QFrame()
        card.setStyleSheet(f"background: {SURFACE_CONTAINER}; border: 1px solid {BORDER}; border-radius: 12px;")
        card.setMinimumHeight(700)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        layout.addWidget(self._build_log_header())
        layout.addWidget(self._build_log_body(), 1)
        return card

    def _build_log_header(self) -> QFrame:
        """Build the log card's header: title, LIVE indicator, and download action."""
        header = QFrame()
        header.setStyleSheet(f"background: #e2e2e2; border-bottom: 1px solid {BORDER};")
        header_layout = QHBoxLayout(header)
        header_layout.setContentsMargins(16, 12, 16, 12)

        log_title = QLabel(f"{icon_text('terminal')}  Model logs")
        log_title.setStyleSheet(CARD_TITLE_STYLE)
        header_layout.addWidget(log_title)
        header_layout.addStretch(1)

        live_badge = QLabel("LIVE")
        live_badge.setStyleSheet(
            f"background: rgba(0,0,0,0.1); color: {PRIMARY}; border-radius: 999px; padding: 4px 8px;"
            " font-size: 10px; font-weight: 700;"
        )
        header_layout.addWidget(live_badge)

        download_button = QPushButton(icon_text("download"))
        download_button.setCursor(Qt.CursorShape.PointingHandCursor)
        download_button.setStyleSheet(f"background: transparent; border: none; color: {TEXT_MUTED}; font-size: 16px;")
        header_layout.addWidget(download_button)

        return header

    def _build_log_body(self) -> QFrame:
        """Build the scrollable, severity-colored log line list."""
        log_body = QFrame()
        log_body.setStyleSheet("background: #1e1e1e;")
        log_layout = QVBoxLayout(log_body)
        log_layout.setContentsMargins(16, 16, 16, 16)
        log_layout.setSpacing(6)

        for timestamp, message, style in _LOG_LINES:
            row = QHBoxLayout()
            row.setSpacing(16)
            time_label = QLabel(timestamp)
            time_label.setStyleSheet("color: #71717a; font-family: Consolas, monospace; font-size: 13px;")
            message_label = QLabel(message)
            message_label.setWordWrap(True)
            message_label.setStyleSheet(
                f"color: {_LOG_MESSAGE_COLOR[style]}; font-family: Consolas, monospace; font-size: 13px;"
            )
            if style == "pulse":
                self._start_pulse_animation(message_label)
            row.addWidget(time_label, 0)
            row.addWidget(message_label, 1)
            log_layout.addLayout(row)

        log_layout.addStretch(1)
        return log_body

    def _start_pulse_animation(self, label: QLabel) -> None:
        """Animate `label`'s opacity in a continuous breathing loop.

        Mirrors a CSS `animate-pulse` treatment for the most recent,
        still-running log line. The animation and its opacity effect are
        stored on `self` so they aren't garbage-collected once this method
        returns.

        Args:
            label: The log message label to animate.
        """
        effect = QGraphicsOpacityEffect(label)
        label.setGraphicsEffect(effect)
        animation = QPropertyAnimation(effect, b"opacity", label)
        animation.setDuration(1400)
        animation.setKeyValueAt(0.0, 1.0)
        animation.setKeyValueAt(0.5, 0.35)
        animation.setKeyValueAt(1.0, 1.0)
        animation.setEasingCurve(QEasingCurve.Type.InOutSine)
        animation.setLoopCount(-1)
        animation.start()
        self._pulse_animation = animation
