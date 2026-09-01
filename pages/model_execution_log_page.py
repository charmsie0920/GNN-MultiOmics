"""Model Execution Log page.

Shown while a model run is in progress. It surfaces two live views side by
side: a per-file pipeline status table and a scrolling terminal-style log
panel that streams the running backend's real stdout.

Main UI components:
    - Shared sidebar (`widgets.navigation.build_sidebar`) with "Model Logs"
      marked as the active nav item.
    - Shared header bar (`widgets.navigation.build_header_bar`) with
      "Model Running" marked as the active workflow tab.
    - An "Active Pipeline Elements" card containing a status table built
      with the shared table helpers in `widgets.tables`.
    - A "Model logs" card containing a dark terminal-style scroll panel that
      appends lines as they arrive from the backend.

Interactions with other pages:
    - `on_upload_clicked` navigates back to `DatasetInitializationPage`.
    - `on_model_visualization_clicked` navigates to `ModelVisualizationPage`.
    - `on_run_complete` is invoked once the watched run finishes
      successfully, and should navigate to the results page.
    All three callbacks are supplied and wired by `UILauncher.py`.

Call `start_watching(run_info)` (with the dict returned by
`client.api_client.start_run`, merged with the upload response) after
constructing this page to begin polling `client.workers.RunStatusPoller`
for live log lines and progress.
"""

from __future__ import annotations

import time
from collections.abc import Callable

from PySide6.QtCore import QEasingCurve, QPropertyAnimation, QTimer, Qt
from PySide6.QtWidgets import (
    QFrame,
    QGraphicsOpacityEffect,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpacerItem,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)

from client.workers import RunStatusPoller
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

# Flat estimate for a run's total duration, used to drive the progress bar
# and pipeline row fill. Replaces the old epoch-based estimate (current_epoch
# / max_epochs), which was inaccurate when epoch counts varied run to run.
_ESTIMATED_RUN_SECONDS = 12 * 60.0

# (filename, status label, state) — the real inputs cross_attention_baseline.py
# reads, in the order it loads them. Size/percent aren't known ahead of time
# so every row starts pending and flips to done together once the run
# completes (no per-line attribution to individual rows).
_PIPELINE_ROWS = [
    ("transcriptomics_pca.csv", "PENDING", "pending"),
    ("genomics_pca.csv", "PENDING", "pending"),
    ("proteomics_pca.csv", "PENDING", "pending"),
    ("gdsc2_response_master.csv", "PENDING", "pending"),
    ("Cross-Attention Fusion Training", "PENDING", "pending"),
]


class ModelExecutionLogPage(QWidget):
    """Live view of an in-progress model run: per-file pipeline status and logs."""

    def __init__(
        self,
        parent: QWidget | None = None,
        on_upload_clicked: Callable[[], None] | None = None,
        on_model_visualization_clicked: Callable[[], None] | None = None,
        on_run_complete: Callable[[], None] | None = None,
    ) -> None:
        """Build the page.

        Args:
            parent: Optional Qt parent widget.
            on_upload_clicked: Invoked when the header's "Upload" tab is
                clicked; should navigate back to the dataset upload page.
            on_model_visualization_clicked: Invoked when the sidebar's
                "Model Visualization" link is clicked; should navigate to
                the model visualization page.
            on_run_complete: Invoked once the watched run finishes
                successfully; should navigate to the results page.
        """
        super().__init__(parent)
        self.setObjectName("ModelExecutionLogPage")
        self._on_upload_clicked = on_upload_clicked
        self._on_model_visualization_clicked = on_model_visualization_clicked
        self._on_run_complete = on_run_complete

        self._poller: RunStatusPoller | None = None
        self._current_run_id: str | None = None

        # Keeps animations alive; QPropertyAnimation is not retained by Qt
        # once its local variable goes out of scope.
        self._pulse_animation: QPropertyAnimation | None = None
        self._last_pulse_label: QLabel | None = None

        self._progress_container: QVBoxLayout | None = None
        self._progress_bar_slot: QVBoxLayout | None = None
        self._progress_label: QLabel | None = None
        self._remaining_label: QLabel | None = None

        self._pipeline_table: QTableWidget | None = None
        self._log_layout: QVBoxLayout | None = None
        self._log_scroll: QScrollArea | None = None

        self._build_ui()

    # -- live updates ---------------------------------------------------

    def start_watching(self, run_info: dict) -> None:
        """Begin polling and displaying live status for a newly started run.

        Args:
            run_info: The merged upload+run-start response from
                `client.workers.UploadWorker` (must include `run_id`; the
                progress bar is driven by a flat `_ESTIMATED_RUN_SECONDS`
                estimate rather than anything in this dict).
        """
        if self._poller is not None:
            self._poller.stop()
            try:
                self._poller.statusUpdate.disconnect(self._on_status_update)
            except (RuntimeError, TypeError):
                pass

        self._current_run_id = run_info["run_id"]
        self._reset_log()
        self._reset_pipeline()
        self._update_progress_ui(0.0, 0.0, "running")

        self._poller = RunStatusPoller(run_info["run_id"], parent=self)
        self._poller.statusUpdate.connect(self._on_status_update)
        self._poller.start()

    def _on_status_update(self, payload: dict) -> None:
        # A poller for a run we've since replaced (e.g. a new upload started
        # while this one was still finishing up) may still have one update
        # in flight when it's told to stop; ignore anything not for the run
        # currently being watched so a late "completed" from an old run
        # can't jump the page to results out from under the user.
        if payload.get("run_id") != self._current_run_id:
            return

        for line in payload.get("new_log_lines", []):
            self._append_log_line(line)

        status = payload.get("status", "running")
        elapsed = payload.get("elapsed_seconds", 0.0)
        fraction = max(0.0, min(1.0, elapsed / _ESTIMATED_RUN_SECONDS))
        self._update_progress_ui(fraction, elapsed, status)

        if status == "running":
            self._update_pipeline_progress(fraction)
        elif status == "completed":
            self._mark_pipeline_done()
            if self._on_run_complete is not None:
                self._on_run_complete()
        elif status == "failed":
            QMessageBox.critical(
                self, "Model Run Failed", payload.get("error_message") or "Unknown error."
            )

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
            # !Hidden: ModelVisualizationPage is a decorative placeholder
            # (fake canvas + hardcoded run stats), hidden from navigation
            # until it's wired to real data.
            # make_sidebar_nav_button(
            #     "Model Visualization",
            #     "insights",
            #     callback=self._on_model_visualization_clicked,
            # ),
            make_sidebar_nav_button("Model Logs", "terminal", active=True),
        ]
        footer_widgets = [make_sidebar_nav_button("Support", "help_outline")]
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
        subtitle = QLabel("Real-time status of the model run against the uploaded dataset.")
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

        bar_slot = QVBoxLayout()
        bar_slot.setContentsMargins(0, 0, 0, 0)
        bar_slot.setSpacing(0)
        bar_slot.addWidget(self._build_overall_progress_bar(percent=0))
        self._progress_bar_slot = bar_slot
        page_header.addLayout(bar_slot)

        labels_row = QHBoxLayout()
        progress_label = QLabel("Overall Progress: 0%")
        remaining_label = QLabel("Est. Time Remaining: --")
        for label in (progress_label, remaining_label):
            label.setStyleSheet(f"font-size: 13px; font-family: Consolas, monospace; color: {TEXT_MUTED};")
        self._progress_label = progress_label
        self._remaining_label = remaining_label
        labels_row.addWidget(progress_label)
        labels_row.addStretch(1)
        labels_row.addWidget(remaining_label)
        page_header.addLayout(labels_row)

        return page_header

    @staticmethod
    def _build_overall_progress_bar(percent: int) -> QFrame:
        """Build the thin, page-level run-completion progress bar."""
        percent = max(0, min(100, percent))
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

    def _update_progress_ui(self, fraction: float, elapsed: float, status: str) -> None:
        """Refresh the progress bar fill and the progress/remaining-time labels.

        Progress is a flat estimate of `elapsed / _ESTIMATED_RUN_SECONDS`
        (15 minutes) rather than the previous epoch-count-based fraction,
        which was inaccurate when epoch counts varied run to run. Capped at
        99% while still running so it never visually completes early; jumps
        to 100% once `status` is `"completed"`.
        """
        fraction = max(0.0, min(1.0, fraction))
        if status == "completed":
            clamped = 100
        else:
            clamped = min(99, round(fraction * 100))

        if self._progress_bar_slot is not None:
            while self._progress_bar_slot.count():
                item = self._progress_bar_slot.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
            self._progress_bar_slot.addWidget(self._build_overall_progress_bar(percent=clamped))

        if self._progress_label is not None:
            self._progress_label.setText(f"Overall Progress: {clamped}%")

        if self._remaining_label is not None:
            if status == "running":
                remaining = max(0.0, _ESTIMATED_RUN_SECONDS - elapsed)
                minutes, seconds = divmod(int(remaining), 60)
                self._remaining_label.setText(f"Est. Time Remaining: {minutes}m {seconds:02d}s")
            elif status == "completed":
                self._remaining_label.setText("Est. Time Remaining: Complete")
            else:
                self._remaining_label.setText("Est. Time Remaining: --")

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
        badge = QLabel("Cross-Attention Fusion")
        badge.setStyleSheet(
            f"background: {SURFACE_CONTAINER}; border: none; color: {TEXT_MUTED}; border-radius: 4px;"
            " padding: 4px 8px; font-size: 12px; font-weight: 700; letter-spacing: 0.05em;"
        )
        header_layout.addWidget(badge)
        layout.addWidget(header)

        table = QTableWidget(len(_PIPELINE_ROWS), 3)
        table.setHorizontalHeaderLabels(["File Identifier", "Current State", "Pipeline Progress"])
        style_data_table(table, header_background=WINDOW_BACKGROUND)
        self._pipeline_table = table
        self._populate_pipeline_table()

        table.setColumnWidth(0, 260)
        table.setColumnWidth(1, 190)
        layout.addWidget(table)
        return card

    def _populate_pipeline_table(self) -> None:
        table = self._pipeline_table
        if table is None:
            return
        for row_index, (filename, status, state) in enumerate(_PIPELINE_ROWS):
            table.setCellWidget(row_index, 0, self._build_file_cell(filename, state))
            icon_name = "sync" if state == "active" else None
            table.setCellWidget(row_index, 1, build_status_badge(status, _STATE_TO_BADGE_TONE[state], icon_name=icon_name))
            table.setCellWidget(row_index, 2, build_mini_progress_bar(100 if state == "done" else 0))

    def _reset_pipeline(self) -> None:
        self._populate_pipeline_table()

    def _update_pipeline_progress(self, fraction: float) -> None:
        """Advance the pipeline rows as a cosmetic proxy for real progress.

        There's no per-file progress signal from the backend, so this ties
        each row's completion to the same elapsed/expected fraction driving
        the overall progress bar. The fraction is split evenly across the
        rows, and the single row currently "in progress" gets its own mini
        bar filled to its within-row fraction — so it fills gradually over
        roughly a fifth of the total run instead of snapping straight to
        full — rather than every row jumping instantly between 0% and 100%.
        `_mark_pipeline_done` immediately fast-forwards every row to "done"
        the moment the run actually finishes (which may be earlier or later
        than this estimate).
        """
        table = self._pipeline_table
        if table is None:
            return
        fraction = max(0.0, min(1.0, fraction))
        row_count = len(_PIPELINE_ROWS)
        scaled = fraction * row_count

        for row_index, (filename, _status, _state) in enumerate(_PIPELINE_ROWS):
            row_percent = max(0.0, min(100.0, (scaled - row_index) * 100.0))
            if row_percent >= 100.0:
                state, status_text = "done", "DONE"
            elif row_percent <= 0.0:
                state, status_text = "pending", "PENDING"
            else:
                state, status_text = "active", "IN PROGRESS"

            table.setCellWidget(row_index, 0, self._build_file_cell(filename, state))
            icon_name = "sync" if state == "active" else None
            table.setCellWidget(row_index, 1, build_status_badge(status_text, _STATE_TO_BADGE_TONE[state], icon_name=icon_name))
            table.setCellWidget(row_index, 2, build_mini_progress_bar(round(row_percent)))

    def _mark_pipeline_done(self) -> None:
        table = self._pipeline_table
        if table is None:
            return
        for row_index, (filename, _status, _state) in enumerate(_PIPELINE_ROWS):
            table.setCellWidget(row_index, 0, self._build_file_cell(filename, "done"))
            table.setCellWidget(row_index, 1, build_status_badge("DONE", _STATE_TO_BADGE_TONE["done"]))
            table.setCellWidget(row_index, 2, build_mini_progress_bar(100))

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

    def _build_log_body(self) -> QScrollArea:
        """Build the self-scrolling, appendable console panel."""
        content = QFrame()
        content.setStyleSheet("background: #1e1e1e;")
        log_layout = QVBoxLayout(content)
        log_layout.setContentsMargins(16, 16, 16, 16)
        log_layout.setSpacing(6)
        log_layout.addStretch(1)
        self._log_layout = log_layout

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setStyleSheet("background: #1e1e1e; border: none;")
        scroll.setWidget(content)
        self._log_scroll = scroll
        return scroll

    def _reset_log(self) -> None:
        if self._log_layout is None:
            return
        while self._log_layout.count() > 1:
            item = self._log_layout.takeAt(0)
            layout = item.layout()
            if layout is not None:
                while layout.count():
                    child = layout.takeAt(0)
                    widget = child.widget()
                    if widget is not None:
                        widget.deleteLater()
        self._last_pulse_label = None

    def _append_log_line(self, text: str) -> None:
        """Append one streamed log line to the console, pulsing the newest."""
        if self._log_layout is None:
            return

        if self._last_pulse_label is not None:
            self._last_pulse_label.setGraphicsEffect(None)

        row = QHBoxLayout()
        row.setSpacing(16)
        time_label = QLabel(time.strftime("[%H:%M:%S]"))
        time_label.setStyleSheet("color: #71717a; font-family: Consolas, monospace; font-size: 13px;")
        message_label = QLabel(text)
        message_label.setWordWrap(True)
        message_label.setStyleSheet("color: #d4d4d8; font-family: Consolas, monospace; font-size: 13px;")
        row.addWidget(time_label, 0)
        row.addWidget(message_label, 1)

        self._log_layout.insertLayout(self._log_layout.count() - 1, row)
        self._start_pulse_animation(message_label)
        self._last_pulse_label = message_label

        if self._log_scroll is not None:
            QTimer.singleShot(0, self._scroll_log_to_bottom)

    def _scroll_log_to_bottom(self) -> None:
        if self._log_scroll is not None:
            bar = self._log_scroll.verticalScrollBar()
            bar.setValue(bar.maximum())

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
