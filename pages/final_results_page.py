"""Final Results page.

The last stop in the workflow: shows the model's predicted drug rankings for
the selected target cell line, alongside supporting visualizations (training
curve, predicted IC50 vs. confidence) and detail panels (top drug details,
sample profile).

Main UI components:
    - Shared sidebar (`widgets.navigation.build_sidebar`) — this page has no
      "Model Visualization"/"Model Logs" nav section, matching the upload
      page's sidebar, since results is a terminal step in the workflow.
    - Shared header bar (`widgets.navigation.build_header_bar`) with
      "Results" marked as the active workflow tab.
    - A "Predicted Drug Results Panel" table, styled with the same shared
      table helpers (`widgets.tables`) as the model execution log page's
      pipeline table, so the two read as one design system. Populated from
      real backend inference via `load_results(run_id)`.
    - Two custom-painted charts fed from real run data: `TrainingCurveWidget`
      (train_loss/val_rmse per epoch, from `on_training_history`) and
      `PredictionScatterWidget` (predicted IC50 vs. confidence, one point
      per ranked drug, from the same rows as the results table). The
      network proximity box on the sample profile card is still a
      placeholder — out of scope for now.
    - Drug detail and sample profile summary cards, populated from the same
      real results.

Interactions with other pages:
    - `on_upload_clicked` navigates back to `DatasetInitializationPage`.
    - `on_model_running_clicked` navigates to `ModelExecutionLogPage`.
    Both callbacks are supplied and wired by `UILauncher.py`, which also
    calls `load_results(run_id)` / `set_sample_id(target_cell_line)` before
    switching to this page.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import (
    QComboBox,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from client.workers import ResultsWorker, TrainingHistoryWorker
from styles.theme import (
    BORDER,
    CARD_CONTAINER_STYLE,
    CARD_TITLE_STYLE,
    LABEL_CAPS_STYLE,
    PAGE_MARGIN,
    PAGE_SUBTITLE_STYLE,
    PAGE_TITLE_STYLE,
    PRIMARY,
    PRIMARY_BUTTON_STYLE,
    SECONDARY_BUTTON_STYLE,
    SECTION_SPACING,
    SURFACE,
    SURFACE_CONTAINER,
    SURFACE_HIGH,
    TEXT,
    TEXT_MUTED,
    WINDOW_BACKGROUND,
    label_style,
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

# Sensitivity ranking label -> status badge tone (see widgets.tables.build_status_badge).
_RANK_TO_BADGE_TONE = {"HIGH SENSITIVITY": "positive", "MEDIUM": "neutral", "LOW": "muted"}


def _format_ic50(value: float) -> str:
    """Format a predicted IC50 (uM) with precision that scales to its magnitude.

    A flat `.1f` rounds any sub-0.05 uM prediction (a real, very sensitive
    result -- GDSC ln_ic50 values span roughly -10 to 10, i.e. ~5e-5 to
    ~22000 uM) down to a misleading "0.0". Below 1 uM, show 3 significant
    figures instead so small-but-real values stay visible.
    """
    if value == 0:
        return "0.0"
    if abs(value) >= 1:
        return f"{value:.1f}"
    return f"{value:.3g}"


@dataclass(frozen=True)
class DrugResult:
    """One row of the predicted drug results table."""

    name: str
    ic50: float
    rank: str
    confidence: float


class TrainingCurveWidget(QWidget):
    """Custom-painted line chart of the run's real per-epoch train loss / val RMSE.

    Fed from `on_training_history` (backend/model_backends/cross_attention.py),
    parsed straight out of `train()`'s own
    "[epoch N] train_loss=... val_rmse=... val_pcc=..." log line -- these are
    the model's actual convergence numbers for this run, not mocked data.

    Each series is normalized to its own 0-1 range (their absolute scales
    aren't comparable -- loss and RMSE are different units), so the numeric
    value range for each is spelled out in its legend label instead of a
    shared y-axis, alongside epoch ticks on the x-axis.
    """

    # Warm/cool complementary pair so the two series stay visually distinct
    # even though they share the same normalized 0-1 vertical space.
    _TRAIN_LOSS_COLOR = "#8a3419"
    _VAL_RMSE_COLOR = "#1c4a7a"

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(185)
        self._history: list[dict] = []

    def set_history(self, history: list[dict]) -> None:
        self._history = history
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        """Paint a dashed gridline background plus the train-loss/val-RMSE curves."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor(WINDOW_BACKGROUND))

        margin_left, margin_right = 10, 10
        margin_top, margin_bottom = 22, 18
        width = max(1, self.width() - margin_left - margin_right)
        height = max(1, self.height() - margin_top - margin_bottom)

        grid_pen = QPen(QColor(SURFACE_CONTAINER))
        grid_pen.setStyle(Qt.PenStyle.DashLine)
        painter.setPen(grid_pen)
        for ratio in (0.25, 0.5, 0.75):
            y = margin_top + int(height * ratio)
            painter.drawLine(margin_left, y, margin_left + width, y)

        if len(self._history) < 2:
            painter.setPen(QColor(TEXT_MUTED))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Waiting for training history...")
            return

        epochs = [point["epoch"] for point in self._history]
        e_lo, e_hi = min(epochs), max(epochs)
        e_span = (e_hi - e_lo) or 1

        def normalize(values: list[float]) -> list[float]:
            lo, hi = min(values), max(values)
            span = (hi - lo) or 1.0
            return [(v - lo) / span for v in values]

        def map_xy(epoch: int, normalized_value: float) -> QPointF:
            px = (epoch - e_lo) / e_span
            py = 1.0 - normalized_value
            return QPointF(margin_left + (px * width), margin_top + (py * height))

        def draw_series(values: list[float], color: str) -> None:
            path = QPainterPath()
            path.moveTo(map_xy(epochs[0], values[0]))
            for epoch, value in zip(epochs[1:], values[1:]):
                path.lineTo(map_xy(epoch, value))
            pen = QPen(QColor(color))
            pen.setWidth(2)
            painter.setPen(pen)
            painter.drawPath(path)

        train_loss_values = [point["train_loss"] for point in self._history]
        val_rmse_values = [point["val_rmse"] for point in self._history]
        draw_series(normalize(train_loss_values), self._TRAIN_LOSS_COLOR)
        draw_series(normalize(val_rmse_values), self._VAL_RMSE_COLOR)

        # Legend with each series' real value range (its absolute numbers,
        # since the lines themselves are drawn normalized).
        painter.setPen(QColor(self._TRAIN_LOSS_COLOR))
        loss_label = f"● Train Loss {min(train_loss_values):.3f}–{max(train_loss_values):.3f}"
        painter.drawText(margin_left, 14, loss_label)
        loss_label_width = painter.fontMetrics().horizontalAdvance(loss_label)
        painter.setPen(QColor(self._VAL_RMSE_COLOR))
        rmse_label = f"● Val RMSE {min(val_rmse_values):.3f}–{max(val_rmse_values):.3f}"
        painter.drawText(margin_left + loss_label_width + 14, 14, rmse_label)

        # Epoch ticks along the x-axis.
        painter.setPen(QColor(TEXT_MUTED))
        painter.drawText(margin_left, margin_top + height + 14, f"Epoch {e_lo}")
        hi_text = f"Epoch {e_hi}"
        hi_text_width = painter.fontMetrics().horizontalAdvance(hi_text)
        painter.drawText(margin_left + width - hi_text_width, margin_top + height + 14, hi_text)


class PredictionScatterWidget(QWidget):
    """Custom-painted scatter chart: predicted IC50 vs. confidence, one point per ranked drug.

    Fed from the same ranked predictions already loaded into the results
    table (`FinalResultsPage._on_results_succeeded`) -- real per-drug output
    from this run, not mocked data. IC50 (uM) is plotted on a log scale
    since predictions span several orders of magnitude (see `_format_ic50`).
    Each point is colored along a low-to-high confidence gradient so the
    color itself carries information rather than being purely decorative.
    """

    # Low-confidence points read as amber/uncertain, high-confidence points
    # as green/trustworthy -- interpolated per point in `_confidence_color`.
    _LOW_CONFIDENCE_COLOR = QColor(133, 77, 24)
    _HIGH_CONFIDENCE_COLOR = QColor(27, 94, 54)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(185)
        self._points: list[tuple[float, float]] = []

    def set_points(self, points: list[tuple[float, float]]) -> None:
        self._points = points
        self.update()

    @classmethod
    def _confidence_color(cls, confidence_percent: float) -> QColor:
        t = max(0.0, min(1.0, confidence_percent / 100.0))
        low, high = cls._LOW_CONFIDENCE_COLOR, cls._HIGH_CONFIDENCE_COLOR
        return QColor(
            int(low.red() + (high.red() - low.red()) * t),
            int(low.green() + (high.green() - low.green()) * t),
            int(low.blue() + (high.blue() - low.blue()) * t),
            190,
        )

    def paintEvent(self, event) -> None:  # noqa: N802
        """Paint a bordered plot area with one point per ranked drug prediction."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor(WINDOW_BACKGROUND))

        margin_left, margin_right = 10, 10
        margin_top, margin_bottom = 10, 18
        width = max(1, self.width() - margin_left - margin_right)
        height = max(1, self.height() - margin_top - margin_bottom)

        painter.setPen(QPen(QColor(SURFACE_CONTAINER), 1))
        painter.drawRect(margin_left, margin_top, width, height)

        if not self._points:
            painter.setPen(QColor(TEXT_MUTED))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Waiting for results...")
            return

        ic50_values = [ic50 for ic50, _confidence in self._points]
        log_ic50 = [math.log10(max(ic50, 1e-6)) for ic50 in ic50_values]
        lo, hi = min(log_ic50), max(log_ic50)
        span = (hi - lo) or 1.0

        painter.setPen(Qt.PenStyle.NoPen)
        for (_ic50, confidence), log_value in zip(self._points, log_ic50):
            px = margin_left + int(((log_value - lo) / span) * width)
            py = margin_top + int((1.0 - (confidence / 100.0)) * height)
            painter.setBrush(self._confidence_color(confidence))
            painter.drawEllipse(QPointF(px, py), 3.6, 3.6)

        # Y-axis (confidence %) ticks.
        painter.setPen(QColor(TEXT_MUTED))
        painter.drawText(margin_left + 4, margin_top + 12, "100%")
        painter.drawText(margin_left + 4, margin_top + height - 4, "0%")

        # X-axis (predicted IC50, uM) ticks -- real min/max of this run's points.
        lo_text = _format_ic50(min(ic50_values))
        painter.drawText(margin_left, margin_top + height + 14, lo_text)
        hi_text = _format_ic50(max(ic50_values))
        hi_text_width = painter.fontMetrics().horizontalAdvance(hi_text)
        painter.drawText(margin_left + width - hi_text_width, margin_top + height + 14, hi_text)


class FinalResultsPage(QWidget):
    """Displays the model's final predicted drug rankings and patient profile."""

    def __init__(
        self,
        parent: QWidget | None = None,
        on_upload_clicked: Callable[[], None] | None = None,
        on_model_running_clicked: Callable[[], None] | None = None,
    ) -> None:
        """Build the page.

        Args:
            parent: Optional Qt parent widget.
            on_upload_clicked: Invoked when the header's "Upload" tab is
                clicked; should navigate back to the dataset upload page.
            on_model_running_clicked: Invoked when the header's "Model
                Running" tab is clicked; should navigate to the model
                execution log page.
        """
        super().__init__(parent)
        self.setObjectName("FinalResultsPage")
        self._on_upload_clicked = on_upload_clicked
        self._on_model_running_clicked = on_model_running_clicked

        self._results_worker: ResultsWorker | None = None
        self._training_history_worker: TrainingHistoryWorker | None = None
        self._table: QTableWidget | None = None
        self._drug_name_label: QLabel | None = None
        self._drug_ic50_value: QLabel | None = None
        self._drug_confidence_value: QLabel | None = None
        self._drug_rank_slot: QVBoxLayout | None = None
        self._sample_subtitle_label: QLabel | None = None
        self._training_curve_widget: TrainingCurveWidget | None = None
        self._scatter_widget: PredictionScatterWidget | None = None
        self._all_rows: list[DrugResult] = []
        self._search_input: QLineEdit | None = None
        self._rank_filter_combo: QComboBox | None = None
        self._sort_combo: QComboBox | None = None

        self._build_ui()

    # -- live results -----------------------------------------------------

    def load_results(self, run_id: str) -> None:
        """Fetch and display the ranked drug predictions and training curve for `run_id`."""
        self._results_worker = ResultsWorker(run_id, parent=self)
        self._results_worker.succeeded.connect(self._on_results_succeeded)
        self._results_worker.failed.connect(self._on_results_failed)
        self._results_worker.start()

        self._training_history_worker = TrainingHistoryWorker(run_id, parent=self)
        self._training_history_worker.succeeded.connect(self._on_training_history_succeeded)
        self._training_history_worker.start()

    def _on_results_succeeded(self, raw_results: list[dict]) -> None:
        rows = [
            DrugResult(
                name=item["drug_name"],
                ic50=item["predicted_ic50_um"],
                rank=item["ranking"],
                confidence=item["confidence_percent"],
            )
            for item in raw_results
        ]
        rows.sort(key=lambda row: row.ic50)
        self._all_rows = rows
        self._populate_drug_details(rows[0] if rows else None)
        if self._scatter_widget is not None:
            self._scatter_widget.set_points([(row.ic50, row.confidence) for row in rows])
        self._apply_filters()

    def _on_results_failed(self, message: str) -> None:
        QMessageBox.critical(self, "Could Not Load Results", message)

    def _on_training_history_succeeded(self, history: list[dict]) -> None:
        # No failure handler wired -- the training curve is a supplementary
        # panel, not worth interrupting the user with an error dialog if it
        # can't be fetched; it just keeps showing its "waiting" state.
        if self._training_curve_widget is not None:
            self._training_curve_widget.set_history(history)

    def _build_ui(self) -> None:
        """Lay out the sidebar, header, and scrollable body content."""
        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        root.addWidget(self._build_sidebar(), 0)

        shell = QFrame()
        shell.setObjectName("RootShell")
        shell_layout = QVBoxLayout(shell)
        shell_layout.setContentsMargins(0, 0, 0, 0)
        shell_layout.setSpacing(0)

        shell_layout.addWidget(self._build_header(), 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        body = QWidget()
        body_layout = QVBoxLayout(body)
        body_layout.setContentsMargins(PAGE_MARGIN, PAGE_MARGIN, PAGE_MARGIN, PAGE_MARGIN)
        body_layout.setSpacing(SECTION_SPACING)

        body_layout.addLayout(self._build_page_header())
        body_layout.addLayout(self._build_panels_grid())
        body_layout.addStretch(1)

        scroll.setWidget(body)
        shell_layout.addWidget(scroll, 1)
        root.addWidget(shell, 1)

    def _build_sidebar(self) -> QFrame:
        """Build the shared sidebar (no model nav section on this page)."""
        footer_widgets = [
            make_sidebar_nav_button("History", "history"),
            make_sidebar_nav_button("Settings", "settings"),
            make_sidebar_nav_button("Support", "help_outline"),
        ]
        return build_sidebar(footer_widgets=footer_widgets, cta_widget=make_primary_cta_button())

    def _build_header(self) -> QFrame:
        """Build the shared header bar with "Results" as the active tab."""
        tabs = [
            make_top_tab("Upload", callback=self._on_upload_clicked),
            make_top_tab("Model Running", callback=self._on_model_running_clicked),
            make_top_tab("Results", active=True),
        ]
        return build_header_bar(tabs)

    def _build_page_header(self) -> QHBoxLayout:
        """Build the title/subtitle row and the download/export action buttons."""
        title_row = QHBoxLayout()
        title_col = QVBoxLayout()
        title = QLabel("Results")
        title.setStyleSheet(PAGE_TITLE_STYLE)
        subtitle = QLabel("Predicted Drug Results and Patient Molecular Profile")
        subtitle.setStyleSheet(PAGE_SUBTITLE_STYLE)
        title_col.addWidget(title)
        title_col.addWidget(subtitle)
        title_row.addLayout(title_col)
        title_row.addStretch(1)

        actions = QHBoxLayout()
        actions.setSpacing(12)
        download_btn = QPushButton(f"{icon_text('download')}  Download Data")
        download_btn.setStyleSheet(SECONDARY_BUTTON_STYLE)
        export_btn = QPushButton(f"{icon_text('summarize')}  Export Clinical Report")
        export_btn.setStyleSheet(PRIMARY_BUTTON_STYLE)
        actions.addWidget(download_btn)
        actions.addWidget(export_btn)
        title_row.addLayout(actions)
        return title_row

    def _build_panels_grid(self) -> QGridLayout:
        """Build the two-column grid of result panels.

        Left column: predicted drug results table plus the IC50/scatter
        charts. Right column: top-drug detail card and molecular profile
        summary. Column widths are weighted 2:1 to give the data-dense left
        column more room.
        """
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(SECTION_SPACING)
        grid.setVerticalSpacing(SECTION_SPACING)

        left_col = QVBoxLayout()
        left_col.setSpacing(SECTION_SPACING)
        left_col.addWidget(self._build_predicted_results_panel())

        visuals = QGridLayout()
        visuals.setContentsMargins(0, 0, 0, 0)
        visuals.setHorizontalSpacing(SECTION_SPACING)
        visuals.addWidget(self._build_ic50_panel(), 0, 0)
        visuals.addWidget(self._build_scatter_panel(), 0, 1)
        left_col.addLayout(visuals)

        right_col = QVBoxLayout()
        right_col.setSpacing(SECTION_SPACING)
        right_col.addWidget(self._build_drug_details_panel())
        right_col.addWidget(self._build_sample_profile_panel())
        right_col.addStretch(1)

        left_host = QWidget()
        left_host.setLayout(left_col)
        right_host = QWidget()
        right_host.setLayout(right_col)

        grid.addWidget(left_host, 0, 0)
        grid.addWidget(right_host, 0, 1)
        grid.setColumnStretch(0, 2)
        grid.setColumnStretch(1, 1)
        return grid

    def _build_predicted_results_panel(self) -> QFrame:
        """Build the "Predicted Drug Results Panel" card and its table.

        The table uses the same shared header/border/row styling
        (`widgets.tables.style_data_table`) and the same status-badge and
        progress-bar cell widgets as the model execution log page's
        pipeline table, so both tables clearly belong to one design system.
        """
        card = QFrame()
        card.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        header = QFrame()
        header.setStyleSheet(f"background: {WINDOW_BACKGROUND}; border-bottom: 1px solid {SURFACE_CONTAINER};")
        header_row = QHBoxLayout(header)
        header_row.setContentsMargins(16, 12, 16, 12)
        header_title = QLabel("Predicted Drug Results Panel")
        header_title.setStyleSheet(f"background: transparent; border: none; {CARD_TITLE_STYLE}")
        info_icon = QLabel(icon_text("info"))
        info_icon.setStyleSheet(f"background: transparent; border: none; font-size: 18px; color: {TEXT_MUTED};")
        header_row.addWidget(header_title)
        header_row.addStretch(1)
        header_row.addWidget(info_icon)
        layout.addWidget(header)

        controls = QFrame()
        controls.setStyleSheet(f"background: {SURFACE}; border-bottom: 1px solid {SURFACE_CONTAINER};")
        controls_row = QHBoxLayout(controls)
        controls_row.setContentsMargins(16, 10, 16, 10)
        controls_row.setSpacing(8)

        search_input = QLineEdit()
        search_input.setPlaceholderText("Search drug name...")
        search_input.textChanged.connect(self._apply_filters)
        self._search_input = search_input

        rank_filter_combo = QComboBox()
        rank_filter_combo.addItems(["All Rankings", "HIGH SENSITIVITY", "MEDIUM", "LOW"])
        rank_filter_combo.currentTextChanged.connect(self._apply_filters)
        self._rank_filter_combo = rank_filter_combo

        sort_combo = QComboBox()
        sort_combo.addItems(
            [
                "IC50 (Low → High)",
                "IC50 (High → Low)",
                "Confidence (High → Low)",
                "Confidence (Low → High)",
                "Drug Name (A → Z)",
            ]
        )
        sort_combo.currentTextChanged.connect(self._apply_filters)
        self._sort_combo = sort_combo

        controls_row.addWidget(search_input, 1)
        controls_row.addWidget(rank_filter_combo)
        controls_row.addWidget(sort_combo)
        layout.addWidget(controls)

        table = QTableWidget(0, 4)
        table.setHorizontalHeaderLabels(["Drug Name", "Predicted IC50 (uM)", "Sensitivity Ranking", "Confidence"])
        style_data_table(table, header_background=WINDOW_BACKGROUND)
        table.setColumnWidth(0, 140)
        table.setColumnWidth(1, 150)
        table.setColumnWidth(2, 170)
        self._table = table
        self._populate_table([])

        layout.addWidget(table)
        return card

    def _populate_table(self, rows: list[DrugResult], empty_message: str = "Waiting for a completed model run...") -> None:
        """(Re)fill the predicted-results table from real `DrugResult` rows.

        Called with an empty list while a run's results haven't loaded yet
        (shows a single `empty_message` placeholder row, default "waiting"),
        and again -- via `_apply_filters` -- once `load_results` fetches the
        real ranked predictions or the user changes the search/filter/sort
        controls (with `empty_message` swapped to a "no matches" message
        when a filter excludes every row).
        """
        table = self._table
        if table is None:
            return

        if not rows:
            table.setRowCount(1)
            waiting_item = QTableWidgetItem(empty_message)
            waiting_item.setForeground(QColor(TEXT_MUTED))
            table.setItem(0, 0, waiting_item)
            table.setSpan(0, 0, 1, 4)
            return

        table.clearSpans()
        table.setRowCount(len(rows))
        for row_index, row in enumerate(rows):
            name_item = QTableWidgetItem(row.name)
            name_item.setForeground(QColor(TEXT))
            table.setItem(row_index, 0, name_item)

            ic50_item = QTableWidgetItem(_format_ic50(row.ic50))
            ic50_item.setForeground(QColor(TEXT_MUTED))
            table.setItem(row_index, 1, ic50_item)

            tone = _RANK_TO_BADGE_TONE.get(row.rank, "neutral")
            table.setCellWidget(row_index, 2, build_status_badge(row.rank, tone))
            table.setCellWidget(row_index, 3, self._build_confidence_cell(row.confidence))

    # Sort dropdown label -> (sort key, reverse). Keeps `_apply_filters`
    # itself free of a long if/elif chain.
    _SORT_OPTIONS: dict[str, tuple[Callable[["DrugResult"], object], bool]] = {
        "IC50 (Low → High)": (lambda row: row.ic50, False),
        "IC50 (High → Low)": (lambda row: row.ic50, True),
        "Confidence (High → Low)": (lambda row: row.confidence, True),
        "Confidence (Low → High)": (lambda row: row.confidence, False),
        "Drug Name (A → Z)": (lambda row: row.name.lower(), False),
    }

    def _apply_filters(self) -> None:
        """Re-derive the table's rows from `self._all_rows` per the current
        search/filter/sort controls, leaving `self._all_rows` itself (and
        anything else driven from it, like the top-drug and scatter panels)
        untouched.
        """
        if self._search_input is None or self._rank_filter_combo is None or self._sort_combo is None:
            return

        search_text = self._search_input.text().strip().lower()
        rank_filter = self._rank_filter_combo.currentText()

        filtered = [
            row
            for row in self._all_rows
            if (not search_text or search_text in row.name.lower())
            and (rank_filter == "All Rankings" or row.rank == rank_filter)
        ]

        sort_key, reverse = self._SORT_OPTIONS.get(self._sort_combo.currentText(), (lambda row: row.ic50, False))
        filtered.sort(key=sort_key, reverse=reverse)

        empty_message = "No drugs match your search/filter." if self._all_rows else "Waiting for a completed model run..."
        self._populate_table(filtered, empty_message=empty_message)

    @staticmethod
    def _build_confidence_cell(confidence: float) -> QWidget:
        """Build the "Confidence" cell: a percentage label over a mini progress bar.

        Reuses `widgets.tables.build_mini_progress_bar` — the same bar used
        for pipeline completion on the model execution log page — so the
        confidence measurement reads visually consistent with the rest of
        the app's progress indicators.

        Args:
            confidence: Confidence percentage (0-100).

        Returns:
            A `QWidget` suitable for `QTableWidget.setCellWidget`.
        """
        wrapper = transparent_cell_widget()
        layout = QVBoxLayout(wrapper)
        layout.setContentsMargins(0, 6, 0, 6)
        layout.setSpacing(6)

        label = QLabel(f"{confidence:.1f}%")
        label.setStyleSheet(
            f"background: transparent; border: none; font-family: Consolas, monospace;"
            f" font-size: 13px; font-weight: 700; color: {TEXT};"
        )
        layout.addWidget(label)
        layout.addWidget(build_mini_progress_bar(int(round(confidence))))
        return wrapper

    def _build_ic50_panel(self) -> QFrame:
        """Build the "Training Curve" chart card (real train_loss/val_rmse per epoch)."""
        self._training_curve_widget = TrainingCurveWidget()
        return self._build_chart_panel("Training Curve (Train Loss / Val RMSE)", self._training_curve_widget)

    def _build_scatter_panel(self) -> QFrame:
        """Build the "Predicted IC50 vs. Confidence" chart card."""
        self._scatter_widget = PredictionScatterWidget()
        return self._build_chart_panel("Predicted IC50 vs. Confidence", self._scatter_widget)

    @staticmethod
    def _build_chart_panel(title_text: str, chart_widget: QWidget) -> QFrame:
        """Build a small card containing a caption and a dashed-border chart frame.

        Args:
            title_text: Caption shown above the chart.
            chart_widget: The chart widget to embed (already sized via its
                own `minimumHeight`).

        Returns:
            A styled `QFrame` card.
        """
        card = QFrame()
        card.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)

        title = QLabel(title_text)
        title.setStyleSheet(LABEL_CAPS_STYLE)
        layout.addWidget(title)

        chart_frame = QFrame()
        chart_frame.setStyleSheet(f"background: {WINDOW_BACKGROUND}; border: 1px dashed {BORDER}; border-radius: 8px;")
        chart_layout = QVBoxLayout(chart_frame)
        chart_layout.setContentsMargins(8, 8, 8, 8)
        chart_layout.addWidget(chart_widget)
        layout.addWidget(chart_frame, 1)
        return card

    def _build_drug_details_panel(self) -> QFrame:
        """Build the top-ranked drug's detail card (name, predicted IC50, ranking, confidence).

        Only shows values the model actually produces — no fabricated
        chemistry (formula/weight/description); the repo has no drug
        metadata source to back those with real data.
        """
        card = QFrame()
        card.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        top = QHBoxLayout()
        title_col = QVBoxLayout()
        drug_name = QLabel("—")
        drug_name.setStyleSheet(label_style(f"font-size: 24px; font-weight: 600; color: {TEXT};"))
        self._drug_name_label = drug_name
        subtitle = QLabel("Top Predicted Match")
        subtitle.setStyleSheet(LABEL_CAPS_STYLE)
        title_col.addWidget(drug_name)
        title_col.addWidget(subtitle)
        top.addLayout(title_col)
        top.addStretch(1)
        layout.addLayout(top)

        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setStyleSheet(f"color: {BORDER};")
        layout.addWidget(separator)

        metadata = QGridLayout()
        metadata.setHorizontalSpacing(16)
        metadata.setVerticalSpacing(4)
        ic50_label = QLabel("Predicted IC50 (uM)")
        confidence_label = QLabel("Confidence")
        ic50_value = QLabel("—")
        confidence_value = QLabel("—")
        for label in (ic50_label, confidence_label):
            label.setStyleSheet(LABEL_CAPS_STYLE)
        for value in (ic50_value, confidence_value):
            value.setStyleSheet(label_style(f"font-family: Consolas, monospace; font-size: 16px; color: {TEXT};"))
        self._drug_ic50_value = ic50_value
        self._drug_confidence_value = confidence_value
        metadata.addWidget(ic50_label, 0, 0)
        metadata.addWidget(confidence_label, 0, 1)
        metadata.addWidget(ic50_value, 1, 0)
        metadata.addWidget(confidence_value, 1, 1)
        layout.addLayout(metadata)

        rank_slot = QVBoxLayout()
        rank_slot.setContentsMargins(0, 4, 0, 0)
        self._drug_rank_slot = rank_slot
        layout.addLayout(rank_slot)

        return card

    def _populate_drug_details(self, top_result: DrugResult | None) -> None:
        """Fill the drug details card from the top-ranked real prediction."""
        if top_result is None:
            return
        if self._drug_name_label is not None:
            self._drug_name_label.setText(top_result.name)
        if self._drug_ic50_value is not None:
            self._drug_ic50_value.setText(_format_ic50(top_result.ic50))
        if self._drug_confidence_value is not None:
            self._drug_confidence_value.setText(f"{top_result.confidence:.1f}%")
        if self._drug_rank_slot is not None:
            while self._drug_rank_slot.count():
                item = self._drug_rank_slot.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
            tone = _RANK_TO_BADGE_TONE.get(top_result.rank, "neutral")
            self._drug_rank_slot.addWidget(build_status_badge(top_result.rank, tone))

    def _build_sample_profile_panel(self) -> QFrame:
        """Build the sample profile card: target cell line + omics modalities used.

        Replaces the old "Molecular Profile" card, which fabricated a
        patient ID and gene mutations that don't exist anywhere in this
        cell-line-based pipeline.
        """
        card = QFrame()
        card.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        title = QLabel("Sample Profile")
        title.setStyleSheet(CARD_TITLE_STYLE)
        subtitle = QLabel("Cell Line: —")
        subtitle.setStyleSheet(label_style(f"font-size: 14px; color: {TEXT_MUTED};"))
        self._sample_subtitle_label = subtitle
        layout.addWidget(title)
        layout.addWidget(subtitle)

        modalities_title = QLabel("Omics Modalities Used")
        modalities_title.setStyleSheet(LABEL_CAPS_STYLE)
        layout.addWidget(modalities_title)

        chips = QHBoxLayout()
        chips.setSpacing(8)
        chips.addWidget(self._chip("Transcriptomics", primary=True))
        chips.addWidget(self._chip("Genomics (Mut/CNV)"))
        chips.addWidget(self._chip("Proteomics"))
        chips.addStretch(1)
        layout.addLayout(chips)

        network_title = QLabel("Network Proximity")
        network_title.setStyleSheet(LABEL_CAPS_STYLE)
        layout.addWidget(network_title)
        layout.addWidget(self._build_network_proximity_placeholder())

        return card

    def set_sample_id(self, sanger_model_id: str) -> None:
        """Update the Sample Profile card's cell-line identifier."""
        if self._sample_subtitle_label is not None:
            self._sample_subtitle_label.setText(f"Cell Line: {sanger_model_id}")

    @staticmethod
    def _build_network_proximity_placeholder() -> QFrame:
        """Build the dashed placeholder box for the network proximity visualization."""
        card = QFrame()
        card.setFixedHeight(120)
        card.setStyleSheet(f"background: {WINDOW_BACKGROUND}; border: 1px dashed {BORDER}; border-radius: 8px;")
        layout = QVBoxLayout(card)
        icon = QLabel(icon_text("hub"))
        icon.setAlignment(Qt.AlignmentFlag.AlignCenter)
        icon.setStyleSheet(label_style(f"font-size: 24px; color: {TEXT_MUTED};"))
        layout.addWidget(icon)
        return card

    @staticmethod
    def _chip(text: str, primary: bool = False) -> QLabel:
        """Build a small rounded label chip (used for mutation tags).

        Args:
            text: Chip label.
            primary: Whether to render as a filled (primary) chip versus an
                outlined (secondary) chip.

        Returns:
            A styled `QLabel`.
        """
        chip = QLabel(text)
        if primary:
            variant_style = f"background: {PRIMARY}; color: {SURFACE}; border: none;"
        else:
            variant_style = f"background: {SURFACE_HIGH}; color: {TEXT}; border: 1px solid {BORDER};"
        chip.setStyleSheet(
            "padding: 3px 8px; border-radius: 3px; font-family: Consolas, monospace; font-size: 12px;" + variant_style
        )
        return chip
