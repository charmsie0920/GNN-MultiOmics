"""Final Results page.

The last stop in the workflow: shows the model's predicted drug rankings for
the selected target cell line, alongside supporting visualizations (IC50
distribution, biomarker correlation scatter — still mocked) and detail
panels (top drug details, sample profile).

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
    - Two custom-painted placeholder charts (`IC50DistributionWidget`,
      `ScatterPlaceholderWidget`) — still mocked; out of scope for now.
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

from collections.abc import Callable
from dataclasses import dataclass

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from client.workers import ResultsWorker
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


class IC50DistributionWidget(QWidget):
    """Custom-painted placeholder line chart for the IC50 distribution panel."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(185)

    def paintEvent(self, event) -> None:  # noqa: N802
        """Paint a dashed gridline background and a shaded distribution curve."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor(WINDOW_BACKGROUND))

        margin = 12
        width = max(1, self.width() - (margin * 2))
        height = max(1, self.height() - (margin * 2))

        grid_pen = QPen(QColor(SURFACE_CONTAINER))
        grid_pen.setStyle(Qt.PenStyle.DashLine)
        painter.setPen(grid_pen)
        for ratio in (0.25, 0.5, 0.75):
            y = margin + int(height * ratio)
            painter.drawLine(margin, y, margin + width, y)

        # Normalized (0-1, 0-1) placeholder curve points; replace with real
        # IC50 distribution samples once backend results are available.
        points = [
            (0.00, 0.90),
            (0.10, 0.88),
            (0.22, 0.42),
            (0.36, 0.30),
            (0.53, 0.20),
            (0.66, 0.72),
            (0.79, 0.54),
            (0.90, 0.66),
            (1.00, 0.60),
        ]

        def map_xy(px: float, py: float) -> QPointF:
            return QPointF(margin + (px * width), margin + (py * height))

        stroke_path = QPainterPath()
        area_path = QPainterPath()
        first = map_xy(points[0][0], points[0][1])
        stroke_path.moveTo(first)
        area_path.moveTo(margin, margin + height)
        area_path.lineTo(first)

        for x, y in points[1:]:
            pt = map_xy(x, y)
            stroke_path.lineTo(pt)
            area_path.lineTo(pt)

        area_path.lineTo(margin + width, margin + height)
        area_path.closeSubpath()
        painter.fillPath(area_path, QColor(116, 120, 120, 60))

        line_pen = QPen(QColor("#1c1b1b"))
        line_pen.setWidth(2)
        painter.setPen(line_pen)
        painter.drawPath(stroke_path)


class ScatterPlaceholderWidget(QWidget):
    """Custom-painted placeholder scatter chart for the biomarker correlation panel."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(185)

    def paintEvent(self, event) -> None:  # noqa: N802
        """Paint a bordered plot area with scattered points and a highlighted cluster."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor(WINDOW_BACKGROUND))

        margin = 14
        width = max(1, self.width() - (margin * 2))
        height = max(1, self.height() - (margin * 2))

        painter.setPen(QPen(QColor(SURFACE_CONTAINER), 1))
        painter.drawRect(margin, margin, width, height)

        # Normalized (0-1, 0-1) placeholder scatter points; replace with real
        # biomarker correlation samples once backend results are available.
        points = [
            (0.10, 0.80),
            (0.18, 0.70),
            (0.24, 0.65),
            (0.31, 0.45),
            (0.38, 0.55),
            (0.47, 0.35),
            (0.54, 0.32),
            (0.63, 0.23),
            (0.72, 0.18),
            (0.81, 0.28),
            (0.88, 0.12),
        ]
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(68, 71, 72, 160))
        for x, y in points:
            px = margin + int(x * width)
            py = margin + int(y * height)
            painter.drawEllipse(QPointF(px, py), 3.2, 3.2)

        painter.setBrush(QColor(26, 28, 28, 55))
        painter.drawEllipse(QPointF(margin + (width / 2), margin + (height / 2)), 22, 22)


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
        self._table: QTableWidget | None = None
        self._drug_name_label: QLabel | None = None
        self._drug_ic50_value: QLabel | None = None
        self._drug_confidence_value: QLabel | None = None
        self._drug_rank_slot: QVBoxLayout | None = None
        self._sample_subtitle_label: QLabel | None = None

        self._build_ui()

    # -- live results -----------------------------------------------------

    def load_results(self, run_id: str) -> None:
        """Fetch and display the ranked drug predictions for `run_id`."""
        self._results_worker = ResultsWorker(run_id, parent=self)
        self._results_worker.succeeded.connect(self._on_results_succeeded)
        self._results_worker.failed.connect(self._on_results_failed)
        self._results_worker.start()

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
        self._populate_table(rows)
        self._populate_drug_details(rows[0] if rows else None)

    def _on_results_failed(self, message: str) -> None:
        QMessageBox.critical(self, "Could Not Load Results", message)

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

    def _populate_table(self, rows: list[DrugResult]) -> None:
        """(Re)fill the predicted-results table from real `DrugResult` rows.

        Called with an empty list while a run's results haven't loaded yet
        (shows a single "waiting" placeholder row), and again once
        `load_results` fetches the real ranked predictions.
        """
        table = self._table
        if table is None:
            return

        if not rows:
            table.setRowCount(1)
            waiting_item = QTableWidgetItem("Waiting for a completed model run...")
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
        """Build the "IC50 Distribution Comparison" chart card."""
        return self._build_chart_panel("IC50 Distribution Comparison", IC50DistributionWidget())

    def _build_scatter_panel(self) -> QFrame:
        """Build the "Biomarker Correlation Scatter" chart card."""
        return self._build_chart_panel("Biomarker Correlation Scatter", ScatterPlaceholderWidget())

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
