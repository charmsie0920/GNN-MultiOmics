"""Final Results page.

The last stop in the workflow: shows the model's predicted drug rankings for
the current patient alongside supporting visualizations (IC50 distribution,
biomarker correlation scatter) and detail panels (top drug details, patient
molecular profile).

Main UI components:
    - Shared sidebar (`widgets.navigation.build_sidebar`) — this page has no
      "Model Visualization"/"Model Logs" nav section, matching the upload
      page's sidebar, since results is a terminal step in the workflow.
    - Shared header bar (`widgets.navigation.build_header_bar`) with
      "Results" marked as the active workflow tab.
    - A "Predicted Drug Results Panel" table, styled with the same shared
      table helpers (`widgets.tables`) as the model execution log page's
      pipeline table, so the two read as one design system.
    - Two custom-painted placeholder charts (`IC50DistributionWidget`,
      `ScatterPlaceholderWidget`).
    - Drug detail and molecular profile summary cards.

Interactions with other pages:
    - `on_upload_clicked` navigates back to `DatasetInitializationPage`.
    - `on_model_running_clicked` navigates to `ModelExecutionLogPage`.
    Both callbacks are supplied and wired by `main.py`.

The entire results payload (drug rankings, chart data, patient profile) is
currently mocked and must be replaced by real backend inference output.
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
    QPushButton,
    QScrollArea,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

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
        # * Entire results payload is currently mocked and must be replaced by backend inference outputs.
        self._on_upload_clicked = on_upload_clicked
        self._on_model_running_clicked = on_model_running_clicked
        self._build_ui()

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
        right_col.addWidget(self._build_molecular_profile_panel())
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

        rows = [
            # * Replace with backend response from /api/v1/results/drug-ranking.
            DrugResult("Afatinib", 15.2, "HIGH SENSITIVITY", 98.4),
            DrugResult("Erlotinib", 25.0, "HIGH SENSITIVITY", 95.1),
            DrugResult("Gefitinib", 32.8, "HIGH SENSITIVITY", 92.0),
            DrugResult("Osimertinib", 58.5, "MEDIUM", 88.2),
            DrugResult("Dacomitinib", 85.0, "LOW", 81.5),
        ]

        table = QTableWidget(len(rows), 4)
        table.setHorizontalHeaderLabels(["Drug Name", "Predicted IC50 (uM)", "Sensitivity Ranking", "Confidence"])
        style_data_table(table, header_background=WINDOW_BACKGROUND)

        for row_index, row in enumerate(rows):
            name_item = QTableWidgetItem(row.name)
            name_item.setForeground(QColor(TEXT))
            table.setItem(row_index, 0, name_item)

            ic50_item = QTableWidgetItem(f"{row.ic50:.1f}")
            ic50_item.setForeground(QColor(TEXT_MUTED))
            table.setItem(row_index, 1, ic50_item)

            tone = _RANK_TO_BADGE_TONE.get(row.rank, "neutral")
            table.setCellWidget(row_index, 2, build_status_badge(row.rank, tone))
            table.setCellWidget(row_index, 3, self._build_confidence_cell(row.confidence))

        table.setColumnWidth(0, 140)
        table.setColumnWidth(1, 150)
        table.setColumnWidth(2, 170)
        layout.addWidget(table)
        return card

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
        """Build the top-ranked drug's detail card (name, structure, description, metadata)."""
        card = QFrame()
        card.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        top = QHBoxLayout()
        title_col = QVBoxLayout()
        drug_name = QLabel("Afatinib")
        drug_name.setStyleSheet(label_style(f"font-size: 24px; font-weight: 600; color: {TEXT};"))
        drug_type = QLabel("Tyrosine Kinase Inhibitor")
        drug_type.setStyleSheet(LABEL_CAPS_STYLE)
        title_col.addWidget(drug_name)
        title_col.addWidget(drug_type)
        top.addLayout(title_col)
        top.addStretch(1)

        open_btn = QPushButton(icon_text("open_in_new"))
        open_btn.setFixedSize(32, 32)
        open_btn.setStyleSheet(
            f"border: 1px solid {BORDER}; border-radius: 4px; background: {SURFACE}; color: {TEXT_MUTED};"
        )
        top.addWidget(open_btn)
        layout.addLayout(top)

        molecule_frame = QFrame()
        molecule_frame.setFixedHeight(180)
        molecule_frame.setStyleSheet(
            "background: qlineargradient(x1:0,y1:0,x2:1,y2:1, stop:0 #f3f3f3, stop:1 #e8e8e8);"
            f" border: 1px solid {BORDER}; border-radius: 6px;"
        )
        molecule_layout = QVBoxLayout(molecule_frame)
        molecule_layout.setContentsMargins(0, 0, 0, 0)
        molecule_label = QLabel("C24H25ClFN5O3")
        molecule_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        molecule_label.setStyleSheet(label_style(f"font-family: Consolas, monospace; font-size: 14px; color: {TEXT};"))
        molecule_layout.addWidget(molecule_label)
        layout.addWidget(molecule_frame)

        description_title = QLabel("Clinical Description")
        description_title.setStyleSheet(LABEL_CAPS_STYLE)
        description = QLabel(
            "An irreversible kinase inhibitor that selectively targets EGFR and HER2 with therapeutic"
            " potential in solid tumors."
        )
        description.setWordWrap(True)
        description.setStyleSheet(label_style(f"font-size: 14px; color: {TEXT};"))
        layout.addWidget(description_title)
        layout.addWidget(description)

        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setStyleSheet(f"color: {BORDER};")
        layout.addWidget(separator)

        metadata = QGridLayout()
        metadata.setHorizontalSpacing(16)
        composition_label = QLabel("Composition")
        composition_value = QLabel("C24H25ClFN5O3")
        weight_label = QLabel("Weight")
        weight_value = QLabel("485.94 g/mol")
        for label in (composition_label, weight_label):
            label.setStyleSheet(LABEL_CAPS_STYLE)
        for value in (composition_value, weight_value):
            value.setStyleSheet(label_style(f"font-family: Consolas, monospace; font-size: 13px; color: {TEXT};"))
        metadata.addWidget(composition_label, 0, 0)
        metadata.addWidget(weight_label, 0, 1)
        metadata.addWidget(composition_value, 1, 0)
        metadata.addWidget(weight_value, 1, 1)
        layout.addLayout(metadata)
        return card

    def _build_molecular_profile_panel(self) -> QFrame:
        """Build the patient molecular profile card (mutations, expression, network proximity)."""
        card = QFrame()
        card.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        title = QLabel("Molecular Profile")
        title.setStyleSheet(CARD_TITLE_STYLE)
        subtitle = QLabel("Patient: PT-8842-AX")
        # * Replace patient identifier and molecular details with backend-provided patient context.
        subtitle.setStyleSheet(label_style(f"font-size: 14px; color: {TEXT_MUTED};"))
        layout.addWidget(title)
        layout.addWidget(subtitle)

        mutations_title = QLabel("Key Mutations Identified")
        mutations_title.setStyleSheet(LABEL_CAPS_STYLE)
        layout.addWidget(mutations_title)

        chips = QHBoxLayout()
        chips.setSpacing(8)
        chips.addWidget(self._chip("EGFR L858R", primary=True))
        chips.addWidget(self._chip("TP53 R175H"))
        chips.addWidget(self._chip("PIK3CA E545K"))
        chips.addStretch(1)
        layout.addLayout(chips)

        expression_title = QLabel("Expression Summary")
        expression_title.setStyleSheet(LABEL_CAPS_STYLE)
        layout.addWidget(expression_title)
        layout.addWidget(self._build_expression_summary())

        network_title = QLabel("Network Proximity")
        network_title.setStyleSheet(LABEL_CAPS_STYLE)
        layout.addWidget(network_title)
        layout.addWidget(self._build_network_proximity_placeholder())

        return card

    @staticmethod
    def _build_expression_summary() -> QFrame:
        """Build the "EGFR Overexpression" row with a value label and progress bar."""
        card = QFrame()
        card.setStyleSheet(f"background: {WINDOW_BACKGROUND}; border: 1px solid {BORDER}; border-radius: 6px;")
        layout = QVBoxLayout(card)
        layout.setContentsMargins(12, 10, 12, 10)
        layout.setSpacing(6)

        row = QHBoxLayout()
        overexpression_label = QLabel("EGFR Overexpression")
        overexpression_label.setStyleSheet(label_style(f"font-size: 14px; color: {TEXT};"))
        row.addWidget(overexpression_label)
        value = QLabel("High (98th pct)")
        value.setStyleSheet(label_style(f"font-family: Consolas, monospace; font-weight: 700; color: {TEXT};"))
        row.addStretch(1)
        row.addWidget(value)
        layout.addLayout(row)

        bar = build_mini_progress_bar(98, height=6)
        layout.addWidget(bar)
        return card

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
