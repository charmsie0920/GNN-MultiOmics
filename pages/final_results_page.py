"""Final Results page.

The last stop in the workflow: shows the model's predicted drug rankings for
the selected target cell line, and explains whichever drug the user picks.

Main UI components:
    - Shared sidebar (`widgets.navigation.build_sidebar`) with "Drug Results"
      (this page) and "Model Analytics" links. The charts live on
      `ModelAnalyticsPage`, keeping this page focused on the drugs.
    - Shared header bar (`widgets.navigation.build_header_bar`) with
      "Results" marked as the active workflow tab.
    - A sample profile strip: the cell line and omics modalities everything
      below is about.
    - A "Predicted Drug Results Panel" table, styled with the same shared
      table helpers (`widgets.tables`) as the model execution log page's
      pipeline table. Populated from real backend inference via
      `load_results(run_id)`.
    - A selected-drug pane beside the table: a dark header with the picked
      drug's prediction and annotation, and "Top genes" / "Pathways" views
      explaining that prediction. The table row and the pane are tied
      together visually (matching black marker + a crossfade on change)
      rather than with explanatory text.

Interactions with other pages:
    - `on_upload_clicked` navigates back to `DatasetInitializationPage`.
    - `on_model_running_clicked` navigates to `ModelExecutionLogPage`.
    - `on_analytics_clicked` navigates to `ModelAnalyticsPage`.
    All callbacks are supplied and wired by `UILauncher.py`, which also
    calls `load_results(run_id)` / `set_sample_id(target_cell_line)` before
    switching to this page.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime

from PySide6.QtCore import QPointF, QPropertyAnimation, QRect, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPalette
from PySide6.QtWidgets import (
    QButtonGroup,
    QComboBox,
    QFileDialog,
    QFrame,
    QGraphicsOpacityEffect,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QStackedWidget,
    QStyle,
    QStyledItemDelegate,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from client.workers import EnrichmentWorker, GeneAttributionWorker, ResultsWorker
from styles.theme import (
    BORDER,
    CARD_CONTAINER_STYLE,
    CONTRIBUTION_RESISTANCE,
    CONTRIBUTION_SENSITISING,
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
    SAMPLE_PROFILE_ACCENT,
    SAMPLE_PROFILE_ACCENT_WIDTH,
    SELECTED_DRUG_SURFACE,
    SELECTED_DRUG_TEXT,
    SELECTED_DRUG_TEXT_MUTED,
    SELECTED_DRUG_TINT,
    SELECTED_DRUG_TINT_STRONG,
    TEXT,
    TEXT_MUTED,
    WINDOW_BACKGROUND,
    style_card_section,
    label_style,
)
from widgets.formatting import format_ic50
from widgets.help import HELP, make_info_icon, set_header_help
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
# Shorter badge text for the table column, so "HIGH" lines up with "MEDIUM"
# and "LOW" in a narrow column. The pane keeps the full label.
_RANK_TABLE_LABELS = {"HIGH SENSITIVITY": "HIGH"}

# Ranking filter dropdown label -> the ranking it keeps (None = every drug).
_RANK_FILTER_OPTIONS: dict[str, str | None] = {
    "All rankings": None,
    "High": "HIGH SENSITIVITY",
    "Medium": "MEDIUM",
    "Low": "LOW",
}

# Widths (px) of the ranking filter and sort dropdowns beside the search box,
# each sized so its longest option just fits. Tweak here if the labels change.
_RANK_FILTER_WIDTH = 118
_SORT_WIDTH = 140

# Enrichr library id -> short label, so the column fits the narrow detail pane.
_LIBRARY_SHORT_NAMES = {
    "GO_Biological_Process_2023": "GO BP",
    "KEGG_2021_Human": "KEGG",
    "Reactome_2022": "Reactome",
}

# Segmented "Top genes | Pathways" toggle in the selected-drug pane.
_SEGMENT_BUTTON_STYLE = f"""
    QPushButton {{
        background: transparent; border: 1px solid transparent; border-radius: 6px;
        padding: 6px 14px; color: {TEXT_MUTED};
        font-size: 12px; font-weight: 700; letter-spacing: 0.05em; text-transform: uppercase;
    }}
    QPushButton:checked {{
        background: {SURFACE}; color: {SELECTED_DRUG_SURFACE}; border: 1px solid {SELECTED_DRUG_SURFACE};
    }}
    QPushButton:hover:!checked {{ color: {SELECTED_DRUG_SURFACE}; }}
"""

# Muted text on the selected-drug pane's grey header.
_ON_DARK_MUTED = SELECTED_DRUG_TEXT_MUTED


@dataclass(frozen=True)
class DrugResult:
    """One row of the predicted drug results table."""

    # Drug names repeat in GDSC (two rows both read "Dactinomycin"), so the id
    # is what identifies the row when asking the backend to explain it.
    drug_id: str
    name: str
    ic50: float
    rank: str
    confidence: float
    # GDSC annotation; free text, often a mechanism rather than a gene, and
    # sometimes empty.
    target: str = ""
    pathway: str = ""


class _SelectedRowMarkerDelegate(QStyledItemDelegate):
    """Drug-name column: the selected row gets a bold name and a black bar on its left.

    The bar is the same dark grey as the selected-drug pane's header beside
    the table, so the picked row and the pane read as one pair without any
    text saying so.
    """

    _BAR_WIDTH = 4

    def initStyleOption(self, option, index) -> None:  # noqa: N802
        super().initStyleOption(option, index)
        if option.state & QStyle.StateFlag.State_Selected:
            option.font.setBold(True)

    def paint(self, painter, option, index) -> None:
        super().paint(painter, option, index)
        if option.state & QStyle.StateFlag.State_Selected:
            rect = option.rect
            painter.fillRect(
                QRect(rect.left(), rect.top(), self._BAR_WIDTH, rect.height()), QColor(SELECTED_DRUG_SURFACE)
            )


class _ContributionBar(QWidget):
    """One gene's signed contribution as a diverging bar around a zero line.

    Negative scores (pushing the prediction toward sensitivity) extend left,
    positive ones (toward resistance) extend right, each in its own colour,
    with the exact value printed above the zero line. Direction is readable
    at a glance down the column, not just from a +/- sign.

    Length is `sqrt(|score| / max_abs)` of the half-width. Contributions span
    orders of magnitude -- a drug's own target routinely scores ~100x the rest,
    as the target edge is the drug node's only link into the graph -- so a
    linear scale leaves one full bar over a column of invisible slivers. The
    earlier log scale over-corrected: 0.8 and 0.01 drew at 100% and 60%,
    hiding real differences. Square root keeps the ranking and the gaps
    visible (0.01 vs 0.8 -> 11% vs 100%) while small genes stay drawable.
    """

    _BAR_HEIGHT = 8
    _BAR_RADIUS = 4
    _MIN_BAR = 2.0  # px, so a tiny non-zero contribution never vanishes
    _SIDE_MARGIN = 10

    def __init__(self, score: float, max_abs: float) -> None:
        super().__init__()
        self._score = score
        self._fraction = math.sqrt(abs(score) / max_abs) if max_abs > 0 else 0.0
        direction = "toward sensitivity" if score < 0 else "toward resistance"
        self.setToolTip(f"{score:+.4g} · {direction}")

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        centre_x = self.width() / 2
        half_width = max(1.0, centre_x - self._SIDE_MARGIN)
        bar_top = self.height() / 2 + 3

        # Exact value in text ink, centred over the zero line.
        font = QFont(painter.font())
        font.setFamily("Consolas")
        font.setPixelSize(12)
        painter.setFont(font)
        painter.setPen(QColor(TEXT))
        text = f"{self._score:+.3g}"
        text_width = painter.fontMetrics().horizontalAdvance(text)
        painter.drawText(QPointF(centre_x - text_width / 2, bar_top - 6), text)

        # Faint full-width track: the scale the bar is measured against.
        track = QRectF(self._SIDE_MARGIN, bar_top, 2 * half_width, self._BAR_HEIGHT)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(SURFACE_CONTAINER))
        painter.drawRoundedRect(track, self._BAR_RADIUS, self._BAR_RADIUS)

        if self._score != 0:
            length = max(self._MIN_BAR, self._fraction * half_width)
            if self._score < 0:
                bar = QRectF(centre_x - length, bar_top, length, self._BAR_HEIGHT)
                square_end = QRectF(centre_x - min(length, self._BAR_RADIUS), bar_top, min(length, self._BAR_RADIUS), self._BAR_HEIGHT)
                color = CONTRIBUTION_SENSITISING
            else:
                bar = QRectF(centre_x, bar_top, length, self._BAR_HEIGHT)
                square_end = QRectF(centre_x, bar_top, min(length, self._BAR_RADIUS), self._BAR_HEIGHT)
                color = CONTRIBUTION_RESISTANCE
            # Rounded at the outer (data) end, square where it meets zero.
            painter.setBrush(QColor(color))
            painter.drawRoundedRect(bar, self._BAR_RADIUS, self._BAR_RADIUS)
            painter.drawRect(square_end)

        # Zero line, drawn last so it sits on top of both bars.
        painter.setPen(QColor(TEXT_MUTED))
        painter.drawLine(QPointF(centre_x, bar_top - 3), QPointF(centre_x, bar_top + self._BAR_HEIGHT + 3))


class FinalResultsPage(QWidget):
    """Displays the model's predicted drug rankings and explains the picked drug."""

    def __init__(
        self,
        parent: QWidget | None = None,
        on_upload_clicked: Callable[[], None] | None = None,
        on_model_running_clicked: Callable[[], None] | None = None,
        on_analytics_clicked: Callable[[], None] | None = None,
    ) -> None:
        """Build the page.

        Args:
            parent: Optional Qt parent widget.
            on_upload_clicked: Invoked when the header's "Upload" tab is
                clicked; should navigate back to the dataset upload page.
            on_model_running_clicked: Invoked when the header's "Model
                Running" tab is clicked; should navigate to the model
                execution log page.
            on_analytics_clicked: Invoked when the sidebar's "Model
                Analytics" link is clicked; should navigate to the charts page.
        """
        super().__init__(parent)
        self.setObjectName("FinalResultsPage")
        self._on_upload_clicked = on_upload_clicked
        self._on_model_running_clicked = on_model_running_clicked
        self._on_analytics_clicked = on_analytics_clicked

        self._results_worker: ResultsWorker | None = None
        # Held on self because an unreferenced QThread is collected mid-flight.
        self._gene_worker: GeneAttributionWorker | None = None
        self._enrichment_worker: EnrichmentWorker | None = None
        self._run_id: str | None = None
        self._selected_drug_id: str | None = None
        self._visible_rows: list[DrugResult] = []
        self._all_rows: list[DrugResult] = []

        self._table: QTableWidget | None = None
        self._search_input: QLineEdit | None = None
        self._rank_filter_combo: QComboBox | None = None
        self._sort_combo: QComboBox | None = None
        self._sample_id_label: QLabel | None = None

        # Selected-drug pane.
        self._pane: QFrame | None = None
        self._pane_fade: QPropertyAnimation | None = None
        self._top_match_chip: QLabel | None = None
        self._drug_name_label: QLabel | None = None
        self._drug_ic50_value: QLabel | None = None
        self._drug_confidence_value: QLabel | None = None
        self._drug_rank_slot: QVBoxLayout | None = None
        self._drug_target_label: QLabel | None = None
        self._drug_pathway_label: QLabel | None = None
        self._drug_position_label: QLabel | None = None
        self._interpretation_host: QWidget | None = None
        self._pathways_button: QPushButton | None = None
        self._gene_table: QTableWidget | None = None
        self._enrichment_table: QTableWidget | None = None
        self._recovery_slot: QVBoxLayout | None = None

        self._build_ui()

    # -- live results -----------------------------------------------------

    def load_results(self, run_id: str) -> None:
        """Fetch and display the ranked drug predictions for `run_id`."""
        self._run_id = run_id
        self._selected_drug_id = None
        self._results_worker = ResultsWorker(run_id, parent=self)
        self._results_worker.succeeded.connect(self._on_results_succeeded)
        self._results_worker.failed.connect(self._on_results_failed)
        self._results_worker.start()

    def _on_results_succeeded(self, raw_results: list[dict]) -> None:
        rows = [
            DrugResult(
                drug_id=str(item["drug_id"]),
                name=item["drug_name"],
                ic50=item["predicted_ic50_um"],
                rank=item["ranking"],
                confidence=item["confidence_percent"],
                target=item.get("putative_target") or "",
                pathway=item.get("pathway_name") or "",
            )
            for item in raw_results
        ]
        rows.sort(key=lambda row: row.ic50)
        self._all_rows = rows
        self._apply_filters()

    def _on_results_failed(self, message: str) -> None:
        QMessageBox.critical(self, "Could Not Load Results", message)

    # -- export -------------------------------------------------------------

    def _on_download_data_clicked(self) -> None:
        """Save the full predicted drug results (`self._all_rows`) as CSV."""
        if not self._all_rows:
            QMessageBox.information(self, "No Results", "No results to export yet.")
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Download Data", "predicted_drug_results.csv", "CSV Files (*.csv)"
        )
        if not path:
            return

        try:
            with open(path, "w", newline="", encoding="utf-8") as file:
                writer = csv.writer(file)
                writer.writerow(["Drug Name", "Predicted IC50 (uM)", "Sensitivity Ranking", "Confidence (%)"])
                for row in self._all_rows:
                    writer.writerow([row.name, row.ic50, row.rank, row.confidence])
        except OSError as exc:
            QMessageBox.critical(self, "Could Not Save File", str(exc))
            return

        QMessageBox.information(self, "Download Complete", f"Saved to {path}")

    def _on_export_report_clicked(self) -> None:
        """Save a plain-text clinical summary report of `self._all_rows`."""
        if not self._all_rows:
            QMessageBox.information(self, "No Results", "No results to export yet.")
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Export Clinical Report", "clinical_report.txt", "Text Files (*.txt)"
        )
        if not path:
            return

        lines = [
            "Predicted Drug Sensitivity Report",
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Total drugs ranked: {len(self._all_rows)}",
            "",
        ]
        for index, row in enumerate(self._all_rows, start=1):
            lines.append(f"{index}. {row.name}")
            lines.append(f"   Predicted IC50: {format_ic50(row.ic50)} uM")
            lines.append(f"   Sensitivity Ranking: {row.rank}")
            lines.append(f"   Confidence: {row.confidence:.1f}%")
            lines.append("")

        try:
            with open(path, "w", encoding="utf-8") as file:
                file.write("\n".join(lines))
        except OSError as exc:
            QMessageBox.critical(self, "Could Not Save File", str(exc))
            return

        QMessageBox.information(self, "Export Complete", f"Saved to {path}")

    # -- layout -------------------------------------------------------------

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
        body_layout.addWidget(self._build_sample_strip())
        body_layout.addLayout(self._build_panels_grid())
        body_layout.addStretch(1)

        scroll.setWidget(body)
        shell_layout.addWidget(scroll, 1)
        root.addWidget(shell, 1)

    def _build_sidebar(self) -> QFrame:
        """Build the shared sidebar, with this page and the analytics page as nav links."""
        nav_widgets = [
            make_sidebar_nav_button("Drug Results", "insights", active=True),
            make_sidebar_nav_button("Model Analytics", "analytics", callback=self._on_analytics_clicked),
        ]
        footer_widgets = [make_sidebar_nav_button("Support", "help_outline")]
        return build_sidebar(nav_widgets, footer_widgets, cta_widget=make_primary_cta_button())

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
        subtitle = QLabel("Predicted drug response for your cell line")
        subtitle.setStyleSheet(PAGE_SUBTITLE_STYLE)
        title_col.addWidget(title)
        title_col.addWidget(subtitle)
        title_row.addLayout(title_col)
        title_row.addStretch(1)

        actions = QHBoxLayout()
        actions.setSpacing(12)
        download_btn = QPushButton(f"{icon_text('download')}  Download Data")
        download_btn.setStyleSheet(SECONDARY_BUTTON_STYLE)
        download_btn.clicked.connect(self._on_download_data_clicked)
        export_btn = QPushButton(f"{icon_text('summarize')}  Export Clinical Report")
        export_btn.setStyleSheet(PRIMARY_BUTTON_STYLE)
        export_btn.clicked.connect(self._on_export_report_clicked)
        actions.addWidget(download_btn)
        actions.addWidget(export_btn)
        title_row.addLayout(actions)
        return title_row

    def _build_sample_strip(self) -> QFrame:
        """Build the sample profile strip: the cell line and omics every result is about.

        A full-width strip above the results rather than a side card, since
        it's the context for everything on the page. Replaces the old
        "Molecular Profile" card, which fabricated a patient ID and gene
        mutations that don't exist anywhere in this cell-line-based pipeline.
        """
        strip = QFrame()
        strip.setObjectName("SampleStrip")
        strip.setStyleSheet(
            f"QFrame#SampleStrip {{ background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 10px; }}"
        )
        row = QHBoxLayout(strip)
        row.setContentsMargins(0, 0, 20, 0)
        row.setSpacing(20)

        accent = QFrame()
        accent.setFixedWidth(SAMPLE_PROFILE_ACCENT_WIDTH)
        accent.setStyleSheet(
            f"background: {SAMPLE_PROFILE_ACCENT}; border: none;"
            " border-top-left-radius: 9px; border-bottom-left-radius: 9px;"
        )
        row.addWidget(accent)

        sample_col = QVBoxLayout()
        sample_col.setContentsMargins(0, 14, 0, 14)
        sample_col.setSpacing(2)
        sample_title = QLabel("SAMPLE PROFILE")
        sample_title.setStyleSheet(
            label_style(f"font-size: 12px; font-weight: 700; letter-spacing: 0.05em; color: {TEXT};")
        )
        sample_id = QLabel("—")
        sample_id.setStyleSheet(
            label_style(
                f"font-family: Consolas, monospace; font-size: 22px; font-weight: 700; color: {TEXT};"
            )
        )
        self._sample_id_label = sample_id
        sample_col.addWidget(sample_title)
        sample_col.addWidget(sample_id)
        row.addLayout(sample_col)

        divider = QFrame()
        divider.setFixedWidth(1)
        divider.setStyleSheet(f"background: {SURFACE_CONTAINER}; border: none;")
        row.addWidget(divider)

        omics_col = QVBoxLayout()
        omics_col.setContentsMargins(0, 14, 0, 14)
        omics_col.setSpacing(6)
        omics_title = QLabel("OMICS MODALITIES USED")
        omics_title.setStyleSheet(LABEL_CAPS_STYLE)
        chips = QHBoxLayout()
        chips.setSpacing(8)
        chips.addWidget(self._chip("Transcriptomics", primary=True, help_key="omics_transcriptomics"))
        chips.addWidget(self._chip("Genomics (Mut/CNV)", help_key="omics_genomics"))
        chips.addWidget(self._chip("Proteomics", help_key="omics_proteomics"))
        chips.addStretch(1)
        omics_col.addWidget(omics_title)
        omics_col.addLayout(chips)
        row.addLayout(omics_col, 1)

        # !Hidden: Network Proximity is an empty placeholder box, out of
        # scope for now — see _build_network_proximity_placeholder below.
        # row.addWidget(self._build_network_proximity_placeholder())

        return strip

    def _build_panels_grid(self) -> QGridLayout:
        """Build the results table and, beside it, the selected-drug pane.

        Equal widths: the pane now carries the gene and pathway tables, so it
        needs as much room as the drug list.
        """
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(SECTION_SPACING)
        grid.addWidget(self._build_predicted_results_panel(), 0, 0)
        grid.addWidget(self._build_selected_drug_pane(), 0, 1)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        return grid

    # -- predicted results table ---------------------------------------------

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
        style_card_section(header, WINDOW_BACKGROUND, top=True)
        header_row = QHBoxLayout(header)
        header_row.setContentsMargins(16, 12, 16, 12)
        header_title = QLabel("Predicted Drug Results Panel")
        header_title.setStyleSheet(f"background: transparent; border: none; {CARD_TITLE_STYLE}")
        header_row.addWidget(header_title)
        header_row.addStretch(1)
        header_row.addWidget(make_info_icon("results_panel"))
        layout.addWidget(header)

        controls = QFrame()
        style_card_section(controls, SURFACE)
        controls_row = QHBoxLayout(controls)
        controls_row.setContentsMargins(16, 10, 16, 10)
        controls_row.setSpacing(8)

        search_input = QLineEdit()
        search_input.setPlaceholderText("Search drug name...")
        search_input.textChanged.connect(self._apply_filters)
        self._search_input = search_input

        # Short option labels keep both dropdowns narrow enough to sit beside
        # the search box now the table shares the page width with the pane.
        rank_filter_combo = QComboBox()
        rank_filter_combo.addItems(list(_RANK_FILTER_OPTIONS))
        rank_filter_combo.setFixedWidth(_RANK_FILTER_WIDTH)
        rank_filter_combo.currentTextChanged.connect(self._apply_filters)
        self._rank_filter_combo = rank_filter_combo

        sort_combo = QComboBox()
        sort_combo.addItems(list(self._SORT_OPTIONS))
        sort_combo.setFixedWidth(_SORT_WIDTH)
        sort_combo.currentTextChanged.connect(self._apply_filters)
        self._sort_combo = sort_combo

        controls_row.addWidget(search_input, 1)
        controls_row.addWidget(rank_filter_combo)
        controls_row.addWidget(sort_combo)
        layout.addWidget(controls)

        table = QTableWidget(0, 4)
        table.setHorizontalHeaderLabels(["Drug Name", "IC50 (µM)", "Sensitivity", "Confidence"])
        set_header_help(table, {1: "col_ic50", 2: "col_rank", 3: "col_confidence"})
        style_data_table(table, header_background=WINDOW_BACKGROUND)
        # `style_data_table` makes tables non-selectable, which is right for the
        # read-only ones. This one picks the drug the pane explains, so
        # selection is re-enabled to give that click feedback.
        table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        table.setCursor(Qt.CursorShape.PointingHandCursor)
        table.cellClicked.connect(self._on_row_clicked)
        table.setItemDelegateForColumn(0, _SelectedRowMarkerDelegate(table))
        # The default highlight is near-black, which swallows the dark
        # confidence label and bar fill painted by the cell widgets. A light
        # tint keeps every cell's own colours readable on the selected row.
        palette = table.palette()
        palette.setColor(QPalette.ColorRole.Highlight, QColor(SELECTED_DRUG_TINT))
        palette.setColor(QPalette.ColorRole.HighlightedText, QColor(TEXT))
        table.setPalette(palette)
        table.setStyleSheet(
            table.styleSheet()
            + f"QTableWidget::item:selected {{ background: {SELECTED_DRUG_TINT}; color: {TEXT}; }}"
        )
        table.setColumnWidth(0, 150)
        table.setColumnWidth(1, 96)
        table.setColumnWidth(2, 108)
        table.setMinimumHeight(480)
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

            ic50_item = QTableWidgetItem(format_ic50(row.ic50))
            ic50_item.setForeground(QColor(TEXT_MUTED))
            table.setItem(row_index, 1, ic50_item)

            tone = _RANK_TO_BADGE_TONE.get(row.rank, "neutral")
            table.setCellWidget(row_index, 2, build_status_badge(_RANK_TABLE_LABELS.get(row.rank, row.rank), tone))
            table.setCellWidget(row_index, 3, self._build_confidence_cell(row.confidence))

    # Sort dropdown label -> (sort key, reverse). Keeps `_apply_filters`
    # itself free of a long if/elif chain.
    _SORT_OPTIONS: dict[str, tuple[Callable[["DrugResult"], object], bool]] = {
        "Lowest IC50": (lambda row: row.ic50, False),
        "Highest IC50": (lambda row: row.ic50, True),
        "Most confident": (lambda row: row.confidence, True),
        "Least confident": (lambda row: row.confidence, False),
        "Name (A–Z)": (lambda row: row.name.lower(), False),
    }

    def _apply_filters(self) -> None:
        """Re-derive the table's rows from `self._all_rows` per the current
        search/filter/sort controls, leaving `self._all_rows` itself
        untouched.
        """
        if self._search_input is None or self._rank_filter_combo is None or self._sort_combo is None:
            return

        search_text = self._search_input.text().strip().lower()
        rank_filter = _RANK_FILTER_OPTIONS.get(self._rank_filter_combo.currentText())

        filtered = [
            row
            for row in self._all_rows
            if (not search_text or search_text in row.name.lower())
            and (rank_filter is None or row.rank == rank_filter)
        ]

        sort_key, reverse = self._SORT_OPTIONS.get(self._sort_combo.currentText(), (lambda row: row.ic50, False))
        filtered.sort(key=sort_key, reverse=reverse)

        empty_message = "No drugs match your search/filter." if self._all_rows else "Waiting for a completed model run..."
        self._visible_rows = filtered
        self._populate_table(filtered, empty_message=empty_message)

        # Explain the top row by default, so the pane is populated on arrival
        # rather than waiting for a click. Re-selects only when the previous
        # choice has been filtered away.
        if filtered and not any(row.drug_id == self._selected_drug_id for row in filtered):
            self._select_drug(filtered[0])
        self._sync_table_selection()

    def _sync_table_selection(self) -> None:
        """Mark the selected drug's row in the table, wherever sorting put it."""
        table = self._table
        if table is None:
            return
        for row_index, row in enumerate(self._visible_rows):
            if row.drug_id == self._selected_drug_id:
                table.selectRow(row_index)
                return
        table.clearSelection()

    def show_drug(self, drug_id: str) -> None:
        """Select the drug with `drug_id` and scroll its row into view.

        Used when a drug is picked elsewhere (the analytics page's scatter).
        Clears the search/filter first if they would hide it.
        """
        row = next((candidate for candidate in self._all_rows if candidate.drug_id == drug_id), None)
        if row is None or self._search_input is None or self._rank_filter_combo is None:
            return
        if not any(candidate.drug_id == drug_id for candidate in self._visible_rows):
            for control in (self._search_input, self._rank_filter_combo):
                control.blockSignals(True)
            self._search_input.clear()
            self._rank_filter_combo.setCurrentIndex(0)
            for control in (self._search_input, self._rank_filter_combo):
                control.blockSignals(False)
            self._apply_filters()
        if row.drug_id != self._selected_drug_id:
            self._select_drug(row)
        table = self._table
        if table is not None and table.currentRow() >= 0:
            table.scrollToItem(table.item(table.currentRow(), 0), QTableWidget.ScrollHint.PositionAtCenter)

    def _on_row_clicked(self, row_index: int, _column: int) -> None:
        if 0 <= row_index < len(self._visible_rows):
            row = self._visible_rows[row_index]
            if row.drug_id != self._selected_drug_id:
                self._select_drug(row)

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

    # -- selected-drug pane -----------------------------------------------------

    def _build_selected_drug_pane(self) -> QFrame:
        """Build the pane explaining the drug picked in the table.

        A dark header (inverted from every other card, so it stands out)
        carries the drug's prediction and annotation; below it, a segmented
        toggle switches between the genes and the pathways behind that
        prediction.
        """
        pane = QFrame()
        pane.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(pane)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._build_drug_header())

        interpretation = QWidget()
        interpretation_layout = QVBoxLayout(interpretation)
        interpretation_layout.setContentsMargins(0, 0, 0, 0)
        interpretation_layout.setSpacing(0)

        toggle_bar = QFrame()
        style_card_section(toggle_bar, SURFACE)
        toggle_row = QHBoxLayout(toggle_bar)
        toggle_row.setContentsMargins(16, 10, 16, 10)
        segment = QFrame()
        segment.setStyleSheet(f"background: {SELECTED_DRUG_TINT_STRONG}; border: none; border-radius: 8px;")
        segment_row = QHBoxLayout(segment)
        segment_row.setContentsMargins(3, 3, 3, 3)
        segment_row.setSpacing(2)
        genes_button = self._segment_button("Top genes", "gene_panel")
        pathways_button = self._segment_button("Pathways", "enrichment_panel")
        self._pathways_button = pathways_button
        segment_row.addWidget(genes_button)
        segment_row.addWidget(pathways_button)
        toggle_row.addWidget(segment)
        toggle_row.addStretch(1)
        interpretation_layout.addWidget(toggle_bar)

        stack = QStackedWidget()
        stack.addWidget(self._build_genes_view())
        stack.addWidget(self._build_pathways_view())
        interpretation_layout.addWidget(stack, 1)

        group = QButtonGroup(pane)
        group.setExclusive(True)
        group.addButton(genes_button, 0)
        group.addButton(pathways_button, 1)
        group.idClicked.connect(stack.setCurrentIndex)
        genes_button.setChecked(True)

        self._interpretation_host = interpretation
        layout.addWidget(interpretation, 1)
        self._pane = pane
        return pane

    @staticmethod
    def _segment_button(text: str, help_key: str) -> QPushButton:
        button = QPushButton(text)
        button.setCheckable(True)
        button.setCursor(Qt.CursorShape.PointingHandCursor)
        button.setStyleSheet(_SEGMENT_BUTTON_STYLE)
        button.setToolTip(HELP[help_key])
        return button

    def _build_drug_header(self) -> QFrame:
        """The pane's dark header: the picked drug's name, prediction and annotation."""
        header = QFrame()
        style_card_section(header, SELECTED_DRUG_SURFACE, top=True, divider=None)
        layout = QVBoxLayout(header)
        layout.setContentsMargins(22, 18, 22, 18)
        layout.setSpacing(10)

        top_row = QHBoxLayout()
        caption = QLabel("SELECTED DRUG")
        caption.setStyleSheet(
            label_style(f"font-size: 12px; font-weight: 700; letter-spacing: 0.08em; color: {_ON_DARK_MUTED};")
        )
        chip = QLabel("★  TOP PREDICTED MATCH")
        chip.setStyleSheet(
            f"background: {SURFACE}; color: {SELECTED_DRUG_SURFACE}; border: none; border-radius: 4px;"
            " padding: 3px 8px;"
            " font-size: 11px; font-weight: 700; letter-spacing: 0.05em;"
        )
        chip.setToolTip("Lowest predicted IC50 of all ranked drugs.")
        chip.setVisible(False)
        self._top_match_chip = chip
        top_row.addWidget(caption)
        top_row.addStretch(1)
        top_row.addWidget(chip)
        layout.addLayout(top_row)

        name = QLabel("—")
        name.setWordWrap(True)
        name.setStyleSheet(label_style(f"font-size: 26px; font-weight: 600; color: {SELECTED_DRUG_TEXT};"))
        self._drug_name_label = name
        layout.addWidget(name)

        metrics = QGridLayout()
        metrics.setHorizontalSpacing(28)
        metrics.setVerticalSpacing(2)
        ic50_value = QLabel("—")
        confidence_value = QLabel("—")
        for column, (caption_text, value, help_key) in enumerate(
            (("IC50 (µM)", ic50_value, "col_ic50"), ("CONFIDENCE", confidence_value, "col_confidence"))
        ):
            metric_caption = QLabel(caption_text)
            metric_caption.setStyleSheet(
                label_style(f"font-size: 11px; font-weight: 700; letter-spacing: 0.05em; color: {_ON_DARK_MUTED};")
            )
            metric_caption.setToolTip(HELP[help_key])
            value.setStyleSheet(
                label_style(
                    f"font-family: Consolas, monospace; font-size: 18px; font-weight: 700; color: {SELECTED_DRUG_TEXT};"
                )
            )
            metrics.addWidget(metric_caption, 0, column)
            metrics.addWidget(value, 1, column)
        ranking_caption = QLabel("RANKING")
        ranking_caption.setStyleSheet(
            label_style(f"font-size: 11px; font-weight: 700; letter-spacing: 0.05em; color: {_ON_DARK_MUTED};")
        )
        ranking_caption.setToolTip(HELP["col_rank"])
        rank_slot = QVBoxLayout()
        rank_slot.setContentsMargins(0, 2, 0, 0)
        metrics.addWidget(ranking_caption, 0, 2)
        metrics.addLayout(rank_slot, 1, 2)
        metrics.setColumnStretch(3, 1)
        self._drug_ic50_value = ic50_value
        self._drug_confidence_value = confidence_value
        self._drug_rank_slot = rank_slot
        layout.addLayout(metrics)

        annotation_style = label_style(f"font-size: 13px; color: {_ON_DARK_MUTED};")
        target = QLabel()
        pathway = QLabel()
        position = QLabel("Pick a drug in the table")
        for label in (target, pathway, position):
            label.setWordWrap(True)
            label.setStyleSheet(annotation_style)
        target.setVisible(False)
        pathway.setVisible(False)
        self._drug_target_label = target
        self._drug_pathway_label = pathway
        self._drug_position_label = position
        layout.addWidget(target)
        layout.addWidget(pathway)
        layout.addWidget(position)
        return header

    def _populate_drug_details(self, row: DrugResult) -> None:
        """Fill the pane's dark header from the picked drug's real prediction."""
        is_top_match = bool(self._all_rows) and row.drug_id == self._all_rows[0].drug_id
        if self._top_match_chip is not None:
            self._top_match_chip.setVisible(is_top_match)
        if self._drug_name_label is not None:
            self._drug_name_label.setText(row.name)
        if self._drug_ic50_value is not None:
            self._drug_ic50_value.setText(format_ic50(row.ic50))
        if self._drug_confidence_value is not None:
            self._drug_confidence_value.setText(f"{row.confidence:.1f}%")
        if self._drug_rank_slot is not None:
            while self._drug_rank_slot.count():
                item = self._drug_rank_slot.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
            tone = _RANK_TO_BADGE_TONE.get(row.rank, "neutral")
            self._drug_rank_slot.addWidget(build_status_badge(row.rank, tone))

        for label, prefix, value in (
            (self._drug_target_label, "Target", row.target),
            (self._drug_pathway_label, "Pathway", row.pathway),
        ):
            if label is not None:
                label.setText(f"{prefix}  ·  {value}")
                label.setVisible(bool(value))
        if self._drug_position_label is not None:
            position = next(
                (index for index, candidate in enumerate(self._all_rows, start=1) if candidate.drug_id == row.drug_id),
                None,
            )
            self._drug_position_label.setText(
                f"#{position} of {len(self._all_rows)} by predicted IC50" if position is not None else ""
            )

    def _fade_in_pane(self) -> None:
        """Briefly fade the pane in, so the eye follows a table click across to it.

        The opacity effect is removed once the fade ends: a lingering graphics
        effect would re-render the whole pane (tables included) offscreen on
        every repaint.
        """
        pane = self._pane
        if pane is None:
            return
        if self._pane_fade is not None:
            self._pane_fade.stop()
        effect = QGraphicsOpacityEffect(pane)
        pane.setGraphicsEffect(effect)
        fade = QPropertyAnimation(effect, b"opacity", self)
        fade.setDuration(180)
        fade.setStartValue(0.35)
        fade.setEndValue(1.0)
        fade.finished.connect(lambda: pane.setGraphicsEffect(None))
        fade.start()
        self._pane_fade = fade

    # -- biological interpretation -----------------------------------------

    def set_supports_interpretation(self, supported: bool) -> None:
        """Show or hide the genes/pathways views for the current run's backend.

        Not every model can attribute a prediction to genes -- one trained on
        PCA-projected omics has no per-gene identity left -- so the backend
        declares the capability and the views disappear entirely rather than
        sitting permanently empty. The drug header stays either way.
        """
        if self._interpretation_host is not None:
            self._interpretation_host.setVisible(supported)

    def _select_drug(self, row: DrugResult) -> None:
        """Show `row` in the pane and fetch its gene attribution and enrichment."""
        self._selected_drug_id = row.drug_id
        self._populate_drug_details(row)
        self._fade_in_pane()
        self._sync_table_selection()
        if self._run_id is None:
            return

        self._populate_gene_table([], empty_message="Attributing...")
        self._populate_recovery({"checked": None})
        self._populate_enrichment_table([], empty_message="Querying Enrichr...")
        if self._pathways_button is not None:
            self._pathways_button.setText("Pathways")

        self._gene_worker = GeneAttributionWorker(self._run_id, row.drug_id, parent=self)
        self._gene_worker.succeeded.connect(self._on_gene_attribution_succeeded)
        self._gene_worker.failed.connect(self._on_gene_attribution_failed)
        self._gene_worker.unsupported.connect(lambda: self.set_supports_interpretation(False))
        self._gene_worker.start()

        self._enrichment_worker = EnrichmentWorker(self._run_id, row.drug_id, parent=self)
        self._enrichment_worker.succeeded.connect(self._on_enrichment_succeeded)
        self._enrichment_worker.failed.connect(self._on_enrichment_failed)
        self._enrichment_worker.start()

    def _on_gene_attribution_succeeded(self, payload: dict) -> None:
        self._populate_gene_table(payload.get("genes", []))
        self._populate_recovery(payload.get("target_recovery", {}))

    def _on_gene_attribution_failed(self, message: str) -> None:
        # Supplementary view: degrade in place rather than interrupting with a
        # modal.
        self._populate_gene_table([], empty_message=message)

    def _on_enrichment_succeeded(self, payload: dict) -> None:
        terms = payload.get("terms", [])
        message = payload.get("message", "")
        self._populate_enrichment_table(
            terms, empty_message=message or "No significantly enriched pathways."
        )
        if self._pathways_button is not None:
            self._pathways_button.setText(f"Pathways · {len(terms)}" if terms else "Pathways")

    def _on_enrichment_failed(self, message: str) -> None:
        self._populate_enrichment_table([], empty_message=message)
        if self._pathways_button is not None:
            self._pathways_button.setText("Pathways")

    def _build_genes_view(self) -> QWidget:
        view = QWidget()
        layout = QVBoxLayout(view)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        recovery_strip = QFrame()
        style_card_section(recovery_strip, SURFACE)
        recovery_strip.setToolTip(HELP["target_recovery"])
        recovery_layout = QVBoxLayout(recovery_strip)
        recovery_layout.setContentsMargins(16, 10, 16, 10)
        self._recovery_slot = recovery_layout
        layout.addWidget(recovery_strip)

        # No separate Direction column: the contribution bar's side (left =
        # sensitising, right = resistance) and its signed value already carry it.
        table = QTableWidget(0, 3)
        table.setHorizontalHeaderLabels(["Gene", "Contribution", "Evidence"])
        set_header_help(table, {1: "col_contribution", 2: "col_evidence"})
        style_data_table(table, header_background=SELECTED_DRUG_TINT)
        # Span the pane's full width at any window size: Gene and Evidence are
        # fixed to fit their content (a symbol; the widest badge, "Driver
        # mutation"), and the contribution bar gets all the rest.
        header = table.horizontalHeader()
        header.setStretchLastSection(False)
        for column, (mode, width) in enumerate(
            (
                (QHeaderView.ResizeMode.Fixed, 84),
                (QHeaderView.ResizeMode.Stretch, None),
                (QHeaderView.ResizeMode.Fixed, 140),
            )
        ):
            header.setSectionResizeMode(column, mode)
            if width is not None:
                table.setColumnWidth(column, width)
        table.setMinimumHeight(300)
        self._gene_table = table
        layout.addWidget(table, 1)

        self._populate_gene_table([])
        self._populate_recovery({"checked": None})
        return view

    def _build_pathways_view(self) -> QWidget:
        table = QTableWidget(0, 4)
        table.setHorizontalHeaderLabels(["Pathway / Term", "Library", "Adj. p", "Overlap"])
        set_header_help(table, {2: "col_adj_p", 3: "col_overlap"})
        style_data_table(table, header_background=SELECTED_DRUG_TINT)
        # Long term names take whatever width is left and elide; the full name
        # is in each row's tooltip.
        header = table.horizontalHeader()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        table.setColumnWidth(1, 88)
        table.setColumnWidth(2, 84)
        table.setColumnWidth(3, 80)
        table.setMinimumHeight(300)
        self._enrichment_table = table
        self._populate_enrichment_table([])
        return table

    def _populate_gene_table(
        self, genes: list[dict], empty_message: str = "Pick a drug to see the genes behind it."
    ) -> None:
        table = self._gene_table
        if table is None:
            return
        if not genes:
            table.setRowCount(1)
            item = QTableWidgetItem(empty_message)
            item.setForeground(QColor(TEXT_MUTED))
            table.setItem(0, 0, item)
            table.setSpan(0, 0, 1, 3)
            return

        table.clearSpans()
        table.setRowCount(len(genes))
        max_abs = max(abs(float(gene["score"])) for gene in genes)
        for row_index, gene in enumerate(genes):
            score = float(gene["score"])

            symbol_item = QTableWidgetItem(str(gene["gene_symbol"]))
            symbol_item.setForeground(QColor(TEXT))
            table.setItem(row_index, 0, symbol_item)

            table.setCellWidget(row_index, 1, self._build_contribution_cell(score, max_abs))

            if gene.get("is_drug_target"):
                badge = build_status_badge("Drug target", "positive")
            elif gene.get("is_driver_mutation"):
                badge = build_status_badge("Driver mutation", "neutral")
            else:
                badge = build_status_badge("Network", "muted")
            table.setCellWidget(row_index, 2, badge)

    @staticmethod
    def _build_contribution_cell(score: float, max_abs: float) -> QWidget:
        """Signed contribution as a diverging bar (see `_ContributionBar`)."""
        wrapper = transparent_cell_widget()
        layout = QVBoxLayout(wrapper)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(_ContributionBar(score, max_abs))
        return wrapper

    def _populate_recovery(self, recovery: dict) -> None:
        """Show whether the drug's GDSC-annotated target came back in the top genes.

        `{"checked": None}` means "not known yet" (no drug picked, or still
        attributing) and shows nothing.
        """
        slot = self._recovery_slot
        if slot is None:
            return
        while slot.count():
            item = slot.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

        if recovery.get("checked") is None:
            slot.addWidget(build_status_badge("Checking known target...", "muted"))
            return
        if not recovery.get("checked"):
            # `putative_target` is frequently a mechanism ("Microtubule
            # destabiliser") rather than a gene, so there is nothing to check
            # against -- say so instead of implying a failed recovery.
            slot.addWidget(build_status_badge("No gene-level target annotated", "muted"))
            return

        recovered = recovery.get("recovered", [])
        if recovered:
            best = min(entry["rank"] for entry in recovered)
            names = ", ".join(entry["gene_symbol"] for entry in recovered)
            slot.addWidget(
                build_status_badge(f"Known target recovered: {names} (rank {best})", "positive")
            )
        else:
            targets = ", ".join(recovery.get("target_genes", []))
            slot.addWidget(build_status_badge(f"Known target not in top genes ({targets})", "muted"))

    def _populate_enrichment_table(
        self, terms: list[dict], empty_message: str = "Pick a drug to see enriched pathways."
    ) -> None:
        table = self._enrichment_table
        if table is None:
            return
        if not terms:
            table.setRowCount(1)
            item = QTableWidgetItem(empty_message)
            item.setForeground(QColor(TEXT_MUTED))
            table.setItem(0, 0, item)
            table.setSpan(0, 0, 1, 4)
            return

        table.clearSpans()
        table.setRowCount(len(terms))
        for row_index, term in enumerate(terms):
            term_item = QTableWidgetItem(str(term["term"]))
            term_item.setForeground(QColor(TEXT))
            # The name is elided in the narrow pane, and the overlapping genes
            # are the evidence behind the term: both are worth having on hover.
            genes = ", ".join(term.get("genes", []))
            term_item.setToolTip(f"{term['term']}\n{genes}" if genes else str(term["term"]))
            table.setItem(row_index, 0, term_item)

            library = str(term["library"])
            library_item = QTableWidgetItem(_LIBRARY_SHORT_NAMES.get(library, library.split("_")[0]))
            library_item.setForeground(QColor(TEXT_MUTED))
            library_item.setToolTip(library.replace("_", " "))
            table.setItem(row_index, 1, library_item)

            p_item = QTableWidgetItem(f"{float(term['adjusted_p_value']):.1e}")
            p_item.setForeground(QColor(TEXT_MUTED))
            table.setItem(row_index, 2, p_item)

            overlap_item = QTableWidgetItem(str(term["overlap"]))
            overlap_item.setForeground(QColor(TEXT_MUTED))
            table.setItem(row_index, 3, overlap_item)

    # -- sample profile ---------------------------------------------------------

    def set_sample_id(self, sanger_model_id: str) -> None:
        """Update the sample strip's cell-line identifier."""
        if self._sample_id_label is not None:
            self._sample_id_label.setText(sanger_model_id)

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
    def _chip(text: str, primary: bool = False, help_key: str | None = None) -> QLabel:
        """Build a small rounded label chip (used for omics modality tags).

        Args:
            text: Chip label.
            primary: Whether to render as a filled (primary) chip versus an
                outlined (secondary) chip.
            help_key: Optional key into `widgets.help.HELP` shown on hover.

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
        if help_key is not None:
            chip.setToolTip(HELP[help_key])
        return chip
