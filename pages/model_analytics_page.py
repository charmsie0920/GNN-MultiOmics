"""Model Analytics page.

The charts that used to sit under the results table, given a page of their
own so the results page stays focused on the drugs. Reached from the
sidebar's "Model Analytics" link on the results page.

Main UI components:
    - Shared sidebar with "Drug Results" / "Model Analytics" links, the
      latter active; and the shared header with "Results" active, since this
      page is a view of the same results.
    - A KPI row of stat tiles, computed from this run's real results and
      training history.
    - A two-column grid of chart cards (`widgets.charts`). Add a chart with
      one `_add_chart_card` call; cards fill the grid row by row.

Like the results page, it loads its own data from `run_id` (see
`load_results`), so the two pages stay decoupled.
"""

from __future__ import annotations

import statistics
from collections.abc import Callable

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QFrame, QGridLayout, QHBoxLayout, QLabel, QScrollArea, QVBoxLayout, QWidget

from client.workers import ResultsWorker, TrainingHistoryWorker
from styles.theme import (
    CARD_CONTAINER_STYLE,
    PAGE_MARGIN,
    PAGE_SUBTITLE_STYLE,
    PAGE_TITLE_STYLE,
    SECTION_SPACING,
    TEXT,
    TEXT_MUTED,
    TEXT_SOFT,
    label_style,
)
from widgets.charts import (
    HIGH_CONFIDENCE_COLOR,
    LOW_CONFIDENCE_COLOR,
    TRAIN_LOSS_COLOR,
    VAL_RMSE_COLOR,
    VALIDATION_POINT_COLOR,
    PredictionScatterWidget,
    TrainingCurveWidget,
)
from widgets.help import HELP, make_info_icon
from widgets.navigation import (
    build_header_bar,
    build_sidebar,
    make_primary_cta_button,
    make_sidebar_nav_button,
    make_top_tab,
)


def _legend_html(items: list[tuple[str, str, str]]) -> str:
    """Legend row as rich text: a coloured key glyph, then its label in text ink."""
    return "&nbsp;&nbsp;&nbsp;&nbsp;".join(
        f'<span style="color:{color};">{glyph}</span>&nbsp;<span style="color:{TEXT_MUTED};">{label}</span>'
        for glyph, color, label in items
    )


class _StatTile(QFrame):
    """One KPI: a sentence-case label, a large value, and a muted detail line."""

    def __init__(self, label: str, help_text: str | None = None) -> None:
        super().__init__()
        self.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 16, 20, 16)
        layout.setSpacing(4)

        title = QLabel(label)
        title.setStyleSheet(label_style(f"font-size: 13px; color: {TEXT_MUTED};"))
        self._value = QLabel("—")
        self._value.setStyleSheet(label_style(f"font-size: 30px; font-weight: 600; color: {TEXT};"))
        self._detail = QLabel(" ")
        self._detail.setStyleSheet(label_style(f"font-size: 12px; color: {TEXT_SOFT};"))
        layout.addWidget(title)
        layout.addWidget(self._value)
        layout.addWidget(self._detail)
        if help_text:
            self.setToolTip(help_text)

    def set_value(self, value: str, detail: str = " ") -> None:
        self._value.setText(value)
        self._detail.setText(detail)


class _ChartCard(QFrame):
    """A chart with its title, one-line subtitle, "ⓘ" help and legend row."""

    def __init__(self, title: str, chart: QWidget, help_key: str) -> None:
        super().__init__()
        self.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 18, 20, 16)
        layout.setSpacing(4)

        title_row = QHBoxLayout()
        title_label = QLabel(title)
        title_label.setStyleSheet(label_style(f"font-size: 16px; font-weight: 600; color: {TEXT};"))
        self.info_icon = make_info_icon(help_key)
        title_row.addWidget(title_label)
        title_row.addStretch(1)
        title_row.addWidget(self.info_icon)
        layout.addLayout(title_row)

        self._subtitle = QLabel(" ")
        self._subtitle.setStyleSheet(label_style(f"font-size: 12px; color: {TEXT_SOFT};"))
        layout.addWidget(self._subtitle)

        self._legend = QLabel(" ")
        self._legend.setTextFormat(Qt.TextFormat.RichText)
        self._legend.setStyleSheet(label_style("font-size: 12px;"))
        layout.addWidget(self._legend)
        layout.addSpacing(6)
        layout.addWidget(chart, 1)

    def set_subtitle(self, text: str) -> None:
        self._subtitle.setText(text)

    def set_legend(self, items: list[tuple[str, str, str]]) -> None:
        """`items`: (glyph, colour, label), e.g. ("●", "#8a3419", "Train loss")."""
        self._legend.setText(_legend_html(items) if items else " ")


class ModelAnalyticsPage(QWidget):
    """KPIs and hoverable charts describing how the model performed on this run."""

    def __init__(
        self,
        parent: QWidget | None = None,
        on_upload_clicked: Callable[[], None] | None = None,
        on_model_running_clicked: Callable[[], None] | None = None,
        on_results_clicked: Callable[[], None] | None = None,
        on_drug_clicked: Callable[[str], None] | None = None,
    ) -> None:
        """Build the page.

        Args:
            parent: Optional Qt parent widget.
            on_upload_clicked: Header "Upload" tab; back to the upload page.
            on_model_running_clicked: Header "Model Running" tab.
            on_results_clicked: Sidebar "Drug Results" link and header
                "Results" tab; back to the results page.
            on_drug_clicked: Called with a drug id when its dot in the IC50
                vs. confidence chart is clicked; should open that drug on
                the results page.
        """
        super().__init__(parent)
        self.setObjectName("ModelAnalyticsPage")
        self._on_upload_clicked = on_upload_clicked
        self._on_model_running_clicked = on_model_running_clicked
        self._on_results_clicked = on_results_clicked
        self._on_drug_clicked = on_drug_clicked

        # Held on self because an unreferenced QThread is collected mid-flight.
        self._results_worker: ResultsWorker | None = None
        self._training_history_worker: TrainingHistoryWorker | None = None
        self._chart_grid: QGridLayout | None = None
        self._chart_count = 0

        self._build_ui()

    # -- data -------------------------------------------------------------

    def load_results(self, run_id: str) -> None:
        """Fetch this run's ranked predictions and training history."""
        for tile in (self._drugs_tile, self._high_tile, self._confidence_tile, self._rmse_tile):
            tile.set_value("—")
        self._scatter.set_points([])
        self._performance.set_history([])

        self._results_worker = ResultsWorker(run_id, parent=self)
        self._results_worker.succeeded.connect(self._on_results_succeeded)
        self._results_worker.start()

        self._training_history_worker = TrainingHistoryWorker(run_id, parent=self)
        self._training_history_worker.succeeded.connect(self._on_training_history_succeeded)
        self._training_history_worker.start()

    def set_sample_id(self, sanger_model_id: str) -> None:
        self._subtitle.setText(f"How the model performed for cell line {sanger_model_id}")

    def _on_results_succeeded(self, raw_results: list[dict]) -> None:
        # No failure handler: the results page already reports a failed load,
        # and this page just keeps its "waiting" state.
        count = len(raw_results)
        high = sum(1 for item in raw_results if item["ranking"] == "HIGH SENSITIVITY")
        self._drugs_tile.set_value(f"{count:,}", "GDSC2 compounds scored")
        self._high_tile.set_value(f"{high:,}", f"of {count:,} drugs · IC50 < 1 µM")
        if raw_results:
            median = statistics.median(item["confidence_percent"] for item in raw_results)
            self._confidence_tile.set_value(f"{median:.0f}%", "across all ranked drugs")
        self._scatter.set_points(
            [
                (item["predicted_ic50_um"], item["confidence_percent"], item["drug_name"], str(item["drug_id"]))
                for item in raw_results
            ]
        )

    def _on_training_history_succeeded(self, history: list[dict]) -> None:
        self._performance.set_history(history)
        mode = self._performance.mode()
        card = self._performance_card
        if mode == "curve":
            card.info_icon.setToolTip(HELP["perf_curve"])
            card.set_subtitle("Per epoch · each metric on its own scale")
            card.set_legend([("●", TRAIN_LOSS_COLOR, "Train loss"), ("●", VAL_RMSE_COLOR, "Validation RMSE")])
        elif mode == "scatter":
            card.info_icon.setToolTip(HELP["perf_checkpoint"])
            card.set_subtitle("Pretrained checkpoint · held-out validation samples")
            card.set_legend([("●", VALIDATION_POINT_COLOR, "Validation sample"), ("┄", TEXT_SOFT, "Perfect prediction")])

        # Always the backend's own figure (the one it logs and scales
        # confidence by), never recomputed here from the plotted points.
        eval_rmse = history[0].get("eval_val_rmse") if history else None
        if eval_rmse is not None:
            self._rmse_tile.set_value(f"{eval_rmse:.4f}", "full validation set")
        else:
            self._rmse_tile.set_value("—", "not reported by this backend")

    # -- layout -----------------------------------------------------------

    def _build_ui(self) -> None:
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

        title = QLabel("Model Analytics")
        title.setStyleSheet(PAGE_TITLE_STYLE)
        self._subtitle = QLabel("How the model performed on this run")
        self._subtitle.setStyleSheet(PAGE_SUBTITLE_STYLE)
        title_col = QVBoxLayout()
        title_col.setSpacing(4)
        title_col.addWidget(title)
        title_col.addWidget(self._subtitle)
        body_layout.addLayout(title_col)

        body_layout.addLayout(self._build_kpi_row())

        self._chart_grid = QGridLayout()
        self._chart_grid.setHorizontalSpacing(SECTION_SPACING)
        self._chart_grid.setVerticalSpacing(SECTION_SPACING)
        self._chart_grid.setColumnStretch(0, 1)
        self._chart_grid.setColumnStretch(1, 1)
        body_layout.addLayout(self._chart_grid)

        self._performance = TrainingCurveWidget()
        self._performance_card = _ChartCard("Model Performance", self._performance, "perf_checkpoint")
        self._add_chart_card(self._performance_card)

        self._scatter = PredictionScatterWidget()
        if self._on_drug_clicked is not None:
            self._scatter.drugClicked.connect(self._on_drug_clicked)
        scatter_card = _ChartCard("Predicted IC50 vs. Confidence", self._scatter, "ic50_vs_confidence")
        scatter_card.set_subtitle("One dot per drug · top-left = strongest candidates")
        scatter_card.set_legend(
            [("●", LOW_CONFIDENCE_COLOR.name(), "Low confidence"), ("●", HIGH_CONFIDENCE_COLOR.name(), "High confidence")]
        )
        self._add_chart_card(scatter_card)

        body_layout.addStretch(1)
        scroll.setWidget(body)
        shell_layout.addWidget(scroll, 1)
        root.addWidget(shell, 1)

    def _build_kpi_row(self) -> QHBoxLayout:
        row = QHBoxLayout()
        row.setSpacing(SECTION_SPACING)
        self._drugs_tile = _StatTile("Drugs ranked")
        self._high_tile = _StatTile("High sensitivity")
        self._confidence_tile = _StatTile("Median confidence")
        self._rmse_tile = _StatTile("Validation RMSE", HELP["kpi_rmse"])
        for tile in (self._drugs_tile, self._high_tile, self._confidence_tile, self._rmse_tile):
            row.addWidget(tile, 1)
        return row

    def _add_chart_card(self, card: QFrame) -> None:
        """Place `card` in the next free slot of the two-column chart grid."""
        row, column = divmod(self._chart_count, 2)
        self._chart_grid.addWidget(card, row, column)
        self._chart_count += 1

    def _build_sidebar(self) -> QFrame:
        nav_widgets = [
            make_sidebar_nav_button("Drug Results", "insights", callback=self._on_results_clicked),
            make_sidebar_nav_button("Model Analytics", "analytics", active=True),
        ]
        footer_widgets = [make_sidebar_nav_button("Support", "help_outline")]
        return build_sidebar(nav_widgets, footer_widgets, cta_widget=make_primary_cta_button())

    def _build_header(self) -> QFrame:
        tabs = [
            make_top_tab("Upload", callback=self._on_upload_clicked),
            make_top_tab("Model Running", callback=self._on_model_running_clicked),
            make_top_tab("Results", active=True, callback=self._on_results_clicked),
        ]
        return build_header_bar(tabs)
