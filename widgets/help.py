"""Hover help: every tooltip string in the app, plus the helpers that attach them.

Keeping the copy here rather than inline in `pages/` means all of it can be
reviewed and reworded in one place, the same way `styles.theme` centralizes
colors. Pages refer to entries by key only.

Copy rules: a few words, or one short sentence at most -- readable in under a
second. Only add a tooltip where the value isn't self-explanatory; a label
that already says what it does doesn't need one.
"""

from __future__ import annotations

from collections.abc import Mapping

from PySide6.QtWidgets import QLabel, QTableWidget

from styles.theme import TEXT_MUTED, label_style
from widgets.icons import icon_text

HELP: dict[str, str] = {
    # -- dataset initialization page --
    "upload": "CSV with sanger_model_id, drug_id and ln_ic50 columns.",
    # -- model execution log page --
    "overall_progress": "Estimated from elapsed time.",
    "file_graph": "Prebuilt graph of cell lines, drugs and proteins.",
    "file_dataset": "Your uploaded dataset.",
    "file_checkpoint": "Trained model weights.",
    "file_inference": "Predicting each drug's response.",
    # -- final results page: predicted drug results --
    "results_panel": "Click a row to see the genes behind it.",
    "col_ic50": "Lower = more sensitive.",
    "col_rank": "High < 1 µM · Medium 1–10 µM · Low ≥ 10 µM",
    "col_confidence": "How consistent the model's repeated predictions are.",
    # -- final results page: charts --
    "perf_curve": "Falling validation RMSE = model improving.",
    "perf_checkpoint": "Dots near the diagonal = accurate predictions.",
    "ic50_vs_confidence": "Top-left = strongest candidates. Hover a dot to see it.",
    "kpi_rmse": "Typical prediction error on held-out data, in ln(IC50).",
    # -- final results page: sample profile --
    "omics_transcriptomics": "Gene expression (RNA) levels measured in the cell line.",
    "omics_genomics": "DNA mutations and gene copy-number gains or losses.",
    "omics_proteomics": "Protein abundance measured in the cell line.",
    # -- final results page: biological interpretation --
    "gene_panel": "Genes that most influenced this prediction.",
    "col_contribution": "Left = toward sensitivity, right = toward resistance (√ scale).",
    "col_evidence": "Drug target, driver mutation, or network link.",
    "target_recovery": "Is the drug's known target among the top genes?",
    "enrichment_panel": "Pathways enriched in the top genes (needs internet).",
    "col_adj_p": "Lower = stronger; < 0.05 is significant.",
    "col_overlap": "Top genes in pathway / pathway size.",
}


def make_info_icon(key: str) -> QLabel:
    """Build the small "ⓘ" glyph placed beside a card title, showing `HELP[key]` on hover.

    The visible icon is what tells users there's help to find; plain hover
    tooltips on unmarked labels are rarely discovered.
    """
    icon = QLabel(icon_text("info"))
    icon.setStyleSheet(label_style(f"font-size: 16px; color: {TEXT_MUTED};"))
    icon.setToolTip(HELP[key])
    return icon


def set_header_help(table: QTableWidget, column_keys: Mapping[int, str]) -> None:
    """Show `HELP[key]` when hovering the header of each given column.

    Call after `setHorizontalHeaderLabels`, which is what creates the header
    items these tooltips are attached to.
    """
    for column, key in column_keys.items():
        item = table.horizontalHeaderItem(column)
        if item is not None:
            item.setToolTip(HELP[key])
