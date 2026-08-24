"""Shared helpers for building and extending the heterogeneous graph.

`03_graph_construction.py` (initial build) and `04_link_cell_lines.py`
(cell-line linkage) both need to translate free-text gene symbols into
STRING's Ensembl protein IDs, so that mapping lives here rather than being
duplicated across the two scripts.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, Set

STRING_ALIASES_PATH = Path("data/raw/string/9606.protein.aliases.v12.0.txt.gz")


def build_gene_symbol_to_ensp(
    referenced_symbols: Set[str], aliases_path: Path = STRING_ALIASES_PATH
) -> Dict[str, str]:
    """Build a gene_symbol -> ENSP map from the STRING alias file, filtered to referenced symbols.

    The full alias file has ~40M rows; loading it wholesale into memory is
    unnecessary when only the referenced symbol vocabulary needs resolving, so
    this does a single streamed pass and only retains rows whose alias is in
    `referenced_symbols`. The first ENSP seen for a symbol is kept (STRING
    lists the same gene symbol from multiple alias sources for the same
    canonical protein in practice; no attempt is made to adjudicate
    conflicting sources beyond first-match).
    """

    symbol_to_ensp: Dict[str, str] = {}
    with gzip.open(aliases_path, "rt", encoding="utf-8") as handle:
        next(handle)  # header: #string_protein_id, alias, source
        for line in handle:
            protein_id, alias, _source = line.rstrip("\n").split("\t")
            if alias in referenced_symbols and alias not in symbol_to_ensp:
                symbol_to_ensp[alias] = protein_id
    return symbol_to_ensp
