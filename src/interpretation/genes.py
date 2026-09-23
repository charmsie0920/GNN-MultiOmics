"""Rank and annotate per-protein attribution scores into gene-level rows.

Deliberately model-agnostic: this takes a score vector and knows nothing about
how it was produced, so a gradient-based scorer, a GNNExplainer mask, or an
attention readout all feed the same downstream enrichment and UI.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Sequence

SYMBOL_MAP_PATH = Path("data/processed/protein_symbol_map.csv")


def load_symbol_map(path: Path = SYMBOL_MAP_PATH) -> dict[str, str]:
    """Load the ENSP -> gene symbol map written by src/data/build_protein_symbol_map.py."""
    if not path.exists():
        raise FileNotFoundError(
            f"Gene symbol map not found at {path}. Build it with "
            "src/data/build_protein_symbol_map.py (after src/data/fetch_string_aliases.py)."
        )
    with path.open(encoding="utf-8") as handle:
        return {row["protein_id"]: row["gene_symbol"] for row in csv.DictReader(handle)}


def annotate_genes(
    scores: Sequence[float],
    node_ids: Sequence[str],
    symbol_map: dict[str, str],
    top_k: int = 50,
    *,
    driver_proteins: Iterable[str] = (),
    target_proteins: Iterable[str] = (),
) -> list[dict]:
    """Rank proteins by absolute attribution and annotate the top `top_k`.

    The score's *sign* carries the interpretation: the model predicts ln(IC50),
    so a negative contribution pushes the prediction down, i.e. toward
    sensitivity. Ranking is on magnitude because a strong resistance signal is
    just as interesting as a strong sensitivity one.

    `driver_proteins` and `target_proteins` are passed in as plain id sets
    rather than derived here, so this stays independent of any particular
    graph's edge semantics.
    """
    drivers = set(driver_proteins)
    targets = set(target_proteins)

    ranked = sorted(range(len(scores)), key=lambda i: abs(scores[i]), reverse=True)

    rows: list[dict] = []
    for index in ranked[:top_k]:
        protein_id = node_ids[index]
        score = float(scores[index])
        rows.append(
            {
                "gene_symbol": symbol_map.get(protein_id) or protein_id,
                "protein_id": protein_id,
                "score": score,
                "direction": "sensitising" if score < 0 else "resistance",
                "is_driver_mutation": protein_id in drivers,
                "is_drug_target": protein_id in targets,
            }
        )
    return rows


def parse_target_symbols(putative_target: str, known_symbols: set[str]) -> tuple[list[str], list[str]]:
    """Split GDSC's `putative_target` into recognised gene symbols and leftovers.

    The column mixes real symbols ("PARP1, PARP2") with free-text mechanism
    descriptions ("Microtubule destabiliser", "DNA crosslinker") and informal
    names that are not HGNC symbols ("MEK1", "BCL-XL", "VEGFR"). Validating each
    token against the actual symbol vocabulary is what keeps the two apart --
    anything unrecognised is reported rather than silently dropped, so a low
    recovery rate can't be mistaken for a model failure when it is really a
    naming mismatch.
    """
    if not putative_target or not putative_target.strip():
        return [], []

    recognised: list[str] = []
    unmatched: list[str] = []
    for raw_token in putative_target.split(","):
        token = raw_token.strip()
        if not token:
            continue
        if token.upper() in known_symbols:
            recognised.append(token.upper())
        else:
            unmatched.append(token)
    return recognised, unmatched


def target_recovery(
    top_genes: Sequence[dict],
    putative_target: str,
    known_symbols: set[str],
) -> dict:
    """Check whether the drug's annotated target appears among the top genes.

    This is the project's cheapest independent validation: the target
    annotation comes from GDSC and was never shown to the model, so recovering
    it is evidence the attribution reflects real pharmacology rather than an
    artefact of the graph's topology.
    """
    target_genes, unmatched = parse_target_symbols(putative_target, known_symbols)

    rank_by_symbol = {row["gene_symbol"].upper(): i + 1 for i, row in enumerate(top_genes)}
    recovered = [
        {"gene_symbol": symbol, "rank": rank_by_symbol[symbol]}
        for symbol in target_genes
        if symbol in rank_by_symbol
    ]

    return {
        "putative_target": putative_target,
        "target_genes": target_genes,
        "unmatched_tokens": unmatched,
        "recovered": recovered,
        "checked": bool(target_genes),
    }
