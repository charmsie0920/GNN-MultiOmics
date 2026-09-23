"""Pathway over-representation analysis for a ranked gene list, via Enrichr.

Model-agnostic: takes gene symbols, returns enriched terms. Results are cached
on disk keyed by the gene list, because the same (cell line, drug) pair is
re-requested every time the user clicks back to a row and Enrichr is a public
service that should not be hit repeatedly for an identical query.

Nothing here raises into the caller. A dead or blocked network must leave the
pathway panel empty with an explanation, never fail the prediction the user
actually asked for.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Sequence

CACHE_DIR = Path("data/processed/enrichment_cache")

DEFAULT_LIBRARIES = (
    "GO_Biological_Process_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
)

# Enrichr rejects very short lists as statistically meaningless, and an
# over-representation test on a handful of genes is not worth showing anyway.
MIN_GENES = 5
MAX_TERMS_PER_LIBRARY = 25


def _cache_key(gene_symbols: Sequence[str], library: str) -> str:
    digest = hashlib.sha1("|".join(sorted(gene_symbols)).encode("utf-8")).hexdigest()
    return f"{digest}-{library}"


def _read_cache(key: str) -> list[dict] | None:
    path = CACHE_DIR / f"{key}.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _write_cache(key: str, terms: list[dict]) -> None:
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        (CACHE_DIR / f"{key}.json").write_text(json.dumps(terms), encoding="utf-8")
    except OSError:
        # A cache that cannot be written is a slow path, not an error.
        pass


def _query_enrichr(gene_symbols: Sequence[str], library: str) -> list[dict]:
    import gseapy

    result = gseapy.enrichr(
        gene_list=list(gene_symbols),
        gene_sets=library,
        organism="human",
        outdir=None,  # keep it in memory; this must not litter the repo
        no_plot=True,
    )
    frame = result.results
    frame = frame.sort_values("Adjusted P-value").head(MAX_TERMS_PER_LIBRARY)

    return [
        {
            "term": str(row["Term"]),
            "library": library,
            "p_value": float(row["P-value"]),
            "adjusted_p_value": float(row["Adjusted P-value"]),
            "combined_score": float(row["Combined Score"]),
            "overlap": str(row["Overlap"]),
            "genes": str(row["Genes"]).split(";") if row["Genes"] else [],
        }
        for _, row in frame.iterrows()
    ]


def enrich(
    gene_symbols: Sequence[str],
    libraries: Sequence[str] = DEFAULT_LIBRARIES,
) -> dict:
    """Run over-representation analysis, returning terms plus a status.

    The status is part of the contract, not a side channel: the UI needs to
    distinguish "no pathways were significant" from "we could not reach
    Enrichr", and those look identical if only a term list comes back.
    """
    genes = [symbol for symbol in dict.fromkeys(gene_symbols) if symbol]
    if len(genes) < MIN_GENES:
        return {
            "terms": [],
            "status": "insufficient_genes",
            "message": f"Need at least {MIN_GENES} mapped genes to test for enrichment, got {len(genes)}.",
        }

    terms: list[dict] = []
    failures: list[str] = []
    for library in libraries:
        key = _cache_key(genes, library)
        cached = _read_cache(key)
        if cached is not None:
            terms.extend(cached)
            continue
        try:
            fetched = _query_enrichr(genes, library)
        except Exception as exc:  # noqa: BLE001 -- any failure degrades to an empty panel
            failures.append(f"{library}: {type(exc).__name__}")
            continue
        _write_cache(key, fetched)
        terms.extend(fetched)

    if failures and not terms:
        return {
            "terms": [],
            "status": "unavailable",
            "message": "Could not reach Enrichr. Pathway enrichment needs an internet connection. "
            + "; ".join(failures),
        }

    terms.sort(key=lambda t: t["adjusted_p_value"])
    return {
        "terms": terms,
        "status": "partial" if failures else "ok",
        "message": ("Some libraries failed: " + "; ".join(failures)) if failures else "",
    }
