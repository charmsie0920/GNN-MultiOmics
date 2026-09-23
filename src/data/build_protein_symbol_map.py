"""Build the ENSP -> gene symbol lookup the interpretation layer displays.

`graph_utils.build_gene_symbol_to_ensp` resolves the *forward* direction
(symbol -> ENSP) during graph construction, filtered to the symbols GDSC
happens to reference and taking whichever alias row appears first. That is
fine for edge building, but the reverse direction needs a canonical symbol per
protein: a gene importance panel showing `9606.ENSP00000269305` instead of
`TP53` is unreadable, and first-match would just as happily yield a PDB
accession.

So this resolves each protein against an explicit source priority and writes a
small CSV keyed by the graph's own protein node ids. The 19 MB alias file is
read once here; nothing at runtime touches it again.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import Optional, Sequence

import torch

ALIASES_PATH = Path("data/raw/string/9606.protein.aliases.v12.0.txt.gz")
GRAPH_PATH = Path("src/graph/hetero_graph.pt")
OUTPUT_PATH = Path("data/processed/protein_symbol_map.csv")

# Most authoritative first. STRING carries ~90 alias sources per release and
# most are accessions (PDB, RefSeq, UniParc), not symbols. `Ensembl_HGNC_symbol`
# is the approved HGNC symbol and covers ~19.2k of the ~19.7k human genes; the
# rest are fallbacks for proteins HGNC has not assigned a current symbol to.
SOURCE_PRIORITY = (
    "Ensembl_HGNC_symbol",
    "Ensembl_HGNC",
    "UniProt_GN_Name",
    "Ensembl_HGNC_prev_symbol",
    "Ensembl_HGNC_alias_symbol",
)


def load_graph_protein_ids(graph_path: Path) -> list[str]:
    """Read the protein node ids straight off the graph, preserving node order."""
    graph = torch.load(graph_path, weights_only=False)
    return [str(protein_id) for protein_id in graph["protein"].node_ids]


def build_symbol_map(protein_ids: Sequence[str], aliases_path: Path) -> dict[str, str]:
    """Resolve each protein id to its highest-priority available gene symbol.

    Streams the gzip rather than loading it: the file is ~40M rows and only the
    graph's proteins matter.
    """
    wanted = set(protein_ids)
    rank = {source: i for i, source in enumerate(SOURCE_PRIORITY)}

    # protein_id -> (source_rank, symbol); a row only wins if it outranks what
    # is already held, so input order does not affect the result.
    best: dict[str, tuple[int, str]] = {}
    with gzip.open(aliases_path, "rt", encoding="utf-8") as handle:
        next(handle)  # header: #string_protein_id, alias, source
        for line in handle:
            protein_id, alias, source = line.rstrip("\n").split("\t")
            if protein_id not in wanted:
                continue
            source_rank = rank.get(source)
            if source_rank is None:
                continue
            current = best.get(protein_id)
            if current is None or source_rank < current[0]:
                best[protein_id] = (source_rank, alias)

    return {protein_id: symbol for protein_id, (_, symbol) in best.items()}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aliases", type=Path, default=ALIASES_PATH)
    parser.add_argument("--graph", type=Path, default=GRAPH_PATH)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args(argv)

    if not args.aliases.exists():
        raise SystemExit(
            f"Alias file not found at {args.aliases}. Run src/data/fetch_string_aliases.py first."
        )

    protein_ids = load_graph_protein_ids(args.graph)
    print(f"Graph protein nodes: {len(protein_ids)}")

    symbol_map = build_symbol_map(protein_ids, args.aliases)
    resolved = len(symbol_map)
    print(f"Resolved to a gene symbol: {resolved} ({resolved / len(protein_ids):.1%})")

    unresolved = [protein_id for protein_id in protein_ids if protein_id not in symbol_map]
    if unresolved:
        print(f"Unresolved (kept, symbol left blank): {len(unresolved)}")
        print(f"  examples: {', '.join(unresolved[:5])}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["protein_id", "gene_symbol"])
        # Every graph protein gets a row, symbol blank when unresolved, so a
        # consumer can distinguish "not in the graph" from "no symbol known".
        for protein_id in protein_ids:
            writer.writerow([protein_id, symbol_map.get(protein_id, "")])

    print(f"Wrote {args.output} ({args.output.stat().st_size / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
