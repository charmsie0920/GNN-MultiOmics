"""Add (cell_line, has_mutation, protein) edges so cell lines aren't a disconnected component.

`03_graph_construction.py` builds cell_line, drug, and protein nodes plus
STRING PPI and drug-target edges -- but nothing connects cell_line to
anything else, so a GNN over that graph cannot propagate any PPI or
drug-target signal into a cell-line embedding. This script closes that gap
by linking each cell line to the proteins whose genes carry a *cancer driver*
mutation in that cell line.

Why driver mutations specifically: the mutation table has ~6,300 mutations
per cell line, the overwhelming majority intronic/non-coding. Linking all of
them would add millions of mostly-noise edges that swamp the 474K
high-confidence PPI edges. The `cancer_driver` flag (~0.3% of rows) is the
biologically motivated sparse subset, and produces exactly the multi-hop path
the architecture is meant to exploit:

    cell_line --has_mutation--> protein --interacts_with--> protein <--targets-- drug

Crucially these edges are derived from mutation status only, never from IC50,
so they introduce no label leakage into the drug-response task.

Run from the repository root, after 03_graph_construction.py:
    python src/data/04_link_cell_lines.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import pandas as pd
import torch
from torch_geometric.data import HeteroData

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.graph_utils import build_gene_symbol_to_ensp  # noqa: E402

MUTATIONS_PATH = Path("data/processed/aligned/mutations_aligned.csv")
GRAPH_PATH = Path("src/graph/hetero_graph.pt")

MUTATION_CHUNKSIZE = 1_000_000
EDGE_TYPE = ("cell_line", "has_mutation", "protein")


def collect_driver_mutations(
    cell_line_ids: Set[str], driver_only: bool
) -> Tuple[Dict[str, Set[str]], int, Set[str]]:
    """Stream the mutation table, returning {cell_line_id: {gene_symbol, ...}}.

    The file is 2.3GB / 14.4M rows, so it's read in chunks with only the three
    needed columns -- never materialized whole, matching the memory discipline
    already used for the STRING PPI file in 03_graph_construction.py.
    """
    mutations_by_cell: Dict[str, Set[str]] = {}
    all_symbols: Set[str] = set()
    total_rows = 0

    reader = pd.read_csv(
        MUTATIONS_PATH,
        usecols=["standard_model_id", "gene_symbol", "cancer_driver"],
        dtype=str,
        chunksize=MUTATION_CHUNKSIZE,
    )
    for chunk in reader:
        total_rows += len(chunk)
        chunk = chunk[chunk["standard_model_id"].isin(cell_line_ids)]
        if driver_only:
            chunk = chunk[chunk["cancer_driver"].str.lower() == "t"]
        chunk = chunk.dropna(subset=["gene_symbol"])
        for model_id, symbol in zip(chunk["standard_model_id"], chunk["gene_symbol"]):
            if not symbol:
                continue
            mutations_by_cell.setdefault(model_id, set()).add(symbol)
            all_symbols.add(symbol)

    return mutations_by_cell, total_rows, all_symbols


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=GRAPH_PATH)
    parser.add_argument(
        "--all-mutations",
        action="store_true",
        help="Link every mutated gene, not just cancer drivers (much denser, noisier).",
    )
    args = parser.parse_args(argv)

    if not args.graph.exists():
        raise FileNotFoundError(f"Missing {args.graph}. Run src/data/03_graph_construction.py first.")
    if not MUTATIONS_PATH.exists():
        raise FileNotFoundError(f"Missing {MUTATIONS_PATH}. Run src/data/ingest_and_align.py first.")

    data: HeteroData = torch.load(args.graph, weights_only=False)
    cell_line_ids: List[str] = list(data["cell_line"].node_ids)
    protein_ids: List[str] = list(data["protein"].node_ids)
    protein_to_index = {pid: i for i, pid in enumerate(protein_ids)}

    print(f"[graph] loaded {args.graph}: {len(cell_line_ids)} cell lines, {len(protein_ids)} proteins")

    driver_only = not args.all_mutations
    mutations_by_cell, total_rows, all_symbols = collect_driver_mutations(
        set(cell_line_ids), driver_only
    )
    print(
        f"[mutations] scanned {total_rows} rows -> "
        f"{sum(len(v) for v in mutations_by_cell.values())} (cell_line, gene) pairs "
        f"across {len(mutations_by_cell)} cell lines, {len(all_symbols)} distinct genes "
        f"(driver_only={driver_only})"
    )

    symbol_to_ensp = build_gene_symbol_to_ensp(all_symbols)
    unresolved = all_symbols - set(symbol_to_ensp)
    print(f"[aliases] {len(symbol_to_ensp)} / {len(all_symbols)} gene symbols resolved to an ENSP")

    src_indices: List[int] = []
    dst_indices: List[int] = []
    new_proteins = 0
    for cell_index, cell_id in enumerate(cell_line_ids):
        for symbol in mutations_by_cell.get(cell_id, ()):
            ensp = symbol_to_ensp.get(symbol)
            if ensp is None:
                continue
            protein_index = protein_to_index.get(ensp)
            if protein_index is None:
                protein_index = len(protein_ids)
                protein_to_index[ensp] = protein_index
                protein_ids.append(ensp)
                new_proteins += 1
            src_indices.append(cell_index)
            dst_indices.append(protein_index)

    if new_proteins:
        # Newly referenced proteins need matching placeholder feature rows.
        embed_dim = data["protein"].x.shape[1]
        data["protein"].x = torch.cat(
            [data["protein"].x, torch.zeros((new_proteins, embed_dim), dtype=data["protein"].x.dtype)],
            dim=0,
        )
        data["protein"].node_ids = protein_ids
        print(f"[protein] added {new_proteins} protein nodes referenced only by mutations")

    edge_index = torch.tensor([src_indices, dst_indices], dtype=torch.long)
    data[EDGE_TYPE].edge_index = edge_index

    if edge_index.numel() > 0:
        assert edge_index[0].max().item() < data["cell_line"].x.shape[0], "cell_line index out of range"
        assert edge_index[1].max().item() < data["protein"].x.shape[0], "protein index out of range"

    linked_cells = len(set(src_indices))
    print(f"[has_mutation] {edge_index.shape[1]} edges linking {linked_cells} / {len(cell_line_ids)} cell lines")
    if linked_cells < len(cell_line_ids):
        print(f"[has_mutation] WARNING: {len(cell_line_ids) - linked_cells} cell lines remain unlinked")

    print("\n=== Updated HeteroData ===")
    for node_type in data.node_types:
        print(f"  node '{node_type}': num_nodes={data[node_type].x.shape[0]}, dim={data[node_type].x.shape[1]}")
    for edge_type in data.edge_types:
        print(f"  edge {edge_type}: num_edges={data[edge_type].edge_index.shape[1]}")

    torch.save(data, args.graph)
    print(f"\nSaved -> {args.graph}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
