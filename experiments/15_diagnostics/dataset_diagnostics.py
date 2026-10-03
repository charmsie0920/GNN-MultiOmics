"""dataset_diagnostics.py — why the models plateau at RMSE 1.30-1.33.

Every number quoted in docs/15_diagnostics_results.md is produced here, so the
claims in docs/results.md ("drug identity alone accounts for 71% of the
reducible error", the 1.4889 per-drug floor) have an auditable source rather
than living only in a chat transcript or a hand-computed constant.

Three independent diagnostics, none of which train anything:

1. **Target variance decomposition.** How much of ln(IC50) is explained by
   compound identity alone, and what the trivial no-ML baselines score on the
   project's own cell-line-grouped test fold.

2. **Graph connectivity audit.** How much of the heterogeneous graph is
   actually reachable from a labelled (cell_line, drug) pair, and how many
   label pairs have a node that is isolated -- for which the GNN degenerates
   to an MLP.

3. **Receptive-field check.** Whether the PPI network can transmit any
   label-derived signal into the cell_line and drug embeddings at the model's
   configured depth. At 2 layers it cannot: PPI neighbours contribute only
   their layer-0 embeddings, which are free parameters carrying no data.

Run from the repository root, after 03_graph_construction.py and
04_link_cell_lines.py:
    python experiments/15_diagnostics/dataset_diagnostics.py
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.data.experiment_utils import (  # noqa: E402
    COL_CELL_LINE,
    COL_DRUG,
    COL_TARGET,
    grouped_split,
    load_targets,
    mean_only_floor,
    per_drug_mean_floor,
)

GRAPH_PATH = Path("src/graph/hetero_graph.pt")
RESULTS_CSV = Path("experiments/15_diagnostics/dataset_diagnostics_results.csv")

# The depth HeteroIC50GNN is configured with. The receptive-field check below
# reports what the PPI edges can and cannot deliver at this value.
NUM_LAYERS = 2

PPI = ("protein", "interacts_with", "protein")
TARGETS = ("drug", "targets", "protein")
MUTATION = ("cell_line", "has_mutation", "protein")


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def load_pairs(data):
    """The exact 111,799-pair population the matrix runs use."""
    cell_ids = list(data["cell_line"].node_ids)
    drug_ids = [str(d) for d in data["drug"].node_ids]
    y_df = load_targets(pd.Index(cell_ids))
    y_df = y_df[y_df[COL_DRUG].astype(str).isin(set(drug_ids))].reset_index(drop=True)

    cell_to_index = {c: i for i, c in enumerate(cell_ids)}
    drug_to_index = {d: i for i, d in enumerate(drug_ids)}
    return (
        y_df,
        y_df[COL_CELL_LINE].map(cell_to_index).to_numpy(np.int64),
        y_df[COL_DRUG].astype(str).map(drug_to_index).to_numpy(np.int64),
        y_df[COL_TARGET].to_numpy(np.float64),
        y_df[COL_CELL_LINE].to_numpy(),
    )


def diagnose_targets(y, drug_index, groups) -> dict:
    section("1. TARGET VARIANCE — how much of ln(IC50) is free?")

    total = float(np.var(y))
    df = pd.DataFrame({"y": y, "d": drug_index, "c": groups})
    within_drug = float(np.var(y - df.groupby("d")["y"].transform("mean").to_numpy()))
    within_cell = float(np.var(y - df.groupby("c")["y"].transform("mean").to_numpy()))

    print(f"  total variance                       {total:8.4f}   (sd {np.sqrt(total):.4f})")
    print(f"  remaining after removing DRUG mean   {within_drug:8.4f}   (sd {np.sqrt(within_drug):.4f})")
    print(f"  remaining after removing CELL mean   {within_cell:8.4f}   (sd {np.sqrt(within_cell):.4f})")
    drug_share = 100 * (1 - within_drug / total)
    cell_share = 100 * (1 - within_cell / total)
    print(f"\n  -> {drug_share:.1f}% of all variation is just WHICH DRUG it is")
    print(f"  -> {cell_share:.1f}% of all variation is just WHICH CELL LINE it is")

    train_idx, val_idx, test_idx = grouped_split(groups)
    global_floor = mean_only_floor(y[train_idx], y[test_idx])
    drug_floor = per_drug_mean_floor(
        y[train_idx], drug_index[train_idx], y[test_idx], drug_index[test_idx]
    )

    print("\n  Baselines on the model's own grouped TEST fold (no ML of any kind):")
    print(f"    predict global training mean     {global_floor:.4f}   <- the floor docs/results.md used to quote")
    print(f"    predict per-drug training mean   {drug_floor:.4f}   <- the floor that actually matters")
    print(f"\n  A model scoring 1.3301 therefore buys {drug_floor - 1.3301:+.4f} RMSE "
          f"({100 * (drug_floor - 1.3301) / drug_floor:+.1f}%) over a lookup table.")

    return {
        "total_variance": total,
        "within_drug_variance": within_drug,
        "within_cell_variance": within_cell,
        "drug_identity_share_pct": drug_share,
        "cell_identity_share_pct": cell_share,
        "global_mean_floor": global_floor,
        "per_drug_mean_floor": drug_floor,
    }


def diagnose_connectivity(data, cell_index, drug_index) -> dict:
    section("2. GRAPH CONNECTIVITY — how much of the graph is reachable?")

    n_cell = data["cell_line"].num_nodes
    n_drug = data["drug"].num_nodes
    n_prot = data["protein"].num_nodes

    tgt_e = data[TARGETS].edge_index.numpy()
    mut_e = data[MUTATION].edge_index.numpy()
    ppi_e = data[PPI].edge_index.numpy()

    drug_deg = np.bincount(tgt_e[0], minlength=n_drug)
    cell_deg = np.bincount(mut_e[0], minlength=n_cell)
    iso_drugs = int((drug_deg == 0).sum())
    iso_cells = int((cell_deg == 0).sum())

    print(f"  drugs with no target edge       {iso_drugs:5d} / {n_drug}  ({100*iso_drugs/n_drug:.1f}%)"
          f"   degree mean {drug_deg.mean():.2f}, max {drug_deg.max()}")
    print(f"  cell lines with no mutation edge{iso_cells:5d} / {n_cell}  ({100*iso_cells/n_cell:.1f}%)"
          f"   degree mean {cell_deg.mean():.2f}, max {cell_deg.max()}")

    iso_drug_set = set(np.flatnonzero(drug_deg == 0).tolist())
    iso_cell_set = set(np.flatnonzero(cell_deg == 0).tolist())
    pairs_iso_drug = int(np.isin(drug_index, list(iso_drug_set)).sum())
    pairs_iso_cell = int(np.isin(cell_index, list(iso_cell_set)).sum())
    n_pairs = len(drug_index)

    print(f"\n  label pairs whose DRUG is isolated      {pairs_iso_drug:6,d} / {n_pairs:,}  "
          f"({100*pairs_iso_drug/n_pairs:.1f}%)")
    print(f"  label pairs whose CELL LINE is isolated {pairs_iso_cell:6,d} / {n_pairs:,}  "
          f"({100*pairs_iso_cell/n_pairs:.1f}%)")
    print("  -> for those pairs SAGEConv aggregates an empty neighbourhood and the")
    print("     model reduces to an MLP on the node's own features.")

    touched_drug = set(np.unique(tgt_e[1]).tolist())
    touched_cell = set(np.unique(mut_e[1]).tolist())
    attached = touched_drug | touched_cell
    orphan = n_prot - len(attached)

    ppi_deg = np.bincount(ppi_e[0], minlength=n_prot)
    print(f"\n  proteins that are a drug target         {len(touched_drug):6d}")
    print(f"  proteins mutated in some cell line      {len(touched_cell):6d}")
    print(f"  proteins attached to NOTHING            {orphan:6d} / {n_prot}  ({100*orphan/n_prot:.1f}%)")
    print(f"  -> {orphan * 128:,} embedding parameters never within one hop of a label.")
    print(f"\n  PPI degree: all {ppi_deg.mean():.1f} | drug targets {ppi_deg[list(touched_drug)].mean():.1f}"
          f" | mutated {ppi_deg[list(touched_cell)].mean():.1f}")
    print("  -> SAGEConv MEAN-aggregates neighbours, so informative ones are divided")
    print("     into an average dominated by unattached proteins.")

    return {
        "isolated_drugs": iso_drugs,
        "isolated_cell_lines": iso_cells,
        "pairs_with_isolated_drug": pairs_iso_drug,
        "pairs_with_isolated_cell": pairs_iso_cell,
        "pct_pairs_isolated_drug": 100 * pairs_iso_drug / n_pairs,
        "proteins_drug_targets": len(touched_drug),
        "proteins_mutated": len(touched_cell),
        "proteins_unattached": orphan,
        "pct_proteins_unattached": 100 * orphan / n_prot,
        "wasted_embedding_params": orphan * 128,
        "ppi_degree_mean": float(ppi_deg.mean()),
        "ppi_degree_mean_targets": float(ppi_deg[list(touched_drug)].mean()),
        "ppi_degree_mean_mutated": float(ppi_deg[list(touched_cell)].mean()),
    }


def diagnose_receptive_field(data) -> dict:
    section(f"3. RECEPTIVE FIELD — what reaches the prediction at {NUM_LAYERS} layers?")

    tgt_e = data[TARGETS].edge_index.numpy()
    mut_e = data[MUTATION].edge_index.numpy()
    ppi_e = data[PPI].edge_index.numpy()

    nbr = defaultdict(set)
    for a, b in zip(ppi_e[0], ppi_e[1]):
        nbr[a].add(b)
    target_set = set(np.unique(tgt_e[1]).tolist())

    cell_muts = defaultdict(set)
    for c, p in zip(mut_e[0], mut_e[1]):
        cell_muts[c].add(p)

    n_cell = data["cell_line"].num_nodes
    direct = sum(1 for c in range(n_cell) if cell_muts.get(c, set()) & target_set)
    via_ppi = 0
    for c in range(n_cell):
        two_hop: set = set()
        for m in cell_muts.get(c, set()):
            two_hop |= nbr[m]
        if two_hop & target_set:
            via_ppi += 1

    print(f"  cell lines whose mutated protein IS itself a drug target   {direct:4d} / {n_cell}")
    print(f"  cell lines reaching a drug target within 1 PPI hop         {via_ppi:4d} / {n_cell}")

    print("\n  Message-passing trace for h_cell:")
    print("    conv 1:  protein <- its PPI neighbours' LAYER-0 embeddings (no data),")
    print("                     <- drugs targeting it      (real fingerprint data)")
    print("                     <- cell lines mutating it  (real omics data)")
    print("    conv 2:  cell_line <- protein (layer 1)")
    print("\n  So after 2 layers h_cell contains drugs that target its OWN mutated")
    print("  proteins -- the direct path. Reaching a drug that targets a PPI")
    print("  NEIGHBOUR of those proteins needs:")
    print("       drug -> protein_B (conv 1), protein_B -> protein_A (conv 2),")
    print("       protein_A -> cell_line (conv 3)")
    print(f"\n  The model has {NUM_LAYERS} layers, so the {ppi_e.shape[1]:,} PPI edges "
          f"({100*ppi_e.shape[1]/(ppi_e.shape[1]+tgt_e.shape[1]+mut_e.shape[1]):.1f}% of all edges)")
    print("  transmit NO label-derived signal into h_cell or h_drug. They contribute")
    print("  only free parameters -- consistent with the measured ~0.038 RMSE cost of")
    print("  the PPI branch in docs/09 and docs/13.")

    return {
        "cells_with_direct_target_mutation": direct,
        "cells_reaching_target_via_1_ppi_hop": via_ppi,
        "num_layers": NUM_LAYERS,
        "ppi_edges": int(ppi_e.shape[1]),
        "ppi_signal_reaches_readout": False,
    }


def main() -> None:
    if not GRAPH_PATH.exists():
        raise FileNotFoundError(f"Missing {GRAPH_PATH}. Run src/data/03_graph_construction.py first.")
    data = torch.load(GRAPH_PATH, weights_only=False)

    _, cell_index, drug_index, y, groups = load_pairs(data)
    print(f"\n[population] {len(y):,} labelled pairs, "
          f"{data['cell_line'].num_nodes} cell lines, {data['drug'].num_nodes} drugs, "
          f"{data['protein'].num_nodes:,} proteins")

    row = {"n_pairs": len(y)}
    row.update(diagnose_targets(y, drug_index, groups))
    row.update(diagnose_connectivity(data, cell_index, drug_index))
    row.update(diagnose_receptive_field(data))

    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(RESULTS_CSV, index=False)
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
