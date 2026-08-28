"""leave_drugs_out_gnn.py — does the graph help on *unseen drugs*?

`leave_drugs_out.py` benchmarks RF / MLP / CrossAttention on the drug-grouped
split and shows Morgan fingerprints beating one-hot by a wide margin (see
docs/leave_drugs_out_results.md). No GNN has ever been run through that split:
every GNN number in the project (docs/gnn_ablation_results.md, E12/E24/E56/E58
in docs/results.md) comes from the *cell-line*-grouped protocol, which says
nothing about generalization to new compounds.

That is the gap this script fills. It runs the tracked `HeteroGNN`
(variant="gcn", with cell_line->protein mutation edges — the E12 configuration,
RMSE 1.3513 on the cell-line split) against `leave_drugs_out_split`, so ~15% of
compounds are never scored during training.

The question it answers: does message passing over the PPI network and
drug-target edges add anything *beyond* what the 2048-bit Morgan fingerprint
already gives a flat model? The fingerprint MLP reaches RMSE 1.8516 / R^2 0.463
on this split. If the GNN cannot beat that, the graph is not contributing to
the proposal's actual generalization claim and is only earning its place on the
cell-line axis.

There is no one-hot arm here, unlike `leave_drugs_out.py`: the graph's drug
nodes carry fingerprint features by construction (03_graph_construction.py), so
this is inherently a fingerprint-only experiment — the same caveat
`gnn_baseline.py` notes.

CAVEAT — transductive vs inductive. This is full-batch message passing over one
fixed graph, so a held-out drug's node still participates in message passing
during training: its fingerprint and its drug->protein target edges are part of
the graph the whole time, even though its IC50 labels never enter the loss. The
MLP/RF baselines are strictly inductive by comparison — they see a held-out
drug's feature vector only at test time. The GNN therefore has a mild
structural advantage here, and a like-for-like inductive comparison would need
the held-out drug nodes masked out of the graph during training. Read the
comparison with that in mind.

Run from the repository root, after 03_graph_construction.py and
04_link_cell_lines.py:
    python "experiments/Leave Drugs Out/leave_drugs_out_gnn.py"
"""

from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, TensorDataset
from torch_geometric.data import HeteroData

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.data.experiment_utils import (  # noqa: E402
    COL_CELL_LINE,
    COL_DRUG,
    COL_TARGET,
    compute_shared_threshold,
    evaluate,
    leave_drugs_out_split,
    load_targets,
    mean_only_floor,
    peak_rss_gb,
    print_metric_block,
)
from src.models.hetero_gnn import HeteroGNN  # noqa: E402

RESULTS_CSV = Path("experiments/Leave Drugs Out/leave_drugs_out_gnn_results.csv")
GRAPH_PATH = Path("src/graph/hetero_graph.pt")
MUTATION_EDGE = ("cell_line", "has_mutation", "protein")


def load_module(name: str, relative_path: str):
    """Load a sibling experiment script by path (they live in spaced dirs)."""
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# Reuse the tracked GNN's hyperparameters and train/predict loop verbatim so
# this differs from the E12 run in exactly one respect: the split.
gnn_mod = load_module("gnn_baseline", "experiments/GNN Ablation/gnn_baseline.py")


def load_graph_and_drug_groups():
    """Like gnn_baseline.load_graph_and_pairs, but also returns drug_id groups."""
    if not GRAPH_PATH.exists():
        raise FileNotFoundError(f"Missing {GRAPH_PATH}. Run src/data/03_graph_construction.py first.")

    data: HeteroData = torch.load(GRAPH_PATH, weights_only=False)
    cell_ids = list(data["cell_line"].node_ids)
    drug_ids = list(data["drug"].node_ids)
    print(f"[graph] cell_line={len(cell_ids)}  drug={len(drug_ids)}  protein={len(data['protein'].node_ids)}")
    for et in data.edge_types:
        print(f"[graph] edge {et}: {data[et].edge_index.shape[1]}")

    y_df = load_targets(pd.Index(cell_ids))

    drug_to_index = {d: i for i, d in enumerate(drug_ids)}
    cell_to_index = {c: i for i, c in enumerate(cell_ids)}
    before = len(y_df)
    y_df = y_df[y_df[COL_DRUG].isin(drug_to_index)].reset_index(drop=True)
    print(f"[pairs] {before - len(y_df)} pairs dropped (drug absent from graph) -> {len(y_df)} remain")

    cell_index = y_df[COL_CELL_LINE].map(cell_to_index).to_numpy(dtype=np.int64)
    drug_index = y_df[COL_DRUG].map(drug_to_index).to_numpy(dtype=np.int64)
    target = y_df[COL_TARGET].to_numpy(dtype=np.float32)
    # The grouping key: drug identity, so the split holds out whole compounds.
    drug_groups = y_df[COL_DRUG].to_numpy()

    return data, cell_index, drug_index, target, drug_groups


def describe_split(drug_groups, cell_index, train_idx, val_idx, test_idx) -> dict:
    """Confirm the held-out drugs really are unseen during training."""
    train_drugs, test_drugs = set(drug_groups[train_idx]), set(drug_groups[test_idx])
    overlap = train_drugs & test_drugs
    print(f"[drugs]   train={len(train_drugs)}  val={len(set(drug_groups[val_idx]))}  "
          f"test={len(test_drugs)}  overlap={len(overlap)}")
    assert not overlap, f"drug leaked across splits: {sorted(overlap)[:5]}"
    # Cell lines DO overlap here by design -- we are holding out drugs, not cells.
    print(f"[cells]   train={len(set(cell_index[train_idx]))}  "
          f"test={len(set(cell_index[test_idx]))} (overlap expected)")
    return {
        "train_drugs": len(train_drugs),
        "test_drugs": len(test_drugs),
        "drug_overlap": len(overlap),
    }


def main() -> None:
    t_start = time.perf_counter()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    threshold = compute_shared_threshold()

    data, cell_index, drug_index, y, drug_groups = load_graph_and_drug_groups()

    print("\n" + "#" * 82)
    print("# LEAVE-DRUGS-OUT | GNN-GCN | drug=fingerprint | +mutation edges")
    print("#" * 82)

    torch.manual_seed(gnn_mod.TORCH_SEED)
    t0 = time.perf_counter()

    x_dict = {nt: data[nt].x.to(device) for nt in data.node_types}
    edge_dict = {et: data[et].edge_index.to(device) for et in data.edge_types}
    if MUTATION_EDGE not in edge_dict:
        raise RuntimeError(f"{MUTATION_EDGE} missing — run src/data/04_link_cell_lines.py first.")
    print(f"[edges] using {len(edge_dict)} edge types: {list(edge_dict.keys())}")

    train_idx, val_idx, test_idx = leave_drugs_out_split(drug_groups)
    info = describe_split(drug_groups, cell_index, train_idx, val_idx, test_idx)
    t_prep = time.perf_counter() - t0

    ds = TensorDataset(
        torch.from_numpy(cell_index[train_idx]),
        torch.from_numpy(drug_index[train_idx]),
        torch.from_numpy(y[train_idx]),
    )
    loader = DataLoader(ds, batch_size=gnn_mod.BATCH_SIZE, shuffle=True, drop_last=True)

    model = HeteroGNN(
        metadata=(list(data.node_types), list(edge_dict.keys())),
        cell_line_dim=data["cell_line"].x.shape[1],
        drug_dim=data["drug"].x.shape[1],
        num_proteins=data["protein"].x.shape[0],
        variant="gcn",
        hidden_dim=gnn_mod.HIDDEN_DIM,
        num_layers=gnn_mod.NUM_LAYERS,
        heads=gnn_mod.HEADS,
        dropout=gnn_mod.DROPOUT,
        head_hidden_dims=gnn_mod.HEAD_HIDDEN_DIMS,
        head_dropout=gnn_mod.HEAD_DROPOUT,
    ).to(device)

    # Lazy (-1, -1) conv dims materialize on the first forward pass.
    model.eval()  # BatchNorm in the head rejects a 1-sample batch in train mode
    with torch.no_grad():
        model(
            x_dict, edge_dict,
            torch.from_numpy(cell_index[train_idx[:2]]).to(device),
            torch.from_numpy(drug_index[train_idx[:2]]).to(device),
        )
    n_params = sum(p.numel() for p in model.parameters())

    t1 = time.perf_counter()
    model, best_epoch = gnn_mod.train(
        model, x_dict, edge_dict, loader,
        cell_index[val_idx], drug_index[val_idx], y[val_idx], device,
    )
    t_fit = time.perf_counter() - t1

    val = evaluate(y[val_idx], gnn_mod.predict(model, x_dict, edge_dict, cell_index[val_idx], drug_index[val_idx], device), threshold)
    test = evaluate(y[test_idx], gnn_mod.predict(model, x_dict, edge_dict, cell_index[test_idx], drug_index[test_idx], device), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])

    print_metric_block("LDO | GNN-GCN | fingerprint", val, test, floor)
    print(f"n_pairs={len(y)}  params={n_params:,}  best_epoch={best_epoch}  "
          f"prep={t_prep:.1f}s  fit={t_fit:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    row = {
        "model": "GNN-GCN",
        "omics": "GE+Mut_CNV+Proteomics",
        "drug_rep": "fingerprint",
        "mutation_edges": True,
        "n_pairs": len(y),
        "mean_only_rmse": floor,
        **info,
        **{f"val_{k}": v for k, v in val.items()},
        **{f"test_{k}": v for k, v in test.items()},
        "params": n_params,
        "best_epoch": best_epoch,
        "fit_seconds": t_fit,
    }

    df = pd.DataFrame([row])
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"LEAVE-DRUGS-OUT GNN COMPLETE — 1 run in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    print(df[["model", "drug_rep", "test_drugs", "mean_only_rmse",
              "test_rmse", "test_pcc", "test_r2"]].round(4).to_string(index=False))
    print("\n--- context: fingerprint arms from leave_drugs_out.py (docs/leave_drugs_out_results.md) ---")
    print("  RF             fingerprint  RMSE 1.8719  R2 0.451")
    print("  MLP            fingerprint  RMSE 1.8516  R2 0.463")
    print("  CrossAttention fingerprint  RMSE 1.9200  R2 0.422")
    print(f"  GNN-GCN        fingerprint  RMSE {row['test_rmse']:.4f}  R2 {row['test_r2']:.3f}  <- this run")
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
