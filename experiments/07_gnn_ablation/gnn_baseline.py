"""gnn_baseline.py — Heterogeneous GNN over the cell_line/drug/protein graph.

Phase 3 of the experiment matrix (docs/plan/experiment_matrix_plan.md), and
the first model in the project to use relational structure at all: message
passing over STRING PPI edges, drug-target edges, and (optionally) the
cell_line->protein driver-mutation edges added by 04_link_cell_lines.py.

Four runs: {GCN, GAT} x {with, without has_mutation edges}. The with/without
axis is the direct analogue of the MoGraphDRP paper's "without graph"
ablation -- it answers whether connecting cell lines into the graph actually
buys anything, or whether the PPI subgraph is decorative.

Drug features here are the 2048-bit Morgan fingerprints already stored on the
graph's drug nodes, so this is inherently a fingerprint-mode experiment;
there is no one-hot variant.

The whole graph fits in memory (~20MB of tensors), so this does full-batch
message passing each step and mini-batches only the (cell_line, drug) label
pairs -- no neighbor sampling needed.

Run from the repository root, after 03_graph_construction.py and
04_link_cell_lines.py:
    python "experiments/07_gnn_ablation/gnn_baseline.py"
"""

from __future__ import annotations

import copy
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset
from torch_geometric.data import HeteroData

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.experiment_utils import (  # noqa: E402
    COL_CELL_LINE,
    COL_DRUG,
    COL_TARGET,
    compute_shared_threshold,
    evaluate,
    grouped_split,
    load_targets,
    mean_only_floor,
    peak_rss_gb,
    print_metric_block,
)
from src.models.hetero_gnn import HeteroGNN  # noqa: E402

GRAPH_PATH = Path("src/graph/hetero_graph.pt")
RESULTS_CSV = Path("experiments/07_gnn_ablation/gnn_results.csv")

TORCH_SEED = 42
HIDDEN_DIM = 128
NUM_LAYERS = 2
HEADS = 4
DROPOUT = 0.2
HEAD_HIDDEN_DIMS = (256, 128)
HEAD_DROPOUT = 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-5
BATCH_SIZE = 1024
MAX_EPOCHS = 200
PATIENCE = 15
LR_PATIENCE = 5

MUTATION_EDGE = ("cell_line", "has_mutation", "protein")


def load_graph_and_pairs(threshold: float):
    """Load the hetero graph and join GDSC IC50 labels onto its node index space."""
    if not GRAPH_PATH.exists():
        raise FileNotFoundError(f"Missing {GRAPH_PATH}. Run src/data/03_graph_construction.py first.")

    data: HeteroData = torch.load(GRAPH_PATH, weights_only=False)
    cell_ids = list(data["cell_line"].node_ids)
    drug_ids = list(data["drug"].node_ids)
    print(f"[graph] cell_line={len(cell_ids)}  drug={len(drug_ids)}  protein={len(data['protein'].node_ids)}")
    for et in data.edge_types:
        print(f"[graph] edge {et}: {data[et].edge_index.shape[1]}")

    y_df = load_targets(pd.Index(cell_ids))

    # Restrict labels to pairs whose drug survived SMILES parsing (drug nodes
    # only exist for parseable drugs), then map both sides onto node indices.
    drug_to_index = {d: i for i, d in enumerate(drug_ids)}
    cell_to_index = {c: i for i, c in enumerate(cell_ids)}
    before = len(y_df)
    y_df = y_df[y_df[COL_DRUG].isin(drug_to_index)].reset_index(drop=True)
    print(f"[pairs] {before - len(y_df)} pairs dropped (drug absent from graph) -> {len(y_df)} remain")

    cell_index = y_df[COL_CELL_LINE].map(cell_to_index).to_numpy(dtype=np.int64)
    drug_index = y_df[COL_DRUG].map(drug_to_index).to_numpy(dtype=np.int64)
    target = y_df[COL_TARGET].to_numpy(dtype=np.float32)
    groups = y_df[COL_CELL_LINE].to_numpy()

    return data, cell_index, drug_index, target, groups


@torch.no_grad()
def predict(model, x_dict, edge_dict, cell_index, drug_index, device) -> np.ndarray:
    model.eval()
    preds = []
    for start in range(0, len(cell_index), BATCH_SIZE):
        end = start + BATCH_SIZE
        ci = torch.from_numpy(cell_index[start:end]).to(device)
        di = torch.from_numpy(drug_index[start:end]).to(device)
        preds.append(model(x_dict, edge_dict, ci, di).cpu().numpy())
    return np.concatenate(preds)


def train(model, x_dict, edge_dict, loader, cell_va, drug_va, y_va, device):
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=LR_PATIENCE
    )
    best_state = copy.deepcopy(model.state_dict())
    best_rmse, best_epoch, stale = float("inf"), 0, 0

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        for ci, di, yb in loader:
            ci, di, yb = ci.to(device), di.to(device), yb.to(device)
            optimizer.zero_grad()
            loss = criterion(model(x_dict, edge_dict, ci, di), yb)
            loss.backward()
            optimizer.step()

        val_pred = predict(model, x_dict, edge_dict, cell_va, drug_va, device)
        val_rmse = float(np.sqrt(np.mean((y_va - val_pred) ** 2)))
        scheduler.step(val_rmse)

        if val_rmse < best_rmse - 1e-4:
            best_rmse, best_epoch, stale = val_rmse, epoch, 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            stale += 1

        if epoch % 5 == 0:
            print(f"  [epoch {epoch:>3}] val_rmse={val_rmse:.4f}  best={best_rmse:.4f} @ {best_epoch}")

        if stale >= PATIENCE:
            print(f"  [early stop] epoch {epoch}, best epoch {best_epoch}, val_rmse={best_rmse:.4f}")
            break

    model.load_state_dict(best_state)
    return model, best_epoch


def run_one(data, cell_index, drug_index, y, groups, variant, use_mutation_edges, threshold, device) -> dict:
    label = f"{variant.upper()} | mutation_edges={'yes' if use_mutation_edges else 'no'}"
    print("\n" + "#" * 82)
    print(f"# GNN | {label}")
    print("#" * 82)

    torch.manual_seed(TORCH_SEED)
    t0 = time.perf_counter()

    x_dict = {nt: data[nt].x.to(device) for nt in data.node_types}
    edge_dict = {
        et: data[et].edge_index.to(device)
        for et in data.edge_types
        if use_mutation_edges or et != MUTATION_EDGE
    }
    print(f"[edges] using {len(edge_dict)} edge types: {list(edge_dict.keys())}")

    train_idx, val_idx, test_idx = grouped_split(groups)
    t_prep = time.perf_counter() - t0

    ds = TensorDataset(
        torch.from_numpy(cell_index[train_idx]),
        torch.from_numpy(drug_index[train_idx]),
        torch.from_numpy(y[train_idx]),
    )
    loader = DataLoader(ds, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)

    model = HeteroGNN(
        metadata=(list(data.node_types), list(edge_dict.keys())),
        cell_line_dim=data["cell_line"].x.shape[1],
        drug_dim=data["drug"].x.shape[1],
        num_proteins=data["protein"].x.shape[0],
        variant=variant,
        hidden_dim=HIDDEN_DIM,
        num_layers=NUM_LAYERS,
        heads=HEADS,
        dropout=DROPOUT,
        head_hidden_dims=HEAD_HIDDEN_DIMS,
        head_dropout=HEAD_DROPOUT,
    ).to(device)

    # SAGEConv/GATConv are built with (-1, -1) input dims, so their weights are
    # lazily created on the first forward pass. Run one dummy batch to
    # materialize them before counting parameters or handing them to Adam.
    model.eval()  # BatchNorm in the head rejects a 1-sample batch in train mode
    with torch.no_grad():
        model(
            x_dict,
            edge_dict,
            torch.from_numpy(cell_index[train_idx[:2]]).to(device),
            torch.from_numpy(drug_index[train_idx[:2]]).to(device),
        )
    n_params = sum(p.numel() for p in model.parameters())

    t1 = time.perf_counter()
    model, best_epoch = train(
        model, x_dict, edge_dict, loader,
        cell_index[val_idx], drug_index[val_idx], y[val_idx], device,
    )
    t_fit = time.perf_counter() - t1

    val = evaluate(y[val_idx], predict(model, x_dict, edge_dict, cell_index[val_idx], drug_index[val_idx], device), threshold)
    test = evaluate(y[test_idx], predict(model, x_dict, edge_dict, cell_index[test_idx], drug_index[test_idx], device), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])

    print_metric_block(f"GNN | {label}", val, test, floor)
    print(f"n_pairs={len(y)}  params={n_params:,}  best_epoch={best_epoch}  "
          f"prep={t_prep:.1f}s  fit={t_fit:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {
        "model": f"GNN-{variant.upper()}",
        "omics": "GE+Mut_CNV+Proteomics",
        "drug_rep": "fingerprint",
        "mutation_edges": use_mutation_edges,
        "n_pairs": len(y),
        "n_edge_types": len(edge_dict),
        "mean_only_rmse": floor,
        **{f"val_{k}": v for k, v in val.items()},
        **{f"test_{k}": v for k, v in test.items()},
        "params": n_params,
        "best_epoch": best_epoch,
        "fit_seconds": t_fit,
    }


def main() -> None:
    t_start = time.perf_counter()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    threshold = compute_shared_threshold()

    data, cell_index, drug_index, y, groups = load_graph_and_pairs(threshold)
    has_mutation_edges = MUTATION_EDGE in data.edge_types
    if not has_mutation_edges:
        print(f"\n[warn] {MUTATION_EDGE} not present in the graph — run "
              f"src/data/04_link_cell_lines.py to add it. Running only the "
              f"'without' arm of the ablation.")

    results = []
    for variant in ("gcn", "gat"):
        for use_mutation in ([False, True] if has_mutation_edges else [False]):
            results.append(
                run_one(data, cell_index, drug_index, y, groups, variant, use_mutation, threshold, device)
            )

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"GNN ABLATION COMPLETE — {len(df)} runs in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    print(df[["model", "mutation_edges", "n_pairs", "test_rmse", "test_pcc", "test_r2", "test_auc"]]
          .sort_values("test_rmse").to_string(index=False))
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
