"""hetero_ic50_gnn_matrix.py — the tuned HeteroIC50GNN, under the matrix protocol.

`src/models/test/hetero_gnn.py::HeteroIC50GNN` was repaired and tuned in
docs/hetero_gnn_test_bugfixes.md, reaching RMSE 1.3035 — better than the
tracked `HeteroGNN` GCN run (E12, 1.3513). But that number came from
`src/models/test/05_train.py`, which builds its own split in
`src/models/test/dataset.py` and reports only MSE/RMSE/Pearson r. It therefore
could not go into docs/results.md, which is auto-generated from per-family CSVs
and needs the full metric set on the shared split.

This script closes that gap: same model, but run through the matrix's own
`load_graph_and_pairs` -> `grouped_split` -> `evaluate` path, exactly as
`gnn_baseline.py` does for `HeteroGNN`. That makes the emitted row directly
comparable to every other row in the table rather than merely similar-looking
— the split function, the binarization threshold, and the metric definitions
are all shared code, not re-derived.

Training uses the configuration that won the tuning sweep: AdamW lr=1e-3
weight_decay=1e-4, ReduceLROnPlateau(factor=0.5, patience=5), early stopping
patience=15 on val RMSE.

Run from the repository root, after 03_graph_construction.py and
04_link_cell_lines.py:
    python "experiments/07_gnn_ablation/hetero_ic50_gnn_matrix.py"
"""

from __future__ import annotations

import copy
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.data.experiment_utils import (  # noqa: E402
    compute_shared_threshold,
    evaluate,
    grouped_split,
    mean_only_floor,
    peak_rss_gb,
    print_metric_block,
)
from src.models.test.hetero_gnn import HeteroIC50GNN  # noqa: E402

RESULTS_CSV = Path("experiments/07_gnn_ablation/hetero_ic50_gnn_results.csv")

# The tuning sweep's winning configuration (docs/hetero_gnn_test_bugfixes.md).
LR = 1e-3
WEIGHT_DECAY = 1e-4
HIDDEN_DIM = 128
BATCH_SIZE = 1024
MAX_EPOCHS = 200
PATIENCE = 15
LR_PATIENCE = 5
TORCH_SEED = 42


def load_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# Reuse the tracked GNN experiment's graph/pair loader so this run sits on an
# identical population to E12/E24.
gnn_mod = load_module("gnn_baseline", "experiments/07_gnn_ablation/gnn_baseline.py")


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
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
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


def main() -> None:
    t_start = time.perf_counter()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    threshold = compute_shared_threshold()

    data, cell_index, drug_index, y, groups = gnn_mod.load_graph_and_pairs(threshold)

    print("\n" + "#" * 82)
    print("# MATRIX | HeteroIC50GNN (tuned) | fingerprint | +mutation edges")
    print("#" * 82)

    torch.manual_seed(TORCH_SEED)
    t0 = time.perf_counter()

    x_dict = {nt: data[nt].x.to(device, dtype=torch.float32) for nt in data.node_types}
    edge_dict = {
        et: data[et].edge_index.to(device, dtype=torch.long) for et in data.edge_types
    }
    # HeteroIC50GNN declares both reverse relations explicitly.
    src, dst = edge_dict[("drug", "targets", "protein")]
    edge_dict[("protein", "rev_targets", "drug")] = torch.stack([dst, src])
    src, dst = edge_dict[("cell_line", "has_mutation", "protein")]
    edge_dict[("protein", "rev_has_mutation", "cell_line")] = torch.stack([dst, src])
    print(f"[edges] using {len(edge_dict)} edge types: {list(edge_dict.keys())}")

    # The shared 70/15/15 cell-line-grouped split every matrix row uses.
    train_idx, val_idx, test_idx = grouped_split(groups)
    t_prep = time.perf_counter() - t0

    ds = TensorDataset(
        torch.from_numpy(cell_index[train_idx]),
        torch.from_numpy(drug_index[train_idx]),
        torch.from_numpy(y[train_idx]),
    )
    loader = DataLoader(ds, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)

    model = HeteroIC50GNN(
        num_proteins=data["protein"].x.shape[0], hidden_dim=HIDDEN_DIM
    ).to(device)
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

    print_metric_block("MATRIX | HeteroIC50GNN (tuned)", val, test, floor)
    print(f"n_pairs={len(y)}  params={n_params:,}  best_epoch={best_epoch}  "
          f"prep={t_prep:.1f}s  fit={t_fit:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    row = {
        "model": "HeteroIC50GNN",
        "omics": "GE+Mut_CNV+Proteomics",
        "drug_rep": "fingerprint",
        "mutation_edges": True,
        "n_pairs": len(y),
        "n_edge_types": len(edge_dict),
        "mean_only_rmse": floor,
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
    print(f"COMPLETE in {time.perf_counter() - t_start:.1f}s — test RMSE {test['rmse']:.4f}")
    print("=" * 100)
    print(f"Saved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
