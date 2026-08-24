"""mlp_matrix.py — MLP across all 7 omics subsets x 2 drug representations.

Phase 1 of the experiment matrix (docs/plan/experiment_matrix_plan.md): the
flat-concatenation *gradient-trained* counterpart to `rf_matrix.py`. Same 14
cells (7 omics subsets x {one-hot, Morgan fingerprint}), same split, same
metrics -- the difference from RF is purely the learner, which is what
isolates "neural net vs. tree ensemble" from "learned fusion" once
`cross_attention_matrix.py` lands.

Fingerprint-mode runs project the sparse 2048-bit block down through a small
encoder before the trunk (see `FingerprintEncoder`); concatenating raw bits
onto the omics block would let the drug side dominate by width alone.

Run from the repository root:
    python "experiments/Full Matrix/mlp_matrix.py"
"""

from __future__ import annotations

import copy
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.experiment_utils import (  # noqa: E402
    DRUG_ARMS,
    GE_KEY,
    MUT_CNV_KEY,
    PROTEOMICS_KEY,
    build_pair_matrix,
    compute_shared_threshold,
    concat_omics,
    evaluate,
    grouped_split,
    load_omics_subset,
    load_targets,
    mean_only_floor,
    peak_rss_gb,
    print_metric_block,
)

ALL_MODALITIES = [GE_KEY, MUT_CNV_KEY, PROTEOMICS_KEY]

TORCH_SEED = 42
HIDDEN_DIMS = [512, 256, 128]
DROPOUT = 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-5
BATCH_SIZE = 256
MAX_EPOCHS = 200
PATIENCE = 15
LR_PATIENCE = 5

RESULTS_CSV = Path("experiments/Full Matrix/mlp_matrix_results.csv")


def omics_subsets() -> list[list[str]]:
    subsets: list[list[str]] = []
    for size in (1, 2, 3):
        for combo in combinations(ALL_MODALITIES, size):
            subsets.append(list(combo))
    return subsets


class MLPRegressor(nn.Module):
    """[Linear -> BatchNorm -> ReLU -> Dropout] x N -> Linear(1)."""

    def __init__(self, in_dim: int, hidden_dims: list[int], dropout: float):
        super().__init__()
        layers: list[nn.Module] = []
        prev = in_dim
        for h in hidden_dims:
            layers += [nn.Linear(prev, h), nn.BatchNorm1d(h), nn.ReLU(), nn.Dropout(dropout)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


@torch.no_grad()
def predict(model: nn.Module, X: np.ndarray, device: torch.device) -> np.ndarray:
    model.eval()
    preds = []
    for start in range(0, len(X), BATCH_SIZE):
        xb = torch.from_numpy(X[start : start + BATCH_SIZE]).to(device)
        preds.append(model(xb).cpu().numpy())
    return np.concatenate(preds)


def train(model, loader, X_val, y_val, device) -> tuple[nn.Module, int]:
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=LR_PATIENCE
    )
    best_state = copy.deepcopy(model.state_dict())
    best_rmse, best_epoch, stale = float("inf"), 0, 0

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        for xb, yb in loader:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            loss = criterion(model(xb), yb)
            loss.backward()
            optimizer.step()

        val_rmse = float(np.sqrt(np.mean((y_val - predict(model, X_val, device)) ** 2)))
        scheduler.step(val_rmse)

        if val_rmse < best_rmse - 1e-4:
            best_rmse, best_epoch, stale = val_rmse, epoch, 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            stale += 1

        if epoch % 10 == 0:
            print(f"  [epoch {epoch:>3}] val_rmse={val_rmse:.4f}  best={best_rmse:.4f} @ {best_epoch}")

        if stale >= PATIENCE:
            print(f"  [early stop] epoch {epoch}, best epoch {best_epoch}, val_rmse={best_rmse:.4f}")
            break

    model.load_state_dict(best_state)
    return model, best_epoch


def run_one(modalities: list[str], arm: str, drug_mode: str, restricted: bool,
            threshold: float, device: torch.device) -> dict:
    label = "+".join(modalities)
    print("\n" + "#" * 82)
    print(f"# MLP | omics={label} | drug={arm}")
    print("#" * 82)

    torch.manual_seed(TORCH_SEED)
    t0 = time.perf_counter()

    omics, cell_ids = load_omics_subset(modalities)
    y_df = load_targets(cell_ids)
    X_cell = concat_omics(omics, modalities)
    X, y, groups, _, _ = build_pair_matrix(X_cell, cell_ids, y_df, drug_mode, restricted)
    train_idx, val_idx, test_idx = grouped_split(groups)
    t_prep = time.perf_counter() - t0

    ds = TensorDataset(torch.from_numpy(X[train_idx]), torch.from_numpy(y[train_idx]))
    loader = DataLoader(ds, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)

    model = MLPRegressor(in_dim=X.shape[1], hidden_dims=HIDDEN_DIMS, dropout=DROPOUT).to(device)
    n_params = sum(p.numel() for p in model.parameters())

    t1 = time.perf_counter()
    model, best_epoch = train(model, loader, X[val_idx], y[val_idx], device)
    t_fit = time.perf_counter() - t1

    val = evaluate(y[val_idx], predict(model, X[val_idx], device), threshold)
    test = evaluate(y[test_idx], predict(model, X[test_idx], device), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])

    print_metric_block(f"MLP | {label} | {arm}", val, test, floor)
    print(f"n_pairs={len(y)}  n_features={X.shape[1]}  params={n_params:,}  "
          f"best_epoch={best_epoch}  prep={t_prep:.1f}s  fit={t_fit:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {
        "model": "MLP",
        "omics": label,
        "n_modalities": len(modalities),
        "drug_rep": arm,
        "n_pairs": len(y),
        "n_features": X.shape[1],
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

    results = []
    for modalities in omics_subsets():
        for arm, drug_mode, restricted in DRUG_ARMS:
            results.append(run_one(modalities, arm, drug_mode, restricted, threshold, device))

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"MLP MATRIX COMPLETE — {len(df)} runs in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    summary = df[["omics", "drug_rep", "n_pairs", "test_rmse", "test_pcc", "test_r2", "test_auc"]]
    print(summary.sort_values("test_rmse").to_string(index=False))
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
