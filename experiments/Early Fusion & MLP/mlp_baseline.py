"""
mlp_baseline.py — Simple MLP baseline for continuous IC50 prediction.

Deep-learning bridge between the flat-ML baselines (Random Forest) and the
Graph-ML models: same early-fusion (concatenated) input, same train/val/test
split, but a feed-forward network trained with gradient descent instead of
an ensemble of trees.
"""

from __future__ import annotations

import copy
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import GroupShuffleSplit

# --- config ---------------------------------------------------------------
DATA_DIR = Path("data/processed")
ALIGNED_DIR = Path("data/processed/aligned")

FUSED_NPY = DATA_DIR / "fused_early.npy"
FUSED_IDS = DATA_DIR / "fused_cell_lines.csv"
TARGET_CSV = ALIGNED_DIR / "gdsc2_response_master.csv"

COL_CELL_LINE = "sanger_model_id"
COL_DRUG = "drug_id"
COL_TARGET = "ln_ic50"

TEST_FRAC = 0.15
VAL_FRAC = 0.15
RANDOM_STATE = 42  # split RNG; pinned so every experiment-matrix run sits on an identical split

DTYPE = np.float32

# --- model / training hyperparameters --------------------------------------
TORCH_SEED = 42  # seeds weight init + batch shuffling; kept separate from the data-split RNG
HIDDEN_DIMS = [512, 256, 128]
DROPOUT = 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-5
BATCH_SIZE = 256
MAX_EPOCHS = 200
PATIENCE = 15  # early-stopping patience, in epochs without val RMSE improvement
LR_PATIENCE = 5  # ReduceLROnPlateau patience
# --------------------------------------------------------------------------


def peak_rss_gb() -> float:
    """Peak resident set size. Returns 0.0 on Windows to prevent crashes."""
    if platform.system() == "Windows":
        return 0.0

    import resource
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw / 1024**3 if platform.system() == "Darwin" else raw / 1024**2


def load_features() -> tuple[np.ndarray, pd.Index]:
    """Load the fused matrix from cache."""
    if not (FUSED_NPY.exists() and FUSED_IDS.exists()):
        raise FileNotFoundError(
            f"Missing {FUSED_NPY.name} or {FUSED_IDS.name}. "
            "Please run load_and_fuse.py first to generate these files."
        )

    X = np.load(FUSED_NPY, mmap_mode=None).astype(DTYPE, copy=False)
    ids = pd.Index(pd.read_csv(FUSED_IDS, header=None).iloc[:, 0].astype(str))
    print(f"[features] loaded cache {FUSED_NPY.name} {X.shape}")

    if len(ids) != X.shape[0]:
        raise ValueError(f"ID count {len(ids)} != feature rows {X.shape[0]}")
    return X, ids


def load_targets(valid_ids: pd.Index) -> pd.DataFrame:
    if not TARGET_CSV.exists():
        raise FileNotFoundError(
            f"Missing target file at {TARGET_CSV}. "
            "Ensure the DE pipeline (ingest_and_align.py) has been run."
        )

    y = pd.read_csv(TARGET_CSV, usecols=[COL_CELL_LINE, COL_DRUG, COL_TARGET])
    n_raw = len(y)

    y[COL_CELL_LINE] = y[COL_CELL_LINE].astype(str)
    y[COL_DRUG] = y[COL_DRUG].astype(str)
    y[COL_TARGET] = pd.to_numeric(y[COL_TARGET], errors="coerce")

    y = y[np.isfinite(y[COL_TARGET])]
    n_finite = len(y)

    y = y[y[COL_CELL_LINE].isin(set(valid_ids))]
    n_matched = len(y)

    dup = y.duplicated([COL_CELL_LINE, COL_DRUG]).sum()
    if dup:
        print(f"[targets] {dup} duplicate (cell_line, drug) rows -> averaging")
        y = y.groupby([COL_CELL_LINE, COL_DRUG], as_index=False)[COL_TARGET].mean()

    print(f"[targets] {n_raw} rows -> {n_finite} finite -> {n_matched} with omics "
          f"-> {len(y)} unique pairs")
    if y.empty:
        raise ValueError("No target rows survived filtering; check ID formatting.")
    return y.reset_index(drop=True)


def build_pair_matrix(
    X_cell: np.ndarray, cell_ids: pd.Index, y: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Gather cell-line features per pair and append a one-hot drug block."""
    row_of = pd.Series(np.arange(len(cell_ids)), index=cell_ids)
    rows = row_of.loc[y[COL_CELL_LINE]].to_numpy()

    drug_codes, drug_levels = pd.factorize(y[COL_DRUG], sort=True)
    n_pairs, n_omics, n_drugs = len(y), X_cell.shape[1], len(drug_levels)

    X = np.empty((n_pairs, n_omics + n_drugs), dtype=DTYPE)
    X[:, :n_omics] = X_cell[rows]                       # fancy-index gather
    X[:, n_omics:] = 0.0
    X[np.arange(n_pairs), n_omics + drug_codes] = 1.0   # one-hot

    target = y[COL_TARGET].to_numpy(dtype=DTYPE)
    groups = y[COL_CELL_LINE].to_numpy()
    names = [f"omics_{i}" for i in range(n_omics)] + [f"drug={d}" for d in drug_levels]

    print(f"[design]  {n_pairs} pairs x {X.shape[1]} features "
          f"({n_omics} omics + {n_drugs} drug one-hot) = {X.nbytes / 1024**2:.1f} MB")
    return X, target, groups, names


def grouped_split(groups: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """70/15/15 by cell line."""
    holdout = VAL_FRAC + TEST_FRAC
    gss1 = GroupShuffleSplit(n_splits=1, test_size=holdout, random_state=RANDOM_STATE)
    train_idx, rest_idx = next(gss1.split(np.zeros(len(groups)), groups=groups))

    gss2 = GroupShuffleSplit(
        n_splits=1, test_size=TEST_FRAC / holdout, random_state=RANDOM_STATE
    )
    rel_val, rel_test = next(
        gss2.split(np.zeros(len(rest_idx)), groups=groups[rest_idx])
    )
    val_idx, test_idx = rest_idx[rel_val], rest_idx[rel_test]

    for name, idx in [("train", train_idx), ("val", val_idx), ("test", test_idx)]:
        print(f"[split]   {name:<5} {len(idx):>7} pairs  "
              f"{len(np.unique(groups[idx])):>5} cell lines  "
              f"({len(idx) / len(groups):.1%})")

    overlap = set(groups[train_idx]) & (set(groups[val_idx]) | set(groups[test_idx]))
    assert not overlap, f"cell line leaked across splits: {sorted(overlap)[:5]}"
    return train_idx, val_idx, test_idx


def evaluate(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[float, float, float]:
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    pcc = float(np.corrcoef(y_true, y_pred)[0, 1]) if np.std(y_pred) > 0 else float("nan")
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    return rmse, pcc, r2


class MLPRegressor(nn.Module):
    """Feed-forward network: [Linear -> BatchNorm -> ReLU -> Dropout] x N -> Linear(1)."""

    def __init__(self, in_dim: int, hidden_dims: list[int], dropout: float):
        super().__init__()
        layers: list[nn.Module] = []
        prev = in_dim
        for h in hidden_dims:
            layers += [
                nn.Linear(prev, h),
                nn.BatchNorm1d(h),
                nn.ReLU(),
                nn.Dropout(dropout),
            ]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


def make_loader(X: np.ndarray, y: np.ndarray, batch_size: int, shuffle: bool) -> DataLoader:
    ds = TensorDataset(torch.from_numpy(X), torch.from_numpy(y))
    # drop_last avoids a size-1 trailing batch, which BatchNorm1d rejects in train mode
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle, drop_last=shuffle)


@torch.no_grad()
def predict(model: nn.Module, X: np.ndarray, device: torch.device, batch_size: int) -> np.ndarray:
    model.eval()
    preds = []
    for start in range(0, len(X), batch_size):
        xb = torch.from_numpy(X[start : start + batch_size]).to(device)
        preds.append(model(xb).cpu().numpy())
    return np.concatenate(preds)


def train(
    model: nn.Module,
    train_loader: DataLoader,
    X_val: np.ndarray,
    y_val: np.ndarray,
    device: torch.device,
) -> tuple[nn.Module, int]:
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=LR_PATIENCE
    )

    best_state = copy.deepcopy(model.state_dict())
    best_val_rmse = float("inf")
    best_epoch = 0
    epochs_since_improve = 0

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        train_loss = 0.0
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            pred = model(xb)
            loss = criterion(pred, yb)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * len(xb)
        train_loss /= len(train_loader.dataset)

        val_pred = predict(model, X_val, device, BATCH_SIZE)
        val_rmse, val_pcc, val_r2 = evaluate(y_val, val_pred)
        scheduler.step(val_rmse)

        improved = val_rmse < best_val_rmse - 1e-4
        if improved:
            best_val_rmse = val_rmse
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            epochs_since_improve = 0
        else:
            epochs_since_improve += 1

        if epoch == 1 or epoch % 5 == 0 or improved:
            marker = " *" if improved else ""
            print(f"[epoch {epoch:>3}] train_loss={train_loss:.4f}  "
                  f"val_rmse={val_rmse:.4f}  val_pcc={val_pcc:.4f}  val_r2={val_r2:.4f}{marker}")

        if epochs_since_improve >= PATIENCE:
            print(f"[early stop] no val improvement for {PATIENCE} epochs "
                  f"(best epoch {best_epoch}, val_rmse={best_val_rmse:.4f})")
            break

    model.load_state_dict(best_state)
    return model, best_epoch


def main() -> None:
    t0 = time.perf_counter()

    torch.manual_seed(TORCH_SEED)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    X_cell, cell_ids = load_features()
    y_df = load_targets(cell_ids)
    X, y, groups, _ = build_pair_matrix(X_cell, cell_ids, y_df)
    del X_cell

    train_idx, val_idx, test_idx = grouped_split(groups)
    t_prep = time.perf_counter() - t0

    train_loader = make_loader(X[train_idx], y[train_idx], BATCH_SIZE, shuffle=True)

    model = MLPRegressor(in_dim=X.shape[1], hidden_dims=HIDDEN_DIMS, dropout=DROPOUT).to(device)
    n_params = sum(p.numel() for p in model.parameters())

    t1 = time.perf_counter()
    model, best_epoch = train(model, train_loader, X[val_idx], y[val_idx], device)
    t_fit = time.perf_counter() - t1

    t2 = time.perf_counter()
    val_rmse, val_pcc, val_r2 = evaluate(y[val_idx], predict(model, X[val_idx], device, BATCH_SIZE))
    test_rmse, test_pcc, test_r2 = evaluate(y[test_idx], predict(model, X[test_idx], device, BATCH_SIZE))
    t_pred = time.perf_counter() - t2

    base_rmse = float(np.sqrt(np.mean((y[test_idx] - y[train_idx].mean()) ** 2)))

    print("\n" + "=" * 58)
    print("MLP BASELINE — IC50 REGRESSION")
    print("=" * 58)
    print(f"{'':<12}{'RMSE':>10}{'PCC':>10}{'R^2':>10}")
    print(f"{'Validation':<12}{val_rmse:>10.4f}{val_pcc:>10.4f}{val_r2:>10.4f}")
    print(f"{'Test':<12}{test_rmse:>10.4f}{test_pcc:>10.4f}{test_r2:>10.4f}")
    print(f"{'Mean-only':<12}{base_rmse:>10.4f}{'--':>10}{'--':>10} ")
    print("-" * 58)
    print(f"hidden_dims={HIDDEN_DIMS}  dropout={DROPOUT}  lr={LR}  "
          f"weight_decay={WEIGHT_DECAY}  batch_size={BATCH_SIZE}")
    print(f"params={n_params:,}  best_epoch={best_epoch}  device={device}")
    print(f"Prep {t_prep:.1f}s | Fit {t_fit:.1f}s | Predict {t_pred:.1f}s | "
          f"Total {time.perf_counter() - t0:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")
    print("=" * 58)


if __name__ == "__main__":
    main()
