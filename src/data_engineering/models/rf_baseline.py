"""
02_rf_baseline.py — Random Forest baseline for continuous IC50 prediction.
"""

from __future__ import annotations

import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupShuffleSplit

# --- config ---------------------------------------------------------------
DATA_DIR = Path("data/processed")
ALIGNED_DIR = Path("data/processed/aligned")

FUSED_NPY = DATA_DIR / "fused_early.npy"          
FUSED_IDS = DATA_DIR / "fused_cell_lines.csv"
TARGET_CSV = ALIGNED_DIR / "gdsc2_response_master.csv" 

# FIX: Point to the Sanger ID column to match the omics data index
COL_CELL_LINE = "sanger_model_id"
COL_DRUG = "drug_id"
COL_TARGET = "ln_ic50"


TEST_FRAC = 0.15
VAL_FRAC = 0.15
RANDOM_STATE = 42

N_ESTIMATORS = 100
MIN_SAMPLES_LEAF = 5
MAX_FEATURES = "sqrt"
N_JOBS = -1  
DTYPE = np.float32
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


def evaluate(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[float, float]:
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    pcc = float(np.corrcoef(y_true, y_pred)[0, 1]) if np.std(y_pred) > 0 else float("nan")
    return rmse, pcc


def main() -> None:
    t0 = time.perf_counter()

    X_cell, cell_ids = load_features()
    y_df = load_targets(cell_ids)
    X, y, groups, _ = build_pair_matrix(X_cell, cell_ids, y_df)
    del X_cell

    train_idx, val_idx, test_idx = grouped_split(groups)
    t_prep = time.perf_counter() - t0

    model = RandomForestRegressor(
        n_estimators=N_ESTIMATORS,
        min_samples_leaf=MIN_SAMPLES_LEAF,
        max_features=MAX_FEATURES,
        n_jobs=N_JOBS,
        random_state=RANDOM_STATE,
    )

    t1 = time.perf_counter()
    model.fit(X[train_idx], y[train_idx])
    t_fit = time.perf_counter() - t1

    t2 = time.perf_counter()
    val_rmse, val_pcc = evaluate(y[val_idx], model.predict(X[val_idx]))
    test_rmse, test_pcc = evaluate(y[test_idx], model.predict(X[test_idx]))
    t_pred = time.perf_counter() - t2

    base_rmse = float(np.sqrt(np.mean((y[test_idx] - y[train_idx].mean()) ** 2)))

    print("\n" + "=" * 58)
    print("RANDOM FOREST BASELINE — IC50 REGRESSION")
    print("=" * 58)
    print(f"{'':<12}{'RMSE':>10}{'PCC':>10}")
    print(f"{'Validation':<12}{val_rmse:>10.4f}{val_pcc:>10.4f}")
    print(f"{'Test':<12}{test_rmse:>10.4f}{test_pcc:>10.4f}")
    print(f"{'Mean-only':<12}{base_rmse:>10.4f}{'--':>10} ")
    print("-" * 58)
    print(f"n_estimators={N_ESTIMATORS}  min_samples_leaf={MIN_SAMPLES_LEAF}  "
          f"max_features={MAX_FEATURES}")
    print(f"Prep {t_prep:.1f}s | Fit {t_fit:.1f}s | Predict {t_pred:.1f}s | "
          f"Total {time.perf_counter() - t0:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")
    print("=" * 58)


if __name__ == "__main__":
    main()