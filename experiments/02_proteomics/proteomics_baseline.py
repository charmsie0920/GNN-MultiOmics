"""
proteomics_baseline.py — Random Forest baseline for continuous IC50 prediction,
using Proteomics as the sole omics modality (no early fusion with GE/Mut_CNV).
Reports RMSE, PCC, AUC, and F1 on validation and test.

Run from the repository root:
    python "experiments/02_proteomics/proteomics_baseline.py"
"""

from __future__ import annotations

import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import f1_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit

# --- config ---------------------------------------------------------------
DATA_DIR = Path("data/processed")
ALIGNED_DIR = Path("data/processed/aligned")

PROTEOMICS_CSV = DATA_DIR / "proteomics_pca.csv"
TARGET_CSV = ALIGNED_DIR / "gdsc2_response_master.csv"

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


def load_single_omics_features(path: Path) -> tuple[np.ndarray, pd.Index]:
    """Load one PCA-compressed omics CSV (cell lines in column 0) as float32."""
    if not path.exists():
        raise FileNotFoundError(f"Missing omics file at {path}.")

    block = pd.read_csv(path, index_col=0)
    ids = block.index
    X = block.to_numpy(dtype=DTYPE, copy=False)
    print(f"[features] loaded {path.name} {X.shape}")
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


def evaluate(
    y_true: np.ndarray, y_pred: np.ndarray, threshold: float
) -> tuple[float, float, float, float]:
    """RMSE/PCC on the continuous target, plus AUC/F1 on a threshold-binarized target."""
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    pcc = float(np.corrcoef(y_true, y_pred)[0, 1]) if np.std(y_pred) > 0 else float("nan")

    y_true_bin = (y_true >= threshold).astype(int)
    y_pred_bin = (y_pred >= threshold).astype(int)
    auc = (
        float(roc_auc_score(y_true_bin, y_pred))
        if len(np.unique(y_true_bin)) > 1
        else float("nan")
    )
    f1 = float(f1_score(y_true_bin, y_pred_bin, zero_division=0))
    return rmse, pcc, auc, f1


def compute_shared_threshold() -> float:
    """Median ln_ic50 over all target rows, used as the AUC/F1 binarization cutoff."""
    y_all = pd.read_csv(TARGET_CSV, usecols=[COL_TARGET])
    y_all[COL_TARGET] = pd.to_numeric(y_all[COL_TARGET], errors="coerce")
    threshold = float(y_all[COL_TARGET].median(skipna=True))
    print(f"[threshold] median ln_ic50 across all targets = {threshold:.4f}")
    return threshold


def main() -> None:
    t0 = time.perf_counter()

    threshold = compute_shared_threshold()

    X_cell, cell_ids = load_single_omics_features(PROTEOMICS_CSV)
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
    val_pred = model.predict(X[val_idx])
    test_pred = model.predict(X[test_idx])
    val_rmse, val_pcc, val_auc, val_f1 = evaluate(y[val_idx], val_pred, threshold)
    test_rmse, test_pcc, test_auc, test_f1 = evaluate(y[test_idx], test_pred, threshold)
    t_pred = time.perf_counter() - t2

    base_rmse = float(np.sqrt(np.mean((y[test_idx] - y[train_idx].mean()) ** 2)))

    print("\n" + "=" * 68)
    print("RANDOM FOREST BASELINE — PROTEOMICS — IC50 REGRESSION")
    print("=" * 68)
    print(f"{'':<12}{'RMSE':>10}{'PCC':>10}{'AUC':>10}{'F1':>10}")
    print(f"{'Validation':<12}{val_rmse:>10.4f}{val_pcc:>10.4f}{val_auc:>10.4f}{val_f1:>10.4f}")
    print(f"{'Test':<12}{test_rmse:>10.4f}{test_pcc:>10.4f}{test_auc:>10.4f}{test_f1:>10.4f}")
    print(f"{'Mean-only':<12}{base_rmse:>10.4f}{'--':>10}{'--':>10}{'--':>10}")
    print("-" * 68)
    print(f"n_estimators={N_ESTIMATORS}  min_samples_leaf={MIN_SAMPLES_LEAF}  "
          f"max_features={MAX_FEATURES}")
    print(f"Prep {t_prep:.1f}s | Fit {t_fit:.1f}s | Predict {t_pred:.1f}s | "
          f"Total {time.perf_counter() - t0:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")
    print("=" * 68)


if __name__ == "__main__":
    main()
