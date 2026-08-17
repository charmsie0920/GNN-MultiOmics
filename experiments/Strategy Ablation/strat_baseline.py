"""
strat_baseline.py — Late Fusion Ensemble Random Forest for continuous IC50
prediction. Trains three specialized RF models in isolation (Genomics,
Transcriptomics, Proteomics), then combines their predictions at inference
time via a validation-RMSE-weighted average. Tests whether three specialized
models beat the single early-fused Control model (rf_baseline.py, RMSE 2.19).

Run from the repository root:
    python "experiments/Strategy Ablation/strat_baseline.py"
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

TARGET_CSV = ALIGNED_DIR / "gdsc2_response_master.csv"

# Late fusion ensemble: train one RF per modality in isolation, then average.
ALL_OMICS = {
    "Genomics": DATA_DIR / "genomics_pca.csv",
    "Transcriptomics": DATA_DIR / "transcriptomics_pca.csv",
    "Proteomics": DATA_DIR / "proteomics_pca.csv",
}

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
    """Median ln_ic50 over all target rows, used as a fixed cutoff for every model."""
    y_all = pd.read_csv(TARGET_CSV, usecols=[COL_TARGET])
    y_all[COL_TARGET] = pd.to_numeric(y_all[COL_TARGET], errors="coerce")
    threshold = float(y_all[COL_TARGET].median(skipna=True))
    print(f"[threshold] median ln_ic50 across all targets = {threshold:.4f}")
    return threshold


def print_block(title: str, val_metrics: tuple, test_metrics: tuple,
                 mean_only_rmse: float | None, extra_lines: list[str]) -> None:
    val_rmse, val_pcc, val_auc, val_f1 = val_metrics
    test_rmse, test_pcc, test_auc, test_f1 = test_metrics

    print("\n" + "=" * 68)
    print(title)
    print("=" * 68)
    print(f"{'':<12}{'RMSE':>10}{'PCC':>10}{'AUC':>10}{'F1':>10}")
    print(f"{'Validation':<12}{val_rmse:>10.4f}{val_pcc:>10.4f}{val_auc:>10.4f}{val_f1:>10.4f}")
    print(f"{'Test':<12}{test_rmse:>10.4f}{test_pcc:>10.4f}{test_auc:>10.4f}{test_f1:>10.4f}")
    if mean_only_rmse is not None:
        print(f"{'Mean-only':<12}{mean_only_rmse:>10.4f}{'--':>10}{'--':>10}{'--':>10}")
    print("-" * 68)
    for line in extra_lines:
        print(line)
    print("=" * 68)


def train_omics_model(
    omics_name: str, X_cell: np.ndarray, cell_ids: pd.Index, y_df: pd.DataFrame,
    train_idx: np.ndarray, val_idx: np.ndarray, test_idx: np.ndarray, threshold: float,
) -> dict:
    """Fit one specialized RF on a single omics block, using the shared split."""
    t0 = time.perf_counter()

    X, y, _, _ = build_pair_matrix(X_cell, cell_ids, y_df)

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
    val_metrics = evaluate(y[val_idx], val_pred, threshold)
    test_metrics = evaluate(y[test_idx], test_pred, threshold)
    t_pred = time.perf_counter() - t2

    mean_only_rmse = float(np.sqrt(np.mean((y[test_idx] - y[train_idx].mean()) ** 2)))

    print_block(
        f"RANDOM FOREST — {omics_name.upper()} (ISOLATED) — IC50 REGRESSION",
        val_metrics, test_metrics, mean_only_rmse,
        [
            f"n_estimators={N_ESTIMATORS}  min_samples_leaf={MIN_SAMPLES_LEAF}  "
            f"max_features={MAX_FEATURES}",
            f"Fit {t_fit:.1f}s | Predict {t_pred:.1f}s | Total {time.perf_counter() - t0:.1f}s "
            f"| Peak RSS: {peak_rss_gb():.2f} GB",
        ],
    )

    return {
        "omics": omics_name,
        "val_rmse": val_metrics[0], "val_pcc": val_metrics[1],
        "val_auc": val_metrics[2], "val_f1": val_metrics[3],
        "test_rmse": test_metrics[0], "test_pcc": test_metrics[1],
        "test_auc": test_metrics[2], "test_f1": test_metrics[3],
        "mean_only_rmse": mean_only_rmse,
        "val_pred": val_pred, "test_pred": test_pred,
    }


def run_ensemble(model_results: list[dict], y: np.ndarray,
                  val_idx: np.ndarray, test_idx: np.ndarray, threshold: float) -> dict:
    """Combine per-model predictions via a validation-RMSE-weighted average."""
    inv_rmse = np.array([1.0 / r["val_rmse"] for r in model_results])
    weights = inv_rmse / inv_rmse.sum()

    val_pred = sum(w * r["val_pred"] for w, r in zip(weights, model_results))
    test_pred = sum(w * r["test_pred"] for w, r in zip(weights, model_results))

    val_metrics = evaluate(y[val_idx], val_pred, threshold)
    test_metrics = evaluate(y[test_idx], test_pred, threshold)

    weight_line = "Weights: " + "  ".join(
        f"{r['omics']}={w:.3f}" for r, w in zip(model_results, weights)
    )

    print_block(
        "LATE FUSION ENSEMBLE (VALIDATION-RMSE-WEIGHTED AVG) — IC50 REGRESSION",
        val_metrics, test_metrics, None, [weight_line],
    )

    return {
        "omics": "Ensemble (weighted)",
        "val_rmse": val_metrics[0], "val_pcc": val_metrics[1],
        "val_auc": val_metrics[2], "val_f1": val_metrics[3],
        "test_rmse": test_metrics[0], "test_pcc": test_metrics[1],
        "test_auc": test_metrics[2], "test_f1": test_metrics[3],
        "mean_only_rmse": float("nan"),
    }


def print_summary(results: list[dict]) -> None:
    print("\n" + "=" * 100)
    print("SUMMARY — LATE FUSION ENSEMBLE ABLATION")
    print("=" * 100)
    header = (
        f"{'Model':<20}{'Val RMSE':>10}{'Val PCC':>10}{'Val AUC':>10}{'Val F1':>9}"
        f"{'Test RMSE':>11}{'Test PCC':>10}{'Test AUC':>10}{'Test F1':>9}{'Mean RMSE':>11}"
    )
    print(header)
    print("-" * 100)
    for r in results:
        mean_rmse_str = f"{r['mean_only_rmse']:>11.4f}" if not np.isnan(r["mean_only_rmse"]) else f"{'--':>11}"
        print(
            f"{r['omics']:<20}{r['val_rmse']:>10.4f}{r['val_pcc']:>10.4f}"
            f"{r['val_auc']:>10.4f}{r['val_f1']:>9.4f}"
            f"{r['test_rmse']:>11.4f}{r['test_pcc']:>10.4f}"
            f"{r['test_auc']:>10.4f}{r['test_f1']:>9.4f}{mean_rmse_str}"
        )
    print("=" * 100)


def main() -> None:
    threshold = compute_shared_threshold()

    # Load all three omics blocks and find the cell lines common to all of them,
    # so every model trains/evaluates on the exact same (cell_line, drug) pairs.
    omics_blocks = {name: load_single_omics_features(path) for name, path in ALL_OMICS.items()}
    shared_ids = None
    for _, ids in omics_blocks.values():
        shared_ids = ids if shared_ids is None else shared_ids.intersection(ids)
    print(f"[shared] {len(shared_ids)} cell lines common to all {len(ALL_OMICS)} omics")

    y_df = load_targets(shared_ids)

    # Shared split: build one omics' pair matrix just to get the shared `groups`
    # array (identical across omics since y_df/order is fixed), split once.
    first_name = next(iter(omics_blocks))
    X_first, ids_first = omics_blocks[first_name]
    _, y_shared, groups, _ = build_pair_matrix(X_first, ids_first, y_df)
    train_idx, val_idx, test_idx = grouped_split(groups)

    model_results = []
    for omics_name, (X_cell, cell_ids) in omics_blocks.items():
        model_results.append(
            train_omics_model(
                omics_name, X_cell, cell_ids, y_df,
                train_idx, val_idx, test_idx, threshold,
            )
        )

    ensemble_result = run_ensemble(model_results, y_shared, val_idx, test_idx, threshold)

    for r in model_results:
        del r["val_pred"], r["test_pred"]

    print_summary(model_results + [ensemble_result])


if __name__ == "__main__":
    main()
