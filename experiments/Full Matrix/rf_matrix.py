"""rf_matrix.py — Random Forest across all 7 omics subsets x 2 drug representations.

Phase 1 of the experiment matrix (docs/plan/experiment_matrix_plan.md): the
flat-concatenation tree baseline, swept across every omics subset (3 single,
3 dual, 1 tri) and both drug representations (one-hot identity vs. 2048-bit
Morgan fingerprint). 14 runs total.

Every run uses the identical GroupShuffleSplit-by-cell-line 70/15/15 split
(`random_state=42`, via src/data/experiment_utils.py) so results are directly
comparable across cells. Note that fingerprint-mode runs cover a smaller pair
population than one-hot runs (drugs whose SMILES never resolved are dropped),
so fingerprint-vs-one-hot comparisons are across slightly different row
counts -- reported per-run so the docs can state it explicitly.

Run from the repository root:
    python "experiments/Full Matrix/rf_matrix.py"
"""

from __future__ import annotations

import sys
import time
from itertools import combinations
from pathlib import Path

import pandas as pd
from sklearn.ensemble import RandomForestRegressor

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.experiment_utils import (  # noqa: E402
    DRUG_ARMS,
    GE_KEY,
    MUT_CNV_KEY,
    PROTEOMICS_KEY,
    RANDOM_STATE,
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

N_ESTIMATORS = 100
MIN_SAMPLES_LEAF = 5
MAX_FEATURES = "sqrt"
N_JOBS = -1

RESULTS_CSV = Path("experiments/Full Matrix/rf_matrix_results.csv")


def omics_subsets() -> list[list[str]]:
    """All 7 non-empty subsets of the 3 modalities, ordered by size then name."""
    subsets: list[list[str]] = []
    for size in (1, 2, 3):
        for combo in combinations(ALL_MODALITIES, size):
            subsets.append(list(combo))
    return subsets


def run_one(modalities: list[str], arm: str, drug_mode: str, restricted: bool, threshold: float) -> dict:
    label = "+".join(modalities)
    print("\n" + "#" * 82)
    print(f"# RF | omics={label} | drug={arm}")
    print("#" * 82)

    t0 = time.perf_counter()
    omics, cell_ids = load_omics_subset(modalities)
    y_df = load_targets(cell_ids)
    X_cell = concat_omics(omics, modalities)
    X, y, groups, _, _ = build_pair_matrix(X_cell, cell_ids, y_df, drug_mode, restricted)
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

    val = evaluate(y[val_idx], model.predict(X[val_idx]), threshold)
    test = evaluate(y[test_idx], model.predict(X[test_idx]), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])

    print_metric_block(f"RF | {label} | {arm}", val, test, floor)
    print(f"n_pairs={len(y)}  n_features={X.shape[1]}  prep={t_prep:.1f}s  fit={t_fit:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {
        "model": "RF",
        "omics": label,
        "n_modalities": len(modalities),
        "drug_rep": arm,
        "n_pairs": len(y),
        "n_features": X.shape[1],
        "mean_only_rmse": floor,
        **{f"val_{k}": v for k, v in val.items()},
        **{f"test_{k}": v for k, v in test.items()},
        "params": model.n_estimators,
        "fit_seconds": t_fit,
    }


def main() -> None:
    t_start = time.perf_counter()
    threshold = compute_shared_threshold()

    results = []
    for modalities in omics_subsets():
        for arm, drug_mode, restricted in DRUG_ARMS:
            results.append(run_one(modalities, arm, drug_mode, restricted, threshold))

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"RF MATRIX COMPLETE — {len(df)} runs in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    summary = df[["omics", "drug_rep", "n_pairs", "test_rmse", "test_pcc", "test_r2", "test_auc"]]
    print(summary.sort_values("test_rmse").to_string(index=False))
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
