"""pair_gate_report.py — where the pair features help, if anywhere (Phase 3 gate).

Overall RMSE cannot answer the gate. A target is directly mutated in only
1.26% of pairs, so even a perfect use of that fact moves all-pairs RMSE by
about 0.005, far inside the ~0.03 seed noise. This script scores the saved
predictions of a gate run on the subsets where the features say something:

    all          every pair
    has_target   the drug has a known target (the pair module's coverage)
    no_target    it does not; the features are all-zero flags here
    direct_hit   a target protein is itself mutated in the cell line
    hop_1        nearest mutated protein is one PPI hop from a target
    hop_2plus    two or more hops

For each subset it reports RMSE per config and seed, the per-drug-mean floor
on the same rows, and the mean signed error (prediction - truth). A base model
whose signed error on `direct_hit` is near zero already knows the effect from
its omics input, and the features have nothing to add.

Each config is then compared with the baseline on the same seeds, with a
paired bootstrap over cell lines (2,000 resamples) for the difference.

Run from the repository root, after the gate run:
    python src/final_model/pair_gate_report.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.experiment_utils import (  # noqa: E402
    COL_CELL_LINE,
    COL_DRUG,
    COL_TARGET,
    DTYPE,
    RANDOM_STATE,
    grouped_split,
)
from src.data.pair_features import build_pair_features, load_pair_rows  # noqa: E402
from src.final_model.run_ablation import predictions_dir  # noqa: E402

GATE_CSV = Path("src/final_model/results/pair_gate_results.csv")
N_BOOTSTRAP = 2000
FOLDS = ("val", "test")


def rmse(error: np.ndarray) -> float:
    return float(np.sqrt(np.mean(error**2))) if len(error) else float("nan")


def bootstrap_delta(sq_err: np.ndarray, sq_err_base: np.ndarray, cells: np.ndarray,
                    rng: np.random.Generator) -> tuple[float, float, float]:
    """95% CI and P(delta >= 0) for RMSE(config) - RMSE(baseline), resampling cell lines.

    `sq_err` and `sq_err_base` are (n_seeds, n_rows) squared errors on the same
    rows. Each resample draws cell lines with replacement, takes every row of
    each drawn line, and averages the RMSE difference over seeds.
    """
    codes, levels = pd.factorize(cells)
    n_cells = len(levels)
    count = np.bincount(codes, minlength=n_cells).astype(np.float64)
    per_cell = lambda e: np.stack([np.bincount(codes, w, minlength=n_cells) for w in e])  # noqa: E731
    sse, sse_base = per_cell(sq_err), per_cell(sq_err_base)

    draws = rng.integers(0, n_cells, size=(N_BOOTSTRAP, n_cells))
    n_rows = count[draws].sum(axis=1)                        # (resamples,)
    with np.errstate(invalid="ignore", divide="ignore"):
        delta = (np.sqrt(sse[:, draws].sum(axis=2) / n_rows)
                 - np.sqrt(sse_base[:, draws].sum(axis=2) / n_rows)).mean(axis=0)
    delta = delta[np.isfinite(delta)]
    low, high = np.percentile(delta, [2.5, 97.5])
    return float(low), float(high), float((delta >= 0).mean())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, default=GATE_CSV,
                        help="results CSV of the gate run; predictions are read from beside it")
    parser.add_argument("--baseline", default="base", help="config the others are compared with")
    args = parser.parse_args()

    runs = pd.read_csv(args.results)
    runs = runs[runs["protocol"] == "grouped"]
    if args.baseline not in set(runs["config"]):
        raise SystemExit(f"no grouped-split '{args.baseline}' rows in {args.results}")

    y_used = load_pair_rows()
    pair = build_pair_features(y_used)
    y = y_used[COL_TARGET].to_numpy(dtype=DTYPE)
    groups = y_used[COL_CELL_LINE].to_numpy()
    train_idx, _, _ = grouped_split(groups)
    drug_mean = y_used.iloc[train_idx].groupby(COL_DRUG)[COL_TARGET].mean()
    floor_pred = y_used[COL_DRUG].map(drug_mean).fillna(y[train_idx].mean()).to_numpy()

    # error[fold][config][seed] = prediction - truth on that fold's rows
    error: dict = {fold: {} for fold in FOLDS}
    fold_idx: dict = {}
    for config, seed in zip(runs["config"], runs["seed"]):
        saved = np.load(predictions_dir(args.results) / f"{config}__grouped__seed{seed}.npz")
        for fold in FOLDS:
            idx = saved[f"{fold}_idx"]
            if not np.array_equal(saved[f"{fold}_true"], y[idx]):
                raise AssertionError(f"{config} seed {seed}: saved {fold} rows are not the rows rebuilt here")
            fold_idx[fold] = idx
            error[fold].setdefault(config, {})[seed] = saved[f"{fold}_pred"] - y[idx]

    subset_rows, paired_rows = [], []
    rng = np.random.default_rng(RANDOM_STATE)
    for fold in FOLDS:
        idx = fold_idx[fold]
        for subset, mask_all in pair.subsets.items():
            mask = mask_all[idx]
            floor = rmse((floor_pred[idx] - y[idx])[mask])
            for config, by_seed in error[fold].items():
                for seed, err in by_seed.items():
                    subset_rows.append({
                        "config": config, "seed": seed, "fold": fold, "subset": subset,
                        "n_pairs": int(mask.sum()), "rmse": rmse(err[mask]),
                        "mean_error": float(err[mask].mean()), "per_drug_mean_rmse": floor,
                    })

            base = error[fold][args.baseline]
            for config, by_seed in error[fold].items():
                seeds = sorted(set(by_seed) & set(base))
                if config == args.baseline or not seeds:
                    continue
                deltas = [rmse(by_seed[s][mask]) - rmse(base[s][mask]) for s in seeds]
                low, high, p_not_better = bootstrap_delta(
                    np.stack([by_seed[s][mask] ** 2 for s in seeds]),
                    np.stack([base[s][mask] ** 2 for s in seeds]),
                    groups[idx][mask], rng,
                )
                paired_rows.append({
                    "config": config, "baseline": args.baseline, "fold": fold, "subset": subset,
                    "n_pairs": int(mask.sum()), "n_seeds": len(seeds),
                    "delta_rmse": float(np.mean(deltas)),
                    "seeds_better": int(np.sum(np.array(deltas) < 0)),
                    "per_seed_delta": " / ".join(f"{d:+.4f}" for d in deltas),
                    "ci_low": low, "ci_high": high, "p_delta_ge_0": p_not_better,
                })

    subsets = pd.DataFrame(subset_rows)
    paired = pd.DataFrame(paired_rows)
    subsets_csv = args.results.with_name(args.results.stem + "_subsets.csv")
    paired_csv = args.results.with_name(args.results.stem + "_paired.csv")
    subsets.to_csv(subsets_csv, index=False)
    paired.to_csv(paired_csv, index=False)

    for fold in FOLDS:
        print("\n" + "=" * 100)
        print(f"{fold.upper()} — RMSE by subset, mean +/- std over seeds (signed error in brackets)")
        print("=" * 100)
        block = subsets[subsets["fold"] == fold]
        configs = list(dict.fromkeys(block["config"]))
        print(f"{'subset':<12}{'pairs':>7}{'drug-mean':>11}" + "".join(f"{c:>34}" for c in configs))
        for subset in pair.subsets:
            rows = block[block["subset"] == subset]
            line = f"{subset:<12}{rows['n_pairs'].iloc[0]:>7}{rows['per_drug_mean_rmse'].iloc[0]:>11.4f}"
            for config in configs:
                r = rows[rows["config"] == config]
                cell = f"{r['rmse'].mean():.4f} +/- {r['rmse'].std():.4f} [{r['mean_error'].mean():+.3f}]"
                line += f"{cell:>34}"
            print(line)

        for config, block in paired[paired["fold"] == fold].groupby("config", sort=False):
            print(f"\n{config} - {args.baseline}, {fold} (negative = better):")
            print(f"  {'subset':<12}{'pairs':>7}{'delta':>10}{'better':>8}{'95% CI':>22}{'P(d>=0)':>9}   per seed")
            for _, r in block.iterrows():
                ci = f"{r['ci_low']:+.4f} to {r['ci_high']:+.4f}"
                print(f"  {r['subset']:<12}{r['n_pairs']:>7}{r['delta_rmse']:>+10.4f}"
                      f"{str(r['seeds_better']) + '/' + str(r['n_seeds']):>8}{ci:>22}"
                      f"{r['p_delta_ge_0']:>9.3f}   {r['per_seed_delta']}")

    print(f"\nSaved -> {subsets_csv}\nSaved -> {paired_csv}")


if __name__ == "__main__":
    main()
