"""seed_variance.py — how much of a result is the architecture, and how much is the seed?

Every row in docs/results.md is a single run at `TORCH_SEED = 42`. That is
fine for ranking configurations that differ by a lot, but the differences
this project now turns on are small:

- molecular graph (1.3021) vs fingerprint (1.3205) -- 0.018 apart
- full architecture (1.3397) vs no PPI graph (1.3021) -- 0.038 apart, while
  their *validation* RMSEs sat 0.003 apart

A gap that only opens at test time is what unstable single-seed results look
like, so neither claim is safe to report until the run-to-run spread is known.

This re-runs the three configurations those claims depend on across several
seeds and reports mean +/- std. The split is **not** reseeded: `grouped_split`
uses the fixed `RANDOM_STATE = 42`, so every run sees the identical
partition and the only thing varying is weight initialization and batch
shuffling. That is deliberate -- it isolates training variance, which is the
quantity needed to judge whether a 0.018 gap is signal.

Because the split is shared, the comparison is **paired**: at a given seed,
two configurations differ only in architecture, so the per-seed difference is
a cleaner estimate than comparing two independent means. Both are reported.

Reuses the existing runners rather than reimplementing them, so a swept run is
the same code path as the original matrix run, with `TORCH_SEED` overridden.

Run from the repository root:
    python "experiments/13_seed_variance/seed_variance.py"
    python "experiments/13_seed_variance/seed_variance.py" --seeds 42 43 44 --configs MG E10
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path

import pandas as pd
import torch

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src.data.drug_graphs import build_drug_graphs  # noqa: E402
from src.data.experiment_utils import (  # noqa: E402
    GE_KEY,
    PROTEOMICS_KEY,
    compute_shared_threshold,
)

RESULTS_CSV = Path("experiments/13_seed_variance/seed_variance_results.csv")
SUMMARY_CSV = Path("experiments/13_seed_variance/seed_variance_summary.csv")

DEFAULT_SEEDS = (42, 43, 44, 45, 46)


def load_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ca_mod = load_module("ca_matrix", "experiments/06_full_matrix/cross_attention_matrix.py")
mlp_mod = load_module("mlp_matrix", "experiments/06_full_matrix/mlp_matrix.py")
mg_mod = load_module("mg_matrix", "experiments/11_molecular_graph/molecular_graph_matrix.py")

# The three configurations the current claims rest on. Labels match
# docs/results.md row IDs where one exists.
CONFIGS: dict[str, dict] = {
    "E04": {
        "label": "MLP | Proteomics | fingerprint",
        "module": mlp_mod,
        "modalities": [PROTEOMICS_KEY],
        "recorded_rmse": 1.2843,
    },
    "E10": {
        "label": "CrossAttention | GE+Proteomics | fingerprint",
        "module": ca_mod,
        "modalities": [GE_KEY, PROTEOMICS_KEY],
        "recorded_rmse": 1.3205,
    },
    "MG": {
        "label": "CrossAttention | GE+Proteomics | molecular_graph",
        "module": mg_mod,
        "modalities": [GE_KEY, PROTEOMICS_KEY],
        "recorded_rmse": 1.3021,
    },
}


def run_config(config_id: str, seed: int, graphs, threshold: float, device) -> dict:
    """Run one configuration at one seed, via its original runner."""
    config = CONFIGS[config_id]
    module = config["module"]

    # Both runners call `torch.manual_seed(TORCH_SEED)` inside `run_one`,
    # reading the module global -- so overriding it here is what reseeds the
    # run, without touching the runner's source.
    previous_seed = module.TORCH_SEED
    module.TORCH_SEED = seed
    try:
        if config_id == "MG":
            row = module.run_one(config["modalities"], graphs, threshold, device)
        else:
            row = module.run_one(
                config["modalities"], "fingerprint", "fingerprint", False, threshold, device
            )
    finally:
        module.TORCH_SEED = previous_seed

    return {"config": config_id, "config_label": config["label"], "seed": seed, **row}


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    """Mean/std/min/max of test RMSE and R2 per configuration."""
    rows = []
    for config_id, group in df.groupby("config", sort=False):
        rows.append({
            "config": config_id,
            "label": CONFIGS[config_id]["label"],
            "n_seeds": len(group),
            "recorded_rmse": CONFIGS[config_id]["recorded_rmse"],
            "rmse_mean": group["test_rmse"].mean(),
            "rmse_std": group["test_rmse"].std(ddof=1),
            "rmse_min": group["test_rmse"].min(),
            "rmse_max": group["test_rmse"].max(),
            "r2_mean": group["test_r2"].mean(),
            "r2_std": group["test_r2"].std(ddof=1),
        })
    return pd.DataFrame(rows).sort_values("rmse_mean").reset_index(drop=True)


def paired_deltas(df: pd.DataFrame) -> pd.DataFrame:
    """Per-seed RMSE differences between configuration pairs.

    All configurations share one split, so at a fixed seed two runs differ
    only in architecture. The spread of that per-seed difference is a tighter
    test than whether two independent means are separated.
    """
    wide = df.pivot(index="seed", columns="config", values="test_rmse")
    rows = []
    for a, b in [("MG", "E10"), ("MG", "E04"), ("E10", "E04")]:
        if a not in wide.columns or b not in wide.columns:
            continue
        delta = wide[a] - wide[b]
        mean, std = delta.mean(), delta.std(ddof=1)
        rows.append({
            "comparison": f"{a} - {b}",
            "mean_delta": mean,
            "std_delta": std,
            "n_seeds": len(delta),
            "wins": int((delta < 0).sum()),  # negative delta = `a` has lower RMSE
            # A difference smaller than its own spread is not separable at this
            # sample size; this is a readability aid, not a significance test.
            "separable": bool(abs(mean) > std) if len(delta) > 1 else False,
        })
    return pd.DataFrame(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """`argv` is explicit so a notebook can call `main([...])` without argparse
    picking up the kernel's own `-f kernel.json` argument."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", nargs="+", type=int, default=list(DEFAULT_SEEDS))
    parser.add_argument(
        "--configs", nargs="+", choices=sorted(CONFIGS), default=sorted(CONFIGS),
        help="Configurations to sweep. Defaults to all three.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    t_start = time.perf_counter()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    print(f"[sweep] configs={args.configs}  seeds={args.seeds}  "
          f"= {len(args.configs) * len(args.seeds)} runs")
    threshold = compute_shared_threshold()

    # Only needed by the molecular-graph config, but parsing SMILES once and
    # sharing keeps every run on an identical drug set.
    graphs = build_drug_graphs() if "MG" in args.configs else None

    results = []
    for seed in args.seeds:
        for config_id in args.configs:
            print(f"\n{'=' * 82}\n[seed {seed}] {CONFIGS[config_id]['label']}\n{'=' * 82}")
            results.append(run_config(config_id, seed, graphs, threshold, device))

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    summary = summarize(df)
    summary.to_csv(SUMMARY_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"SEED VARIANCE COMPLETE — {len(df)} runs in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    print("\nPer-configuration spread (test RMSE):")
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    if len(args.seeds) > 1:
        print("\nPaired per-seed differences (negative = first config better):")
        print(paired_deltas(df).to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    print(f"\nSaved -> {RESULTS_CSV}\n      -> {SUMMARY_CSV}")


if __name__ == "__main__":
    main()
