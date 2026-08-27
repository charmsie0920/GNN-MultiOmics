"""split_comparison.py — how much of MoGraphDRP's advantage is the split protocol?

Every result in this project is measured under a cell-line-grouped split, so a
test cell line is never seen in training. MoGraphDRP reports RMSE 0.6622 under
a *random* 80/10/10 split over (cell line, drug) pairs, where a cell line's
~279 measurements scatter across all three folds -- meaning at test time their
model has already seen ~223 other drug responses for that same cell line.

Those two numbers measure different tasks (imputation vs. generalization to an
unseen cell line), so they cannot be compared directly. This script produces
the missing like-for-like number by running our best model under BOTH
protocols, everything else held fixed.

It also tests the leakage hypothesis from docs/ensemble_refinement_results.md:
XGBoost residual refinement gave the paper +19.7% but cost us -0.9%. If the
explanation is cell-line leakage -- a refiner can learn "this cell line reads
0.3 high" only if it has seen that cell line -- then refinement should start
helping under the random split. That is a falsifiable prediction, and this
script checks it.

Run from the repository root:
    python "experiments/Split Protocol/split_comparison.py"
"""

from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, TensorDataset

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.data.experiment_utils import (  # noqa: E402
    GE_KEY,
    PROTEOMICS_KEY,
    build_pair_tensors,
    compute_shared_threshold,
    evaluate,
    grouped_split,
    leakage_report,
    load_omics_subset,
    load_targets,
    mean_only_floor,
    peak_rss_gb,
    print_metric_block,
    random_pair_split,
)
from src.models.ensemble_refinement import fit_refiner, refine  # noqa: E402

RESULTS_CSV = Path("experiments/Split Protocol/split_comparison_results.csv")

# Best cell in the matrix (docs/results.md E01).
MODALITIES = [GE_KEY, PROTEOMICS_KEY]
DRUG_MODE = "onehot"


def load_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ca_mod = load_module("ca_matrix", "experiments/Full Matrix/cross_attention_matrix.py")
ens_mod = load_module("ens_base", "experiments/Ensemble Refinement/ensemble_baseline.py")


def run_protocol(protocol: str, gathered, drug_block, y, groups, n_drug_features, threshold, device) -> dict:
    """Train the best model under one split protocol, then also try XGBoost refinement."""
    print("\n" + "#" * 82)
    print(f"# SPLIT PROTOCOL: {protocol}")
    print("#" * 82)

    torch.manual_seed(ca_mod.TORCH_SEED)
    keys = list(gathered.keys())

    if protocol == "grouped_by_cell_line":
        train_idx, val_idx, test_idx = grouped_split(groups)
    elif protocol == "random_pairs":
        train_idx, val_idx, test_idx = random_pair_split(len(y))
    else:
        raise ValueError(protocol)

    leak = leakage_report(groups, train_idx, test_idx)
    print(f"[leakage] test cell lines: {leak['test_cell_lines']}, "
          f"of which also in train: {leak['also_in_train']}")
    print(f"[leakage] test rows whose cell line was seen in training: "
          f"{leak['test_rows_with_seen_cell_line_pct']:.1f}%")
    print(f"[leakage] mean training rows per test cell line: "
          f"{leak['mean_train_rows_per_test_cell_line']:.0f}")

    omics_tr = {k: v[train_idx] for k, v in gathered.items()}
    omics_va = {k: v[val_idx] for k, v in gathered.items()}
    omics_te = {k: v[test_idx] for k, v in gathered.items()}

    tensors = [torch.from_numpy(omics_tr[k]) for k in keys] + [
        torch.from_numpy(drug_block[train_idx]),
        torch.from_numpy(y[train_idx]),
    ]
    loader = DataLoader(
        TensorDataset(*tensors), batch_size=ca_mod.BATCH_SIZE, shuffle=True, drop_last=True
    )

    model = ca_mod.CrossAttentionRegressor(MODALITIES, n_drug_features, DRUG_MODE).to(device)
    t0 = time.perf_counter()
    model, best_epoch = ca_mod.train(
        model, loader, omics_va, drug_block[val_idx], y[val_idx], keys, device
    )
    t_fit = time.perf_counter() - t0

    capture = ens_mod.ActivationCapture(model.head)

    def predict_with_activations(omics_split, drug_split):
        pred = ca_mod.predict(model, omics_split, drug_split, keys, device)
        return pred, capture.drain()

    pred_tr, act_tr = predict_with_activations(omics_tr, drug_block[train_idx])
    pred_va, act_va = predict_with_activations(omics_va, drug_block[val_idx])
    pred_te, act_te = predict_with_activations(omics_te, drug_block[test_idx])
    capture.close()

    base_val = evaluate(y[val_idx], pred_va, threshold)
    base_test = evaluate(y[test_idx], pred_te, threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])
    print_metric_block(f"{protocol} | base (no XGBoost)", base_val, base_test, floor)

    refiner = fit_refiner(act_tr, pred_tr, y[train_idx])
    ref_val = evaluate(y[val_idx], refine(refiner, act_va, pred_va), threshold)
    ref_test = evaluate(y[test_idx], refine(refiner, act_te, pred_te), threshold)
    print_metric_block(f"{protocol} | + XGBoost refinement", ref_val, ref_test, floor)

    delta = base_test["rmse"] - ref_test["rmse"]
    pct = 100 * delta / base_test["rmse"]
    print(f"\n[refinement] test RMSE {base_test['rmse']:.4f} -> {ref_test['rmse']:.4f} "
          f"({delta:+.4f}, {pct:+.1f}%)")
    print(f"[timing] fit={t_fit:.1f}s  best_epoch={best_epoch}")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {
        "protocol": protocol,
        "model": "CrossAttention",
        "omics": "+".join(MODALITIES),
        "drug_rep": DRUG_MODE,
        "n_pairs": len(y),
        "n_train": len(train_idx),
        "n_test": len(test_idx),
        "mean_only_rmse": floor,
        **{f"leak_{k}": v for k, v in leak.items()},
        **{f"base_test_{k}": v for k, v in base_test.items()},
        **{f"refined_test_{k}": v for k, v in ref_test.items()},
        "refinement_delta": delta,
        "refinement_pct": pct,
        "best_epoch": best_epoch,
        "fit_seconds": t_fit,
    }


def main() -> None:
    t_start = time.perf_counter()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    threshold = compute_shared_threshold()

    omics, cell_ids = load_omics_subset(MODALITIES)
    y_df = load_targets(cell_ids)
    gathered, drug_block, y, groups, n_drug_features, _ = build_pair_tensors(
        omics, cell_ids, y_df, DRUG_MODE, False
    )

    results = [
        run_protocol("grouped_by_cell_line", gathered, drug_block, y, groups,
                     n_drug_features, threshold, device),
        run_protocol("random_pairs", gathered, drug_block, y, groups,
                     n_drug_features, threshold, device),
    ]

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    grouped = df[df.protocol == "grouped_by_cell_line"].iloc[0]
    random_ = df[df.protocol == "random_pairs"].iloc[0]

    print("\n" + "=" * 100)
    print(f"SPLIT PROTOCOL COMPARISON COMPLETE — {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    print(df[["protocol", "base_test_rmse", "base_test_pcc", "base_test_r2",
              "refined_test_rmse", "refinement_pct"]].round(4).to_string(index=False))
    print()
    print(f"Protocol effect (grouped -> random): RMSE "
          f"{grouped['base_test_rmse']:.4f} -> {random_['base_test_rmse']:.4f} "
          f"({grouped['base_test_rmse'] - random_['base_test_rmse']:+.4f})")
    print(f"MoGraphDRP published (random split): RMSE 0.6622 / PCC 0.9689 / R2 0.9388")
    print()
    print("Leakage hypothesis for XGBoost refinement:")
    print(f"  grouped split : {grouped['refinement_pct']:+.2f}%")
    print(f"  random split  : {random_['refinement_pct']:+.2f}%")
    verdict = ("SUPPORTED — refinement helps only when cell lines leak"
               if random_["refinement_pct"] > 0 >= grouped["refinement_pct"]
               else "NOT SUPPORTED — refinement behaves the same under both protocols")
    print(f"  verdict: {verdict}")
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
