"""leave_drugs_out.py — can the model predict for compounds it has never seen?

Every other experiment in this project holds out *cell lines*, so a test drug
has always been screened during training. That protocol cannot distinguish the
two drug representations on the axis that actually matters for a discovery
tool: one-hot identity and Morgan fingerprints both work fine when the drug is
in the training vocabulary.

This script holds out *drugs* instead. ~44 of the 295 compounds never appear in
training, and the model must score them cold. The two representations behave
completely differently here:

  - one-hot: a held-out drug's column is never activated during training, so its
    weight stays at initialization. The model has no usable information about
    the compound and can only fall back on cell-line features. It is not that
    one-hot predicts badly -- it structurally cannot represent the drug.
  - fingerprint: a held-out drug still has 2048 structural bits computed from
    its SMILES, so if it shares substructures with training compounds the model
    can transfer what it learned about them.

This is the experiment that justifies (or refutes) the Morgan fingerprint
choice in the final architecture -- see docs/cross_attention_ablation_results.md
§3, where fingerprints look *worse* under the cell-line-grouped split.

Both arms run on the fingerprint-resolvable population so the comparison is not
confounded by differing row sets (same control as the main matrix).

Run from the repository root:
    python "experiments/Leave Drugs Out/leave_drugs_out.py"
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
from sklearn.ensemble import RandomForestRegressor
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.data.experiment_utils import (  # noqa: E402
    COL_DRUG,
    GE_KEY,
    PROTEOMICS_KEY,
    RANDOM_STATE,
    build_pair_matrix,
    build_pair_tensors,
    compute_shared_threshold,
    concat_omics,
    evaluate,
    leave_drugs_out_split,
    load_omics_subset,
    load_targets,
    mean_only_floor,
    peak_rss_gb,
    print_metric_block,
)

RESULTS_CSV = Path("experiments/Leave Drugs Out/leave_drugs_out_results.csv")

MODALITIES = [GE_KEY, PROTEOMICS_KEY]  # best omics subset from the matrix (E01)
ARMS = ["onehot", "fingerprint"]


def load_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ca_mod = load_module("ca_matrix", "experiments/Full Matrix/cross_attention_matrix.py")
mlp_mod = load_module("mlp_matrix", "experiments/Full Matrix/mlp_matrix.py")


def describe_split(y_used: pd.DataFrame, train_idx, val_idx, test_idx) -> dict:
    """Confirm the held-out drugs really are unseen, and report how many."""
    drugs = y_used[COL_DRUG].to_numpy()
    train_drugs, test_drugs = set(drugs[train_idx]), set(drugs[test_idx])
    overlap = train_drugs & test_drugs
    print(f"[drugs]   train={len(train_drugs)}  val={len(set(drugs[val_idx]))}  "
          f"test={len(test_drugs)}  overlap={len(overlap)}")
    return {
        "train_drugs": len(train_drugs),
        "test_drugs": len(test_drugs),
        "drug_overlap": len(overlap),
    }


def run_rf(arm: str, threshold: float) -> dict:
    print("\n" + "#" * 82)
    print(f"# LEAVE-DRUGS-OUT | RF | drug={arm}")
    print("#" * 82)

    omics, cell_ids = load_omics_subset(MODALITIES)
    y_df = load_targets(cell_ids)
    X_cell = concat_omics(omics, MODALITIES)
    X, y, _, _, y_used = build_pair_matrix(X_cell, cell_ids, y_df, arm, restrict_to_fingerprintable=True)

    drug_groups = y_used[COL_DRUG].to_numpy()
    train_idx, val_idx, test_idx = leave_drugs_out_split(drug_groups)
    info = describe_split(y_used, train_idx, val_idx, test_idx)

    model = RandomForestRegressor(
        n_estimators=100, min_samples_leaf=5, max_features="sqrt", n_jobs=-1, random_state=RANDOM_STATE
    )
    t0 = time.perf_counter()
    model.fit(X[train_idx], y[train_idx])
    t_fit = time.perf_counter() - t0

    val = evaluate(y[val_idx], model.predict(X[val_idx]), threshold)
    test = evaluate(y[test_idx], model.predict(X[test_idx]), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])
    print_metric_block(f"LDO | RF | {arm}", val, test, floor)
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {"model": "RF", "drug_rep": arm, "n_pairs": len(y), "mean_only_rmse": floor,
            **info, **{f"test_{k}": v for k, v in test.items()},
            **{f"val_{k}": v for k, v in val.items()}, "fit_seconds": t_fit}


def run_mlp(arm: str, threshold: float, device) -> dict:
    print("\n" + "#" * 82)
    print(f"# LEAVE-DRUGS-OUT | MLP | drug={arm}")
    print("#" * 82)

    torch.manual_seed(mlp_mod.TORCH_SEED)
    omics, cell_ids = load_omics_subset(MODALITIES)
    y_df = load_targets(cell_ids)
    X_cell = concat_omics(omics, MODALITIES)
    X, y, _, _, y_used = build_pair_matrix(X_cell, cell_ids, y_df, arm, restrict_to_fingerprintable=True)

    drug_groups = y_used[COL_DRUG].to_numpy()
    train_idx, val_idx, test_idx = leave_drugs_out_split(drug_groups)
    info = describe_split(y_used, train_idx, val_idx, test_idx)

    ds = TensorDataset(torch.from_numpy(X[train_idx]), torch.from_numpy(y[train_idx]))
    loader = DataLoader(ds, batch_size=mlp_mod.BATCH_SIZE, shuffle=True, drop_last=True)

    model = mlp_mod.MLPRegressor(
        in_dim=X.shape[1], hidden_dims=mlp_mod.HIDDEN_DIMS, dropout=mlp_mod.DROPOUT
    ).to(device)
    t0 = time.perf_counter()
    model, best_epoch = mlp_mod.train(model, loader, X[val_idx], y[val_idx], device)
    t_fit = time.perf_counter() - t0

    val = evaluate(y[val_idx], mlp_mod.predict(model, X[val_idx], device), threshold)
    test = evaluate(y[test_idx], mlp_mod.predict(model, X[test_idx], device), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])
    print_metric_block(f"LDO | MLP | {arm}", val, test, floor)
    print(f"params={sum(p.numel() for p in model.parameters()):,} best_epoch={best_epoch}")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {"model": "MLP", "drug_rep": arm, "n_pairs": len(y), "mean_only_rmse": floor,
            **info, **{f"test_{k}": v for k, v in test.items()},
            **{f"val_{k}": v for k, v in val.items()},
            "best_epoch": best_epoch, "fit_seconds": t_fit}


def run_cross_attention(arm: str, threshold: float, device) -> dict:
    print("\n" + "#" * 82)
    print(f"# LEAVE-DRUGS-OUT | CrossAttention | drug={arm}")
    print("#" * 82)

    torch.manual_seed(ca_mod.TORCH_SEED)
    omics, cell_ids = load_omics_subset(MODALITIES)
    y_df = load_targets(cell_ids)
    gathered, drug_block, y, _, n_drug_features, y_used = build_pair_tensors(
        omics, cell_ids, y_df, arm, restrict_to_fingerprintable=True
    )
    keys = list(gathered.keys())

    drug_groups = y_used[COL_DRUG].to_numpy()
    train_idx, val_idx, test_idx = leave_drugs_out_split(drug_groups)
    info = describe_split(y_used, train_idx, val_idx, test_idx)

    omics_tr = {k: v[train_idx] for k, v in gathered.items()}
    omics_va = {k: v[val_idx] for k, v in gathered.items()}
    omics_te = {k: v[test_idx] for k, v in gathered.items()}

    tensors = [torch.from_numpy(omics_tr[k]) for k in keys] + [
        torch.from_numpy(drug_block[train_idx]), torch.from_numpy(y[train_idx])
    ]
    loader = DataLoader(TensorDataset(*tensors), batch_size=ca_mod.BATCH_SIZE,
                        shuffle=True, drop_last=True)

    model = ca_mod.CrossAttentionRegressor(MODALITIES, n_drug_features, arm).to(device)
    t0 = time.perf_counter()
    model, best_epoch = ca_mod.train(
        model, loader, omics_va, drug_block[val_idx], y[val_idx], keys, device
    )
    t_fit = time.perf_counter() - t0

    val = evaluate(y[val_idx], ca_mod.predict(model, omics_va, drug_block[val_idx], keys, device), threshold)
    test = evaluate(y[test_idx], ca_mod.predict(model, omics_te, drug_block[test_idx], keys, device), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])
    print_metric_block(f"LDO | CrossAttention | {arm}", val, test, floor)
    print(f"params={sum(p.numel() for p in model.parameters()):,} best_epoch={best_epoch}")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {"model": "CrossAttention", "drug_rep": arm, "n_pairs": len(y), "mean_only_rmse": floor,
            **info, **{f"test_{k}": v for k, v in test.items()},
            **{f"val_{k}": v for k, v in val.items()},
            "best_epoch": best_epoch, "fit_seconds": t_fit}


def main() -> None:
    t_start = time.perf_counter()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    threshold = compute_shared_threshold()

    results = []
    for arm in ARMS:
        results.append(run_rf(arm, threshold))
    for arm in ARMS:
        results.append(run_mlp(arm, threshold, device))
    for arm in ARMS:
        results.append(run_cross_attention(arm, threshold, device))

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"LEAVE-DRUGS-OUT COMPLETE — {len(df)} runs in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    print(df[["model", "drug_rep", "test_drugs", "mean_only_rmse",
              "test_rmse", "test_pcc", "test_r2"]].round(4).to_string(index=False))

    print("\n--- fingerprint advantage on UNSEEN drugs ---")
    for m in df["model"].unique():
        oh = df[(df.model == m) & (df.drug_rep == "onehot")].iloc[0]
        fp = df[(df.model == m) & (df.drug_rep == "fingerprint")].iloc[0]
        floor = oh["mean_only_rmse"]
        print(f"  {m:<15} onehot={oh['test_rmse']:.4f}  fingerprint={fp['test_rmse']:.4f}  "
              f"delta={oh['test_rmse'] - fp['test_rmse']:+.4f}   (mean-only floor={floor:.4f})")
        print(f"  {'':<15} onehot beats floor by {floor - oh['test_rmse']:+.4f}, "
              f"fingerprint by {floor - fp['test_rmse']:+.4f}")
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
