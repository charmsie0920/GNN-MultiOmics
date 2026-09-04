"""ensemble_baseline.py — XGBoost residual refinement on the matrix's best models.

Phase 4 of the experiment matrix (docs/plan/experiment_matrix_plan.md), and
the last architectural component from MoGraphDRP we hadn't tested: does
bolting a gradient-boosted residual corrector onto a trained neural model
recover further accuracy, as the paper reports (their Table 7: ~19.7% RMSE
reduction)?

Applied to the two winners from Phases 1-3, chosen after those results landed:
  - CrossAttention, GE+Proteomics, one-hot  (test RMSE 1.2442 — best overall)
  - GNN-GCN, tri-omics, +mutation edges     (test RMSE 1.3513 — best graph model)

Both base models are retrained here from the same seed/split as their matrix
run, then their penultimate head activation (the analogue of the paper's
128-dim interaction vector `f_ij`) and initial prediction are handed to the
refiner. The refiner is fit on the training fold only.

Run from the repository root:
    python "experiments/08_ensemble_refinement/ensemble_baseline.py"
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
    MUT_CNV_KEY,
    PROTEOMICS_KEY,
    build_pair_tensors,
    compute_shared_threshold,
    evaluate,
    grouped_split,
    load_omics_subset,
    load_targets,
    mean_only_floor,
    peak_rss_gb,
    print_metric_block,
)
from src.models.ensemble_refinement import feature_importance, fit_refiner, refine  # noqa: E402

RESULTS_CSV = Path("experiments/08_ensemble_refinement/ensemble_results.csv")


def load_module(name: str, relative_path: str):
    """Import a module from a path (the experiment dirs have spaces in their names)."""
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ca_mod = load_module("ca_matrix", "experiments/06_full_matrix/cross_attention_matrix.py")
gnn_mod = load_module("gnn_base", "experiments/07_gnn_ablation/gnn_baseline.py")


class ActivationCapture:
    """Forward hook capturing the input to the head's final Linear layer.

    That input is the last learned representation before the scalar readout --
    our stand-in for MoGraphDRP's compressed interaction vector `f_ij`.
    """

    def __init__(self, head: torch.nn.Module):
        self.buffer: list[np.ndarray] = []
        final_linear = [m for m in head if isinstance(m, torch.nn.Linear)][-1]
        self.handle = final_linear.register_forward_hook(self._hook)

    def _hook(self, _module, inputs, _output):
        self.buffer.append(inputs[0].detach().cpu().numpy())

    def drain(self) -> np.ndarray:
        out = np.concatenate(self.buffer, axis=0)
        self.buffer.clear()
        return out

    def close(self):
        self.handle.remove()


def run_cross_attention(threshold: float, device) -> dict:
    """Retrain the best cross-attention cell, then refine its predictions."""
    modalities = [GE_KEY, PROTEOMICS_KEY]
    print("\n" + "#" * 82)
    print("# ENSEMBLE | base=CrossAttention | omics=GE+Proteomics | drug=onehot")
    print("#" * 82)

    torch.manual_seed(ca_mod.TORCH_SEED)
    omics, cell_ids = load_omics_subset(modalities)
    y_df = load_targets(cell_ids)
    gathered, drug_block, y, groups, n_drug_features, _ = build_pair_tensors(
        omics, cell_ids, y_df, "onehot", False
    )
    keys = list(gathered.keys())
    train_idx, val_idx, test_idx = grouped_split(groups)

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

    model = ca_mod.CrossAttentionRegressor(modalities, n_drug_features, "onehot").to(device)
    t0 = time.perf_counter()
    model, best_epoch = ca_mod.train(
        model, loader, omics_va, drug_block[val_idx], y[val_idx], keys, device
    )
    t_base = time.perf_counter() - t0

    capture = ActivationCapture(model.head)

    def predict_with_activations(omics_split, drug_split):
        pred = ca_mod.predict(model, omics_split, drug_split, keys, device)
        return pred, capture.drain()

    pred_tr, act_tr = predict_with_activations(omics_tr, drug_block[train_idx])
    pred_va, act_va = predict_with_activations(omics_va, drug_block[val_idx])
    pred_te, act_te = predict_with_activations(omics_te, drug_block[test_idx])
    capture.close()

    return _refine_and_report(
        base_name="CrossAttention",
        config="GE+Proteomics | onehot",
        y=y, train_idx=train_idx, val_idx=val_idx, test_idx=test_idx,
        pred_tr=pred_tr, pred_va=pred_va, pred_te=pred_te,
        act_tr=act_tr, act_va=act_va, act_te=act_te,
        threshold=threshold, best_epoch=best_epoch, t_base=t_base,
    )


def run_gnn(threshold: float, device) -> dict:
    """Retrain the best GNN cell (GCN + mutation edges), then refine its predictions."""
    print("\n" + "#" * 82)
    print("# ENSEMBLE | base=GNN-GCN | omics=tri | drug=fingerprint | +mutation_edges")
    print("#" * 82)

    torch.manual_seed(gnn_mod.TORCH_SEED)
    data, cell_index, drug_index, y, groups = gnn_mod.load_graph_and_pairs(threshold)

    x_dict = {nt: data[nt].x.to(device) for nt in data.node_types}
    edge_dict = {et: data[et].edge_index.to(device) for et in data.edge_types}
    train_idx, val_idx, test_idx = grouped_split(groups)

    ds = TensorDataset(
        torch.from_numpy(cell_index[train_idx]),
        torch.from_numpy(drug_index[train_idx]),
        torch.from_numpy(y[train_idx]),
    )
    loader = DataLoader(ds, batch_size=gnn_mod.BATCH_SIZE, shuffle=True, drop_last=True)

    model = gnn_mod.HeteroGNN(
        metadata=(list(data.node_types), list(edge_dict.keys())),
        cell_line_dim=data["cell_line"].x.shape[1],
        drug_dim=data["drug"].x.shape[1],
        num_proteins=data["protein"].x.shape[0],
        variant="gcn",
        hidden_dim=gnn_mod.HIDDEN_DIM,
        num_layers=gnn_mod.NUM_LAYERS,
        heads=gnn_mod.HEADS,
        dropout=gnn_mod.DROPOUT,
        head_hidden_dims=gnn_mod.HEAD_HIDDEN_DIMS,
        head_dropout=gnn_mod.HEAD_DROPOUT,
    ).to(device)

    model.eval()
    with torch.no_grad():  # materialize lazy conv weights before the optimizer sees them
        model(x_dict, edge_dict,
              torch.from_numpy(cell_index[train_idx[:2]]).to(device),
              torch.from_numpy(drug_index[train_idx[:2]]).to(device))

    t0 = time.perf_counter()
    model, best_epoch = gnn_mod.train(
        model, x_dict, edge_dict, loader,
        cell_index[val_idx], drug_index[val_idx], y[val_idx], device,
    )
    t_base = time.perf_counter() - t0

    capture = ActivationCapture(model.head)

    def predict_with_activations(idx):
        pred = gnn_mod.predict(model, x_dict, edge_dict, cell_index[idx], drug_index[idx], device)
        return pred, capture.drain()

    pred_tr, act_tr = predict_with_activations(train_idx)
    pred_va, act_va = predict_with_activations(val_idx)
    pred_te, act_te = predict_with_activations(test_idx)
    capture.close()

    return _refine_and_report(
        base_name="GNN-GCN",
        config="tri-omics | fingerprint | +mutation_edges",
        y=y, train_idx=train_idx, val_idx=val_idx, test_idx=test_idx,
        pred_tr=pred_tr, pred_va=pred_va, pred_te=pred_te,
        act_tr=act_tr, act_va=act_va, act_te=act_te,
        threshold=threshold, best_epoch=best_epoch, t_base=t_base,
    )


def _refine_and_report(
    base_name, config, y, train_idx, val_idx, test_idx,
    pred_tr, pred_va, pred_te, act_tr, act_va, act_te,
    threshold, best_epoch, t_base,
) -> dict:
    """Fit the refiner on the train fold and report before/after on val + test."""
    base_val = evaluate(y[val_idx], pred_va, threshold)
    base_test = evaluate(y[test_idx], pred_te, threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])

    print_metric_block(f"{base_name} | {config} | WITHOUT XGBoost", base_val, base_test, floor)

    t0 = time.perf_counter()
    refiner = fit_refiner(act_tr, pred_tr, y[train_idx])
    t_refine = time.perf_counter() - t0

    ref_val = evaluate(y[val_idx], refine(refiner, act_va, pred_va), threshold)
    ref_test = evaluate(y[test_idx], refine(refiner, act_te, pred_te), threshold)

    print_metric_block(f"{base_name} | {config} | WITH XGBoost", ref_val, ref_test, floor)

    delta = base_test["rmse"] - ref_test["rmse"]
    pct = 100 * delta / base_test["rmse"]
    print(f"\n[refinement] test RMSE {base_test['rmse']:.4f} -> {ref_test['rmse']:.4f} "
          f"({delta:+.4f}, {pct:+.1f}%)")
    print(f"[refinement] interaction vector dim={act_tr.shape[1]}  "
          f"base_fit={t_base:.1f}s  refine_fit={t_refine:.1f}s")
    print(f"[refinement] top features: {feature_importance(refiner, top_n=6)}")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {
        "base_model": base_name,
        "config": config,
        "n_pairs": len(y),
        "interaction_dim": act_tr.shape[1],
        "mean_only_rmse": floor,
        **{f"base_test_{k}": v for k, v in base_test.items()},
        **{f"refined_test_{k}": v for k, v in ref_test.items()},
        **{f"base_val_{k}": v for k, v in base_val.items()},
        **{f"refined_val_{k}": v for k, v in ref_val.items()},
        "rmse_delta": delta,
        "rmse_pct_change": pct,
        "best_epoch": best_epoch,
    }


def main() -> None:
    t_start = time.perf_counter()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    threshold = compute_shared_threshold()

    results = [run_cross_attention(threshold, device), run_gnn(threshold, device)]

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"ENSEMBLE REFINEMENT COMPLETE — {len(df)} base models in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    print(df[["base_model", "config", "base_test_rmse", "refined_test_rmse",
              "rmse_delta", "rmse_pct_change"]].round(4).to_string(index=False))
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
