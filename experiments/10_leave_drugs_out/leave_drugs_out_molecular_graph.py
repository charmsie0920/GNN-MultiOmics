"""leave_drugs_out_molecular_graph.py — do molecular graphs beat fingerprints on unseen compounds?

`leave_drugs_out.py` compared one-hot identity against Morgan fingerprints on
drugs held out entirely from training, and found the gap the cell-line-grouped
matrix could not show: one-hot collapses to R^2 0.024 while fingerprints hold at
0.422. This adds the third representation -- the atom-level molecular graph from
[`11_molecular_graph_results.md`](../../docs/11_molecular_graph_results.md) --
on the same protocol.

**Why this is the test that matters for structure.** Under the cell-line-grouped
split, drug identity alone accounts for 71% of the reducible error and every
test drug appears in training with thousands of measurements, so a model can
memorize each compound's response profile and the drug *representation* barely
matters. That is exactly what was observed: across five seeds the molecular
graph was indistinguishable from (marginally worse than) a Morgan fingerprint
(docs/13_seed_variance_results.md).

Holding out whole drugs removes the option to memorize. Faced with a molecule it
has never seen, a model can only reason from structure. A fingerprint is a *bag*
of substructures; a molecular graph keeps the topology. If that extra
information is worth anything, this is where it shows.

Baseline to beat, from `leave_drugs_out_results.csv`: CrossAttention +
GE+Proteomics + fingerprint. Omics, fusion, head, split and seed are all held
fixed; only the drug encoder changes.

Run from the repository root:
    python "experiments/10_leave_drugs_out/leave_drugs_out_molecular_graph.py"
"""

from __future__ import annotations

import argparse
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
from src.data.drug_graphs import (  # noqa: E402
    assert_fingerprint_population_parity,
    build_drug_graphs,
    build_graph_pair_tensors,
)
from src.data.experiment_utils import (  # noqa: E402
    COL_DRUG,
    GE_KEY,
    PROTEOMICS_KEY,
    compute_shared_threshold,
    evaluate,
    leave_drugs_out_split,
    load_omics_subset,
    load_targets,
    mean_only_floor,
    peak_rss_gb,
    print_metric_block,
)

RESULTS_CSV = Path("experiments/10_leave_drugs_out/leave_drugs_out_molecular_graph_results.csv")
MODALITIES = [GE_KEY, PROTEOMICS_KEY]  # matches the fingerprint arm being compared against


def load_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


mg_mod = load_module("mg_matrix", "experiments/11_molecular_graph/molecular_graph_matrix.py")


def describe_split(y_used: pd.DataFrame, train_idx, val_idx, test_idx) -> dict:
    """Confirm the held-out drugs really are unseen, and report how many."""
    drugs = y_used[COL_DRUG].to_numpy()
    train_drugs, test_drugs = set(drugs[train_idx]), set(drugs[test_idx])
    overlap = train_drugs & test_drugs
    if overlap:
        raise AssertionError(
            f"{len(overlap)} drugs appear in both train and test, e.g. {sorted(overlap)[:5]}. "
            f"The leave-drugs-out protocol requires disjoint drug sets."
        )
    print(
        f"[ldo] train {len(train_drugs)} drugs / {len(train_idx)} pairs | "
        f"test {len(test_drugs)} drugs / {len(test_idx)} pairs | 0 drugs shared"
    )
    return {
        "n_train_drugs": len(train_drugs),
        "n_test_drugs": len(test_drugs),
        "n_train_pairs": len(train_idx),
        "n_test_pairs": len(test_idx),
    }


def run_molecular_graph(graphs, threshold: float, device, seed: int) -> dict:
    print("\n" + "#" * 82)
    print(f"# LEAVE-DRUGS-OUT | CrossAttention | molecular_graph | seed={seed}")
    print("#" * 82)

    torch.manual_seed(seed)
    omics, cell_ids = load_omics_subset(MODALITIES)
    y_df = load_targets(cell_ids)
    gathered, codes, y, _, batched, y_used = build_graph_pair_tensors(
        omics, cell_ids, y_df, graphs
    )
    keys = list(gathered.keys())

    train_idx, val_idx, test_idx = leave_drugs_out_split(y_used[COL_DRUG].to_numpy())
    info = describe_split(y_used, train_idx, val_idx, test_idx)

    omics_tr = {k: v[train_idx] for k, v in gathered.items()}
    omics_va = {k: v[val_idx] for k, v in gathered.items()}
    omics_te = {k: v[test_idx] for k, v in gathered.items()}

    tensors = [torch.from_numpy(omics_tr[k]) for k in keys] + [
        torch.from_numpy(codes[train_idx]), torch.from_numpy(y[train_idx])
    ]
    loader = DataLoader(TensorDataset(*tensors), batch_size=mg_mod.BATCH_SIZE,
                        shuffle=True, drop_last=True)

    # The full drug table is encoded every forward pass, held-out molecules
    # included. That is not leakage: their *structures* are inputs, while their
    # IC50 labels never enter training -- exactly as a fingerprint for an unseen
    # compound would be computed from its SMILES.
    model = mg_mod.CrossAttentionGraphRegressor(MODALITIES, batched).to(device)
    t0 = time.perf_counter()
    model, best_epoch = mg_mod.train(
        model, loader, omics_va, codes[val_idx], y[val_idx], keys, device
    )
    t_fit = time.perf_counter() - t0

    val = evaluate(y[val_idx], mg_mod.predict(model, omics_va, codes[val_idx], keys, device), threshold)
    test = evaluate(y[test_idx], mg_mod.predict(model, omics_te, codes[test_idx], keys, device), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])

    print_metric_block("LDO | CrossAttention | molecular_graph", val, test, floor)
    print(f"params={sum(p.numel() for p in model.parameters()):,} best_epoch={best_epoch}")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {"model": "CrossAttention", "drug_rep": "molecular_graph", "seed": seed,
            "n_pairs": len(y), "mean_only_rmse": floor, **info,
            **{f"test_{k}": v for k, v in test.items()},
            **{f"val_{k}": v for k, v in val.items()},
            "params": sum(p.numel() for p in model.parameters()),
            "best_epoch": best_epoch, "fit_seconds": t_fit}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """`argv` is explicit so a notebook can call `main([...])` without argparse
    picking up the kernel's own `-f kernel.json` argument."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--seeds", nargs="+", type=int, default=[mg_mod.TORCH_SEED],
        help=(
            "Seeds to repeat over. A single run is not enough to separate arms "
            "0.02 apart -- see docs/13_seed_variance_results.md."
        ),
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    t_start = time.perf_counter()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    threshold = compute_shared_threshold()

    graphs = build_drug_graphs()
    assert_fingerprint_population_parity(graphs)

    results = [run_molecular_graph(graphs, threshold, device, seed) for seed in args.seeds]

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"LEAVE-DRUGS-OUT MOLECULAR GRAPH COMPLETE — {len(df)} runs in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    print(df[["seed", "n_test_drugs", "test_rmse", "test_r2", "test_pcc"]].to_string(index=False))
    if len(df) > 1:
        print(f"\nmean test RMSE {df['test_rmse'].mean():.4f} +/- {df['test_rmse'].std(ddof=1):.4f}")
    print("\nCompare against the fingerprint arm in leave_drugs_out_results.csv "
          "(R^2 0.422) and the one-hot arm (R^2 0.024).")
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
