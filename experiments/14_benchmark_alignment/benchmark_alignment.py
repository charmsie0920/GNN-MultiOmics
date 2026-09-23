"""benchmark_alignment.py — the aligned model under both protocols, side by side.

Answers the question the project could not previously answer fairly: how much
of the gap to MoGraphDRP's published RMSE is **protocol**, how much is
**architecture**, and how much is the data this project does not have?

Two things happen here that have not happened together before:

1. Each configuration is evaluated under **both** the benchmark's random
   pair split and this project's cell-line-grouped split, so their published
   number has a like-for-like counterpart and the honest generalization number
   is reported beside it.
2. The configurations differ from the MoGraphDRP-aligned baseline by **one
   component at a time** (`src/models/mographdrp_aligned.py`), so each
   difference is attributable instead of being one large unexplained gap.

The previous protocol comparison (docs/09_split_protocol_comparison.md) used
row E01 -- the *one-hot* configuration. Under a random split, one-hot plus
memorization is exactly the setup this project criticizes the benchmark for,
so that number could not support a fair comparison. Every configuration here
represents a drug by its structure.

Reference points, all under the random split:

| MoGraphDRP published (with XGBoost)      | 0.6622 |
| MoGraphDRP without XGBoost (Table 6)     | 0.9497 |
| BANDRP, best prior work (their Table 2)  | 0.9305 |

Run from the repository root:
    python "experiments/14_benchmark_alignment/benchmark_alignment.py"
    python "experiments/14_benchmark_alignment/benchmark_alignment.py" --configs aligned --protocols random
"""

from __future__ import annotations

import argparse
import copy
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.drug_graphs import (  # noqa: E402
    assert_fingerprint_population_parity,
    build_drug_graphs,
    collate_drug_graphs,
)
from src.data.experiment_utils import (  # noqa: E402
    COL_DRUG,
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
from src.models.mographdrp_aligned import (  # noqa: E402
    BATCH_SIZE,
    DROPOUT,
    LR,
    MAX_EPOCHS,
    MoGraphDRPAligned,
)

RESULTS_CSV = Path("experiments/14_benchmark_alignment/benchmark_alignment_results.csv")

MODALITIES = [GE_KEY, PROTEOMICS_KEY]
TORCH_SEED = 42
WEIGHT_DECAY = 1e-5
PATIENCE = 15
LR_PATIENCE = 5

# Each entry differs from `aligned` by exactly one switch, so the delta against
# `aligned` measures that one component. `project_current` reproduces E10.
CONFIGS: dict[str, dict] = {
    "aligned": dict(fusion="concat", head="bilinear", drug_mode="both",
                    note="MoGraphDRP-aligned baseline"),
    "aligned_mlp_head": dict(fusion="concat", head="mlp", drug_mode="both",
                             note="-bilinear head (isolates their sec 2.3)"),
    "aligned_fp_only": dict(fusion="concat", head="bilinear", drug_mode="fingerprint",
                            note="-molecular graph (isolates their sec 2.2.3)"),
    "aligned_graph_only": dict(fusion="concat", head="bilinear", drug_mode="graph",
                               note="-fingerprint (isolates their sec 2.2.3)"),
    "cross_attention": dict(fusion="cross_attention", head="bilinear", drug_mode="both",
                            note="+this project's fusion (their sec 2.1 rejects learned fusion)"),
    "project_current": dict(fusion="cross_attention", head="mlp", drug_mode="fingerprint",
                            note="this project's E10 configuration"),
}
PROTOCOLS = ("grouped", "random")


def load_pairs(graphs):
    """Fingerprints and molecular-graph codes on one identical row set.

    Both drug representations are needed simultaneously for `drug_mode="both"`,
    and they must describe the same pairs in the same order. Rather than
    trusting two independent builders to agree, the fingerprint path defines
    the rows and the graph codes are derived from *its* `y_used`.
    """
    omics, cell_ids = load_omics_subset(MODALITIES)
    y_df = load_targets(cell_ids)
    gathered, fingerprints, y, groups, n_fp, y_used = build_pair_tensors(
        omics, cell_ids, y_df, "fingerprint"
    )

    codes, levels = pd.factorize(y_used[COL_DRUG], sort=True)
    missing = [d for d in levels if d not in graphs]
    if missing:
        raise AssertionError(
            f"{len(missing)} drugs have a fingerprint but no molecular graph, "
            f"e.g. {missing[:5]} — the two arms would not be comparable."
        )
    batched = collate_drug_graphs(graphs, list(levels))
    return gathered, fingerprints, codes.astype(np.int64), y, groups, batched, n_fp, y_used


@torch.no_grad()
def predict(model, omics, fingerprints, codes, keys, device) -> np.ndarray:
    model.eval()
    preds = []
    for start in range(0, len(codes), BATCH_SIZE):
        end = start + BATCH_SIZE
        ob = {k: torch.from_numpy(omics[k][start:end]).to(device) for k in keys}
        fb = torch.from_numpy(fingerprints[start:end]).to(device)
        cb = torch.from_numpy(codes[start:end]).to(device)
        preds.append(model(ob, fb, cb).cpu().numpy())
    return np.concatenate(preds)


def train(model, loader, omics_va, fp_va, codes_va, y_va, keys, device) -> tuple[nn.Module, int]:
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=LR_PATIENCE
    )
    best_state = copy.deepcopy(model.state_dict())
    best_rmse, best_epoch, stale = float("inf"), 0, 0

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        for *omics_batches, fb, cb, yb in loader:
            od = {k: t.to(device) for k, t in zip(keys, omics_batches)}
            fb, cb, yb = fb.to(device), cb.to(device), yb.to(device)
            optimizer.zero_grad()
            loss = criterion(model(od, fb, cb), yb)
            loss.backward()
            optimizer.step()

        val_rmse = float(np.sqrt(np.mean(
            (y_va - predict(model, omics_va, fp_va, codes_va, keys, device)) ** 2
        )))
        scheduler.step(val_rmse)

        if val_rmse < best_rmse - 1e-4:
            best_rmse, best_epoch, stale = val_rmse, epoch, 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            stale += 1

        if epoch % 10 == 0:
            print(f"  [epoch {epoch:>3}] val_rmse={val_rmse:.4f}  best={best_rmse:.4f} @ {best_epoch}")

        if stale >= PATIENCE:
            print(f"  [early stop] epoch {epoch}, best epoch {best_epoch}, val_rmse={best_rmse:.4f}")
            break

    model.load_state_dict(best_state)
    return model, best_epoch


def run_one(config_id, protocol, data, batched, n_fp, threshold, device, seed) -> dict:
    gathered, fingerprints, codes, y, groups = data
    cfg = CONFIGS[config_id]
    print("\n" + "#" * 82)
    print(f"# {config_id} | protocol={protocol} | seed={seed}")
    print(f"# {cfg['note']}")
    print("#" * 82)

    torch.manual_seed(seed)
    keys = list(gathered.keys())

    if protocol == "grouped":
        train_idx, val_idx, test_idx = grouped_split(groups)
    else:
        train_idx, val_idx, test_idx = random_pair_split(len(y))
    leak = leakage_report(groups, train_idx, test_idx)

    omics_tr = {k: v[train_idx] for k, v in gathered.items()}
    omics_va = {k: v[val_idx] for k, v in gathered.items()}
    omics_te = {k: v[test_idx] for k, v in gathered.items()}

    tensors = [torch.from_numpy(omics_tr[k]) for k in keys] + [
        torch.from_numpy(fingerprints[train_idx]),
        torch.from_numpy(codes[train_idx]),
        torch.from_numpy(y[train_idx]),
    ]
    loader = DataLoader(TensorDataset(*tensors), batch_size=BATCH_SIZE,
                        shuffle=True, drop_last=True)

    model = MoGraphDRPAligned(
        modalities=MODALITIES, drug_mode=cfg["drug_mode"], fusion=cfg["fusion"],
        head=cfg["head"], batched=batched, omics_in_dim=gathered[keys[0]].shape[1],
        dropout=DROPOUT,
    ).to(device)
    n_params = sum(p.numel() for p in model.parameters())

    t0 = time.perf_counter()
    model, best_epoch = train(
        model, loader, omics_va, fingerprints[val_idx], codes[val_idx], y[val_idx], keys, device
    )
    t_fit = time.perf_counter() - t0

    val = evaluate(y[val_idx], predict(model, omics_va, fingerprints[val_idx], codes[val_idx], keys, device), threshold)
    test = evaluate(y[test_idx], predict(model, omics_te, fingerprints[test_idx], codes[test_idx], keys, device), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])

    print_metric_block(f"{config_id} | {protocol}", val, test, floor)
    print(f"params={n_params:,}  best_epoch={best_epoch}  fit={t_fit:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {
        "config": config_id, "note": cfg["note"], "protocol": protocol, "seed": seed,
        "fusion": cfg["fusion"], "head": cfg["head"], "drug_rep": cfg["drug_mode"],
        "omics": "+".join(MODALITIES), "n_pairs": len(y), "mean_only_rmse": floor,
        **{f"leak_{k}": v for k, v in leak.items()},
        **{f"val_{k}": v for k, v in val.items()},
        **{f"test_{k}": v for k, v in test.items()},
        "params": n_params, "best_epoch": best_epoch, "fit_seconds": t_fit,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """`argv` is explicit so a notebook can call `main([...])` without argparse
    picking up the kernel's own `-f kernel.json` argument."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configs", nargs="+", choices=sorted(CONFIGS), default=sorted(CONFIGS))
    parser.add_argument("--protocols", nargs="+", choices=PROTOCOLS, default=list(PROTOCOLS))
    parser.add_argument("--seeds", nargs="+", type=int, default=[TORCH_SEED])
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    t_start = time.perf_counter()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    print(f"[sweep] {len(args.configs)} configs x {len(args.protocols)} protocols "
          f"x {len(args.seeds)} seeds = {len(args.configs) * len(args.protocols) * len(args.seeds)} runs")
    threshold = compute_shared_threshold()

    graphs = build_drug_graphs()
    assert_fingerprint_population_parity(graphs)
    gathered, fingerprints, codes, y, groups, batched, n_fp, _ = load_pairs(graphs)
    data = (gathered, fingerprints, codes, y, groups)

    results = [
        run_one(config_id, protocol, data, batched, n_fp, threshold, device, seed)
        for seed in args.seeds
        for protocol in args.protocols
        for config_id in args.configs
    ]

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"BENCHMARK ALIGNMENT COMPLETE — {len(df)} runs in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    pivot = df.pivot_table(index="config", columns="protocol", values="test_rmse", aggfunc="mean")
    print("\nTest RMSE by protocol:")
    print(pivot.to_string(float_format=lambda v: f"{v:.4f}"))
    print("\nReference (random split): MoGraphDRP 0.6622 | without XGBoost 0.9497 | BANDRP 0.9305")
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
