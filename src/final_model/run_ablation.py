"""run_ablation.py — the additive ablation ladder on the MoGraphDRP-style model.

The ladder is **Base + one component at a time**:

    base                    the simplest pipeline that functions with none of
                            the improvements: per-omics branches -> concat,
                            Morgan fingerprint encoder, concat -> MLP head,
                            raw ln(IC50) targets, GE + Mut_CNV
    base+<component>        exactly one switch changed, so the delta against
                            `base` measures that component alone
    aligned                 base+mol_graph+bilinear, the MoGraphDRP reproduction
    full                    base + every component

Components are declared once in `COMPONENTS`; adding a new one there adds its
`base+<name>` rung to the default ladder and to `full`. Any combination can be
requested by name, e.g. `--configs base+bilinear+proteomics`.

Every configuration can be evaluated under **both** the benchmark's random
pair split and this project's cell-line-grouped split. The grouped split is the
headline; the random split exists only to sanity-check the reproduction against
the published figures below.

Reference points, all under the random split:

| MoGraphDRP published (with XGBoost)      | 0.6622 |
| MoGraphDRP without XGBoost (Table 6)     | 0.9497 |
| BANDRP, best prior work (their Table 2)  | 0.9305 |

Results are appended one row per finished run, so an interrupted sweep loses
nothing and re-running skips what is already done (`--rerun` to repeat).

Run from the repository root:
    python "src/final_model/run_ablation.py" --list
    python "src/final_model/run_ablation.py" --protocols grouped --seeds 42 43 44
    python "src/final_model/run_ablation.py" --configs base base+bilinear
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
    MUT_CNV_KEY,
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
    per_drug_mean_floor,
    print_metric_block,
    random_pair_split,
)
from src.data.target_scaling import PerDrugTargetScaler  # noqa: E402
from src.final_model.model import (  # noqa: E402
    BATCH_SIZE,
    DROPOUT,
    LR,
    MAX_EPOCHS,
    MoGraphDRPAligned,
)

RESULTS_CSV = Path("src/final_model/results/ablation_results.csv")

ALL_MODALITIES = (GE_KEY, MUT_CNV_KEY, PROTEOMICS_KEY)
TORCH_SEED = 42
WEIGHT_DECAY = 1e-5
PATIENCE = 15
LR_PATIENCE = 5

# The simplest pipeline that functions with none of the improvements. Two omics
# are the minimum for "fusion" to mean anything, and GE + Mut_CNV is the pair
# closest to the benchmark's own inputs, which leaves proteomics as a component.
BASE: dict = dict(
    fusion="concat",
    head="mlp",
    drug_mode="fingerprint",
    modalities=(GE_KEY, MUT_CNV_KEY),
    standardize_targets=False,
)

# One entry per improvement: the single BASE field it overrides, and what it is.
COMPONENTS: dict[str, tuple[dict, str]] = {
    "cross_attention": (dict(fusion="cross_attention"),
                        "cross-attention omics fusion instead of concatenation"),
    "mol_graph": (dict(drug_mode="both"),
                  "molecular-graph GCN alongside the fingerprint (their sec 2.2.3)"),
    "bilinear": (dict(head="bilinear"),
                 "multi-head bilinear cell x drug head instead of concat (their sec 2.3)"),
    "proteomics": (dict(modalities=ALL_MODALITIES),
                   "proteomics as a third omics branch"),
    "std_targets": (dict(standardize_targets=True),
                    "per-drug standardized ln(IC50) targets"),
}
# Each component must own a different field, otherwise two of them could not be
# combined and `base+a` vs `base+b` would not be independent switches.
assert len({k for override, _ in COMPONENTS.values() for k in override}) == len(COMPONENTS)

ALIASES: dict[str, str] = {
    "aligned": "base+mol_graph+bilinear",
    "full": "+".join(["base", *COMPONENTS]),
}
DEFAULT_LADDER = ["base", *[f"base+{name}" for name in COMPONENTS], "aligned", "full"]
PROTOCOLS = ("grouped", "random")


def resolve_config(config_id: str) -> dict:
    """`base`, `base+a+b`, or an alias -> the full set of model/training switches."""
    parts = ALIASES.get(config_id, config_id).split("+")
    unknown = [p for p in parts[1:] if p not in COMPONENTS]
    if parts[0] != "base" or unknown:
        raise ValueError(
            f"unknown config {config_id!r}: expected 'base', 'base+<component>[+...]' "
            f"or one of {sorted(ALIASES)}; components are {list(COMPONENTS)}"
        )
    cfg = dict(BASE)
    for name in parts[1:]:
        cfg.update(COMPONENTS[name][0])
    cfg["components"] = "+".join(parts[1:]) or "none"
    return cfg


def print_ladder(config_ids: list[str]) -> None:
    print(f"{'config':<34}{'fusion':<17}{'head':<10}{'drug':<13}{'std_y':<7}omics")
    print("-" * 104)
    for config_id in config_ids:
        cfg = resolve_config(config_id)
        print(f"{config_id:<34}{cfg['fusion']:<17}{cfg['head']:<10}{cfg['drug_mode']:<13}"
              f"{str(cfg['standardize_targets']):<7}{'+'.join(cfg['modalities'])}")
    print("\ncomponents:")
    for name, (_, description) in COMPONENTS.items():
        print(f"  {name:<17}{description}")


def load_pairs(graphs):
    """Fingerprints and molecular-graph codes on one identical row set.

    Both drug representations are needed simultaneously for `drug_mode="both"`,
    and they must describe the same pairs in the same order. Rather than
    trusting two independent builders to agree, the fingerprint path defines
    the rows and the graph codes are derived from *its* `y_used`.

    All three modalities are loaded once; each configuration then takes the
    subset it needs. The omics files share one cell-line index, so the pair
    rows -- and therefore the split -- are identical for every configuration.
    """
    omics, cell_ids = load_omics_subset(ALL_MODALITIES)
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
def predict(model, omics, fingerprints, codes, idx, keys, device) -> np.ndarray:
    """Raw model output for the rows in `idx` (z-scores when targets are standardized)."""
    model.eval()
    preds = []
    for start in range(0, len(idx), BATCH_SIZE):
        rows = idx[start:start + BATCH_SIZE]
        ob = {k: torch.from_numpy(omics[k][rows]).to(device) for k in keys}
        fb = torch.from_numpy(fingerprints[rows]).to(device)
        cb = torch.from_numpy(codes[rows]).to(device)
        preds.append(model(ob, fb, cb).cpu().numpy())
    return np.concatenate(preds)


def train(model, loader, predict_ln, val_idx, y_va, keys, device, max_epochs) -> tuple[nn.Module, int]:
    """`predict_ln` returns predictions in ln(IC50) units, so the stopping
    criterion is identical whether or not the targets were standardized."""
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=LR_PATIENCE
    )
    best_state = copy.deepcopy(model.state_dict())
    best_rmse, best_epoch, stale = float("inf"), 0, 0

    for epoch in range(1, max_epochs + 1):
        model.train()
        for *omics_batches, fb, cb, yb in loader:
            od = {k: t.to(device) for k, t in zip(keys, omics_batches)}
            fb, cb, yb = fb.to(device), cb.to(device), yb.to(device)
            optimizer.zero_grad()
            loss = criterion(model(od, fb, cb), yb)
            loss.backward()
            optimizer.step()

        val_rmse = float(np.sqrt(np.mean((y_va - predict_ln(val_idx)) ** 2)))
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


def run_one(config_id, protocol, data, batched, threshold, device, seed, max_epochs) -> dict:
    gathered, fingerprints, codes, y, groups = data
    cfg = resolve_config(config_id)
    keys = list(cfg["modalities"])
    print("\n" + "#" * 82)
    print(f"# {config_id} | protocol={protocol} | seed={seed}")
    print(f"# components: {cfg['components']} | omics: {'+'.join(keys)}")
    print("#" * 82)

    torch.manual_seed(seed)

    if protocol == "grouped":
        train_idx, val_idx, test_idx = grouped_split(groups)
    else:
        train_idx, val_idx, test_idx = random_pair_split(len(y))
    leak = leakage_report(groups, train_idx, test_idx)

    # Drug codes are the grouping key for per-drug statistics. Fitted on
    # training rows only; `y` itself is never modified, so every metric below
    # is computed in ln(IC50) units either way.
    if cfg["standardize_targets"]:
        scaler = PerDrugTargetScaler().fit(y[train_idx], codes[train_idx])
        y_fit = scaler.transform(y, codes)
    else:
        scaler, y_fit = None, y

    tensors = [torch.from_numpy(gathered[k][train_idx]) for k in keys] + [
        torch.from_numpy(fingerprints[train_idx]),
        torch.from_numpy(codes[train_idx]),
        torch.from_numpy(y_fit[train_idx]),
    ]
    loader = DataLoader(TensorDataset(*tensors), batch_size=BATCH_SIZE,
                        shuffle=True, drop_last=True)

    model = MoGraphDRPAligned(
        modalities=keys, drug_mode=cfg["drug_mode"], fusion=cfg["fusion"],
        head=cfg["head"], batched=batched, omics_in_dim=gathered[keys[0]].shape[1],
        dropout=DROPOUT,
    ).to(device)
    n_params = sum(p.numel() for p in model.parameters())

    def predict_ln(idx: np.ndarray) -> np.ndarray:
        raw = predict(model, gathered, fingerprints, codes, idx, keys, device)
        return raw if scaler is None else scaler.inverse_transform(raw, codes[idx])

    t0 = time.perf_counter()
    model, best_epoch = train(model, loader, predict_ln, val_idx, y[val_idx], keys, device, max_epochs)
    t_fit = time.perf_counter() - t0

    val = evaluate(y[val_idx], predict_ln(val_idx), threshold)
    test = evaluate(y[test_idx], predict_ln(test_idx), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])
    drug_floor = per_drug_mean_floor(y[train_idx], codes[train_idx], y[test_idx], codes[test_idx])

    print_metric_block(f"{config_id} | {protocol}", val, test, floor, drug_floor)
    print(f"params={n_params:,}  best_epoch={best_epoch}  fit={t_fit:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {
        "config": config_id, "components": cfg["components"], "protocol": protocol, "seed": seed,
        "fusion": cfg["fusion"], "head": cfg["head"], "drug_rep": cfg["drug_mode"],
        "omics": "+".join(keys), "standardize_targets": cfg["standardize_targets"],
        "n_pairs": len(y), "mean_only_rmse": floor, "per_drug_mean_rmse": drug_floor,
        "gain_over_drug_lookup": drug_floor - test["rmse"],
        **{f"leak_{k}": v for k, v in leak.items()},
        **{f"val_{k}": v for k, v in val.items()},
        **{f"test_{k}": v for k, v in test.items()},
        "params": n_params, "best_epoch": best_epoch, "max_epochs": max_epochs,
        "fit_seconds": t_fit,
    }


def print_summary(df: pd.DataFrame) -> None:
    """Mean +/- std over seeds per config, with the delta against `base`."""
    for protocol, block in df.groupby("protocol"):
        stats = block.groupby("config")["test_rmse"].agg(["mean", "std", "count"])
        base = stats.loc["base", "mean"] if "base" in stats.index else float("nan")
        print(f"\nTest RMSE, {protocol} split:")
        print(f"  {'config':<34}{'mean':>9}{'std':>9}{'n':>4}{'vs base':>10}")
        for config_id, row in stats.sort_values("mean").iterrows():
            print(f"  {config_id:<34}{row['mean']:>9.4f}{row['std']:>9.4f}{int(row['count']):>4}"
                  f"{row['mean'] - base:>+10.4f}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """`argv` is explicit so a notebook can call `main([...])` without argparse
    picking up the kernel's own `-f kernel.json` argument."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--configs", nargs="+", default=DEFAULT_LADDER,
                        help="'base', 'base+<component>[+...]', 'aligned' or 'full'")
    parser.add_argument("--protocols", nargs="+", choices=PROTOCOLS, default=list(PROTOCOLS))
    parser.add_argument("--seeds", nargs="+", type=int, default=[TORCH_SEED])
    parser.add_argument("--max-epochs", type=int, default=MAX_EPOCHS)
    parser.add_argument("--out", type=Path, default=RESULTS_CSV,
                        help="results CSV to append to (use a scratch path for smoke tests)")
    parser.add_argument("--rerun", action="store_true",
                        help="repeat runs already present in the results CSV")
    parser.add_argument("--list", action="store_true", help="print the ladder and exit")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    for config_id in args.configs:
        resolve_config(config_id)  # fail on a typo before loading any data
    if args.list:
        print_ladder(args.configs)
        return
    t_start = time.perf_counter()

    done: set[tuple] = set()
    if args.out.exists() and not args.rerun:
        previous = pd.read_csv(args.out)
        done = set(zip(previous["config"], previous["protocol"], previous["seed"]))
    todo = [
        (config_id, protocol, seed)
        for seed in args.seeds
        for protocol in args.protocols
        for config_id in args.configs
        if (config_id, protocol, seed) not in done
    ]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    print(f"[sweep] {len(args.configs)} configs x {len(args.protocols)} protocols "
          f"x {len(args.seeds)} seeds: {len(todo)} to run, "
          f"{len(args.configs) * len(args.protocols) * len(args.seeds) - len(todo)} already in {args.out}")
    if not todo:
        print_summary(pd.read_csv(args.out))
        return
    threshold = compute_shared_threshold()

    graphs = build_drug_graphs()
    assert_fingerprint_population_parity(graphs)
    gathered, fingerprints, codes, y, groups, batched, _, _ = load_pairs(graphs)
    data = (gathered, fingerprints, codes, y, groups)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    for config_id, protocol, seed in todo:
        row = run_one(config_id, protocol, data, batched, threshold, device, seed, args.max_epochs)
        pd.DataFrame([row]).to_csv(args.out, mode="a", header=not args.out.exists(), index=False)

    print("\n" + "=" * 100)
    print(f"ABLATION LADDER: {len(todo)} runs in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    print_summary(pd.read_csv(args.out))
    print("\nReference (random split): MoGraphDRP 0.6622 | without XGBoost 0.9497 | BANDRP 0.9305")
    print(f"\nSaved -> {args.out}")


if __name__ == "__main__":
    main()
