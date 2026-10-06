"""pair_attention_check.py — did the attention learn to find a mutated target?

Nothing in `PairAttention` tells it that a drug's target and a cell line's
mutated protein are the same protein. Both are looked up in one embedding
table, so the match can be learned, but only from the direct-hit training
pairs (989, 1.27% of the training fold). This script reads the trained
attention from the saved checkpoints and measures whether it was.

On every direct-hit pair, for each target that is itself mutated in the cell
line, it records the attention that target gives to:

    match    the key that is the same protein
    none     the "no mutation" token, a key of every cell line
    chance   1 / (number of keys), what uniform attention would give

and whether the matching key is the most-attended one (`top1`). The same is
computed for untrained modules (random initialisations), which is the level
the trained numbers have to beat. Weights are averaged over the 4 heads;
`best_head` is the head with the highest mean attention on the match.

It needs only the attention module's weights, not the rest of the model.

Run from the repository root, after the Stage A run:
    python src/final_model/pair_attention_check.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.experiment_utils import COL_CELL_LINE, grouped_split  # noqa: E402
from src.data.pair_features import build_pair_features, build_pair_sets, load_pair_rows  # noqa: E402
from src.final_model.run_ablation import checkpoints_dir  # noqa: E402
from src.models.pair_graph_drp import PairAttention  # noqa: E402

GATE_CSV = Path("src/final_model/results/pair_gate_results.csv")
CONFIG = "base+pair_attention"
N_UNTRAINED = 20
PREFIX = "pair_attention."


@torch.no_grad()
def match_attention(module: PairAttention, cells: torch.Tensor, drugs: torch.Tensor) -> dict:
    """Attention of each mutated target, over the given direct-hit pairs.

    One entry per (pair, target) where the target is among the pair's mutated
    proteins; a pair with two mutated targets contributes two.
    """
    module.eval()
    _, w = module(cells, drugs, return_weights=True)
    same = (w.query_index[:, :, None] == w.key_index[:, None, :]) \
        & w.query_mask[:, :, None] & w.key_mask[:, None, :]
    pair, target, key = same.nonzero(as_tuple=True)
    per_head = w.weights[pair, :, target, :]                       # (n, heads, keys)
    mean_head = per_head.mean(dim=1)                               # (n, keys)
    rows = torch.arange(len(pair))
    return {
        "match": mean_head[rows, key].numpy(),
        "none": mean_head[:, 0].numpy(),
        "chance": (1.0 / w.key_mask[pair].sum(dim=1)).numpy(),
        "top1": (mean_head.argmax(dim=1) == key).numpy(),
        "match_per_head": per_head[rows, :, key].numpy(),          # (n, heads)
        "pairs": int(len(torch.unique(pair))),
    }


def summarise(stats: dict) -> dict:
    return {
        "n_pairs": stats["pairs"], "n_matches": len(stats["match"]),
        "match": float(stats["match"].mean()), "none": float(stats["none"].mean()),
        "chance": float(stats["chance"].mean()), "top1": float(stats["top1"].mean()),
        "best_head": float(stats["match_per_head"].mean(axis=0).max()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, default=GATE_CSV,
                        help="results CSV of the run; checkpoints are read from beside it")
    args = parser.parse_args()

    runs = pd.read_csv(args.results)
    seeds = sorted(runs.loc[(runs["config"] == CONFIG) & (runs["protocol"] == "grouped"), "seed"])
    if not seeds:
        raise SystemExit(f"no grouped-split '{CONFIG}' rows in {args.results}")

    y_used = load_pair_rows()
    sets = build_pair_sets(y_used)
    direct = build_pair_features(y_used).subsets["direct_hit"]
    folds = dict(zip(("train", "val", "test"), grouped_split(y_used[COL_CELL_LINE].to_numpy())))
    cells, drugs = torch.from_numpy(sets.cell_codes), torch.from_numpy(sets.drug_codes)

    rows = []
    for fold, idx in folds.items():
        hit = idx[direct[idx]]
        for seed in seeds:
            module = PairAttention(sets)
            saved = torch.load(checkpoints_dir(args.results) / f"{CONFIG}__grouped__seed{seed}.pt")
            module.load_state_dict({k[len(PREFIX):]: v for k, v in saved.items()
                                    if k.startswith(PREFIX)})
            rows.append({"fold": fold, "model": "trained", "seed": seed,
                         **summarise(match_attention(module, cells[hit], drugs[hit]))})
        for seed in range(N_UNTRAINED):
            torch.manual_seed(seed)
            rows.append({"fold": fold, "model": "untrained", "seed": seed,
                         **summarise(match_attention(PairAttention(sets), cells[hit], drugs[hit]))})

    table = pd.DataFrame(rows)
    out = args.results.with_name(args.results.stem + "_attention_check.csv")
    table.to_csv(out, index=False)

    print("\nAttention of a mutated target, direct-hit pairs "
          "(mean over heads; mean +/- std over seeds or initialisations)")
    print(f"{'fold':<7}{'model':<11}{'n':>4}{'pairs':>7}{'match':>18}{'none':>18}"
          f"{'chance':>9}{'top1':>18}{'best head':>18}")
    for (fold, model), block in table.groupby(["fold", "model"], sort=False):
        cell = lambda c: f"{block[c].mean():.3f} +/- {block[c].std():.3f}"  # noqa: E731
        print(f"{fold:<7}{model:<11}{len(block):>4}{block['n_pairs'].iloc[0]:>7}"
              f"{cell('match'):>18}{cell('none'):>18}{block['chance'].mean():>9.3f}"
              f"{cell('top1'):>18}{cell('best_head'):>18}")
    trained = table[table["model"] == "trained"]
    print("\nPer seed, trained:")
    print(trained.pivot(index="seed", columns="fold", values="match")[list(folds)]
          .round(3).to_string())
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    main()
