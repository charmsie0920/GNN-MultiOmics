"""pair_feature_check.py — does the trained final model use the direct-hit signal?

The final model's only pair-specific input is five hand-built features
(`src/data/pair_features.py`). This script measures what the trained model
does with the one that carries signal: "a target of the drug is itself
mutated in the cell line" (a direct hit, 1.26% of pairs).

For a direct-hit pair it predicts twice: with the real features, and with the
hit removed. "Removed" describes the same pair as if its nearest mutation
were one hop from a target (`direct` 1 -> 0, `proximity` 1.0 -> 0.5, the
other three features unchanged), an input that 54,844 real pairs have. The
difference, real minus removed, is the model's *effect* of the hit in
ln(IC50) units: negative means the hit makes the model predict the pair more
sensitive. The omics inputs are untouched, so the model still sees the
mutation through Mut_CNV; only the pair features change.

Three tests, with pass rules written in `CLAUDE.md` (Phase 6) before the
weights existed:

    use        RMSE on direct-hit pairs of held-out cell lines, real features
               against hit removed
    drugs      across drugs, does the model's effect track the measured
               difference between hit and non-hit held-out cell lines?
    braf       across BRAF-mutant cell lines on the three BRAF inhibitors,
               is the effect most sensitising where the measured ln(IC50)
               is lowest? (the test of experiment 15)

The same numbers are computed for untrained models of the same architecture
(random initialisations), which is the level the trained ones have to beat.

It is a sensitivity test on hand-built inputs, not a gene-level attribution.

Run from the repository root, after the weights exist:
    python src/final_model/run_ablation.py --configs final --protocols grouped \
        --seeds 42 43 44 45 46 --save-weights \
        --out src/final_model/results/final_weights_results.csv
    python src/final_model/pair_feature_check.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.experiment_utils import COL_CELL_LINE, COL_DRUG, grouped_split  # noqa: E402
from src.data.pair_features import (  # noqa: E402
    FEATURE_NAMES,
    GRAPH_PATH,
    _load_pair_graph,
    build_pair_features,
)
from src.data.target_scaling import PerDrugTargetScaler  # noqa: E402
from src.final_model.run_ablation import (  # noqa: E402
    DROPOUT,
    build_drug_graphs,
    checkpoints_dir,
    load_pairs,
    predict,
    resolve_config,
)
from src.models.pair_graph_drp import PairGraphDRP  # noqa: E402

WEIGHTS_CSV = Path("src/final_model/results/final_weights_results.csv")
SYMBOLS_CSV = Path("data/processed/protein_symbol_map.csv")
CONFIG = "final"
N_UNTRAINED = 50
MIN_HITS_PER_DRUG = 3
HOP_1_PROXIMITY = 0.5
CASE_GENE = "BRAF"


def rmse(error: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(error))))


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    return float(stats.spearmanr(a, b).statistic)


def remove_hit(features: np.ndarray) -> np.ndarray:
    """The feature matrix with every direct hit described as a 1-hop pair."""
    direct, proximity = FEATURE_NAMES.index("direct"), FEATURE_NAMES.index("proximity")
    removed = features.copy()
    hit = removed[:, direct] == 1
    removed[hit, direct] = 0.0
    removed[hit, proximity] = HOP_1_PROXIMITY
    return removed


def case_pairs(y_used: pd.DataFrame, gene: str) -> np.ndarray:
    """Row mask: the drug targets `gene` and the cell line carries a mutation in it."""
    graph = _load_pair_graph(y_used, GRAPH_PATH)
    symbols = pd.read_csv(SYMBOLS_CSV).set_index("gene_symbol")["protein_id"]
    node = [str(p) for p in graph.data["protein"].node_ids].index(symbols[gene])
    drug_hit = np.array([node in graph.targets_of.get(d, ()) for d in graph.drug_levels])
    cell_hit = np.array([node in graph.mutated_of.get(c, ()) for c in graph.cell_levels])
    return drug_hit[graph.drug_codes] & cell_hit[graph.cell_codes]


def per_drug(effect, z, std, direct, idx, drugs) -> pd.DataFrame:
    """Per drug with enough direct hits among the rows `idx`: the model's mean
    effect and the measured hit - non-hit difference, both in units of the
    drug's training std."""
    table = pd.DataFrame({"drug": drugs[idx], "hit": direct[idx], "z": z[idx],
                          "effect_z": effect[idx] / std[idx]})
    counts = table.groupby("drug")["hit"].sum()
    table = table[table["drug"].isin(counts[counts >= MIN_HITS_PER_DRUG].index)]
    hits = table[table["hit"]].groupby("drug")
    return pd.DataFrame({
        "n_hits": hits.size(),
        "measured": hits["z"].mean() - table[~table["hit"]].groupby("drug")["z"].mean(),
        "model": hits["effect_z"].mean(),
    })


def measure(pred: dict, y, z, std, direct, case, folds, drugs) -> list[dict]:
    """The three tests for one model. `pred` holds ln(IC50) predictions on all
    rows of interest, with the real features and with the hit removed."""
    effect = pred["real"] - pred["removed"]
    held_out = np.concatenate([folds["val"], folds["test"]])
    rows = []

    for fold in ("val", "test"):
        hit = folds[fold][direct[folds[fold]]]
        rows.append({
            "test": "use", "subset": fold, "n": len(hit),
            "rmse_real": rmse(pred["real"][hit] - y[hit]),
            "rmse_removed": rmse(pred["removed"][hit] - y[hit]),
            "bias_real": float(np.mean(pred["real"][hit] - y[hit])),
            "bias_removed": float(np.mean(pred["removed"][hit] - y[hit])),
            "effect": float(effect[hit].mean()),
        })

    # Across drugs: the model's effect against the measured hit - non-hit
    # difference, both in units of the drug's training std.
    for subset, idx in (("held_out", held_out), ("test", folds["test"])):
        table = per_drug(effect, z, std, direct, idx, drugs)
        rows.append({"test": "drugs", "subset": subset, "n": len(table),
                     "spearman": spearman(table["model"], table["measured"])})

    # The case study: across pairs, effect against measured ln(IC50).
    in_fold = {"all": np.concatenate(list(folds.values())), "train": folds["train"],
               "held_out": held_out}
    for subset, idx in in_fold.items():
        hit = idx[case[idx]]
        rows.append({
            "test": "braf", "subset": subset, "n": len(hit),
            "spearman": spearman(effect[hit], y[hit]),
            # The same within drug, so a difference between the three drugs'
            # potency cannot produce it.
            "spearman_within_drug": spearman(effect[hit] / std[hit], z[hit]),
            "effect": float(effect[hit].mean()),
            "share_sensitising": float(np.mean(effect[hit] < 0)),
            # Not part of a rule. Whether the hit makes these pairs better
            # predicted, and how far the effect only mirrors what the model
            # predicts without it (the caveat of experiment 15).
            "rmse_real": rmse(pred["real"][hit] - y[hit]),
            "rmse_removed": rmse(pred["removed"][hit] - y[hit]),
            "bias_real": float(np.mean(pred["real"][hit] - y[hit])),
            "bias_removed": float(np.mean(pred["removed"][hit] - y[hit])),
            "spearman_effect_vs_removed": spearman(effect[hit], pred["removed"][hit]),
            "spearman_removed_vs_measured": spearman(pred["removed"][hit], y[hit]),
            "spearman_real_vs_measured": spearman(pred["real"][hit], y[hit]),
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, default=WEIGHTS_CSV,
                        help="results CSV of the run; checkpoints are read from beside it")
    parser.add_argument("--untrained", type=int, default=N_UNTRAINED,
                        help="number of untrained models scored as the null")
    args = parser.parse_args()

    runs = pd.read_csv(args.results)
    seeds = sorted(runs.loc[(runs["config"] == CONFIG) & (runs["protocol"] == "grouped"), "seed"])
    if not seeds:
        raise SystemExit(f"no grouped-split '{CONFIG}' rows in {args.results}")

    cfg = resolve_config(CONFIG)
    keys = list(cfg["modalities"])
    gathered, fingerprints, codes, y, groups, batched, _, y_used = load_pairs(build_drug_graphs())
    pair = build_pair_features(y_used)
    direct = pair.subsets["direct_hit"]
    case = case_pairs(y_used, CASE_GENE)
    assert not (case & ~direct).any(), "a case pair that is not a direct hit"
    folds = dict(zip(("train", "val", "test"), grouped_split(groups)))
    drugs = y_used[COL_DRUG].to_numpy()

    # `final` trains on per-drug z-scores; the same training-fold statistics
    # turn every model's output back into ln(IC50).
    scaler = PerDrugTargetScaler().fit(y[folds["train"]], codes[folds["train"]])
    z = scaler.transform(y, codes)
    std = scaler.params_for(codes)[1]

    inputs = {"real": pair.matrix, "removed": remove_hit(pair.matrix)}
    # Predictions are only needed where the two inputs differ.
    hit_rows = np.flatnonzero(direct)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    def build() -> PairGraphDRP:
        return PairGraphDRP(
            modalities=keys, drug_mode=cfg["drug_mode"], fusion=cfg["fusion"], head=cfg["head"],
            batched=batched, omics_in_dim=gathered[keys[0]].shape[1], dropout=DROPOUT,
            pair_dim=pair.matrix.shape[1],
        ).to(device)

    def score(model: PairGraphDRP) -> dict:
        pred = {}
        for name, matrix in inputs.items():
            raw = predict(model, gathered, fingerprints, codes, hit_rows, keys, device, matrix)
            full = np.full(len(y), np.nan, dtype=np.float64)
            full[hit_rows] = scaler.inverse_transform(raw, codes[hit_rows])
            pred[name] = full
        return pred

    held_out = np.concatenate([folds["val"], folds["test"]])
    rows, by_drug = [], []
    for seed in seeds:
        model = build()
        model.load_state_dict(torch.load(checkpoints_dir(args.results)
                                         / f"{CONFIG}__grouped__seed{seed}.pt"))
        pred = score(model)
        rows += [{"model": "trained", "seed": seed, **r}
                 for r in measure(pred, y, z, std, direct, case, folds, drugs)]
        by_drug.append(per_drug(pred["real"] - pred["removed"], z, std, direct, held_out, drugs))
    for seed in range(args.untrained):
        torch.manual_seed(seed)
        rows += [{"model": "untrained", "seed": seed, **r}
                 for r in measure(score(build()), y, z, std, direct, case, folds, drugs)]

    table = pd.DataFrame(rows)
    out = args.results.with_name(args.results.stem + "_pair_feature_check.csv")
    table.to_csv(out, index=False)
    # Per drug on held-out cell lines, the model column averaged over seeds.
    drug_table = by_drug[0][["n_hits", "measured"]].assign(
        model=np.mean([t["model"] for t in by_drug], axis=0),
        model_std=np.std([t["model"] for t in by_drug], axis=0, ddof=1),
    ).sort_values("measured")
    drug_out = args.results.with_name(args.results.stem + "_pair_feature_check_drugs.csv")
    drug_table.to_csv(drug_out)
    trained, untrained = table[table["model"] == "trained"], table[table["model"] == "untrained"]

    print("\n1. Use: direct-hit pairs of held-out cell lines, trained models (mean +/- std over seeds)")
    print(f"{'fold':<6}{'pairs':>6}{'RMSE real':>18}{'RMSE removed':>18}{'seeds better':>14}"
          f"{'bias real':>11}{'bias removed':>14}{'effect':>9}")
    for fold, block in trained[trained["test"] == "use"].groupby("subset", sort=False):
        print(f"{fold:<6}{block['n'].iloc[0]:>6}"
              f"{block['rmse_real'].mean():>10.4f} +/- {block['rmse_real'].std():.3f}"
              f"{block['rmse_removed'].mean():>10.4f} +/- {block['rmse_removed'].std():.3f}"
              f"{int((block['rmse_real'] < block['rmse_removed']).sum()):>10} / {len(block)}"
              f"{block['bias_real'].mean():>+11.3f}{block['bias_removed'].mean():>+14.3f}"
              f"{block['effect'].mean():>+9.3f}")

    print(f"\n2 and 3. Spearman correlations, trained seeds against {args.untrained} untrained models")
    print(f"{'test':<7}{'subset':<10}{'column':<22}{'n':>5}{'trained':>18}{'seeds > 0':>11}"
          f"{'untrained':>20}{'controls reaching':>19}")
    for (test, subset), block in trained[trained["test"] != "use"].groupby(["test", "subset"], sort=False):
        null = untrained[(untrained["test"] == test) & (untrained["subset"] == subset)]
        for column in ("spearman", "spearman_within_drug"):
            if block[column].isna().all():
                continue
            mean = block[column].mean()
            print(f"{test:<7}{subset:<10}{column:<22}{int(block['n'].iloc[0]):>5}"
                  f"{mean:>+9.3f} +/- {block[column].std():.3f}"
                  f"{int((block[column] > 0).sum()):>7} / {len(block)}"
                  f"{null[column].mean():>+11.3f} +/- {null[column].std():.3f}"
                  f"{int((null[column] >= mean).sum()):>13} / {len(null)}")

    print("\nBRAF case study: mean effect and share of pairs where the hit is sensitising")
    for subset, block in trained[trained["test"] == "braf"].groupby("subset", sort=False):
        null = untrained[(untrained["test"] == "braf") & (untrained["subset"] == subset)]
        print(f"  {subset:<9} n={int(block['n'].iloc[0]):<4} trained effect {block['effect'].mean():+.3f}, "
              f"share {block['share_sensitising'].mean():.3f} | untrained effect "
              f"{null['effect'].mean():+.3f} +/- {null['effect'].std():.3f}, "
              f"share {null['share_sensitising'].mean():.3f} +/- {null['share_sensitising'].std():.3f}")
    block = trained[(trained["test"] == "braf") & (trained["subset"] == "held_out")]
    print(f"  held-out pairs, trained: RMSE {block['rmse_real'].mean():.3f} with the hit, "
          f"{block['rmse_removed'].mean():.3f} removed; Spearman of the prediction with measured "
          f"ln(IC50) {block['spearman_real_vs_measured'].mean():+.3f} with, "
          f"{block['spearman_removed_vs_measured'].mean():+.3f} removed; "
          f"of the effect with the prediction without it "
          f"{block['spearman_effect_vs_removed'].mean():+.3f}")
    print(f"\nSaved -> {out}\n         {drug_out}")


if __name__ == "__main__":
    main()
