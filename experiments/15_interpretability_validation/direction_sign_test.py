"""Check 1 -- does the trained model get the *direction* of BRAF right?

The graph tells the model that a cell line carries a BRAF driver mutation and
that PLX-4720 / SB590885 / Dabrafenib target BRAF. It does not say whether that
mutation makes those drugs work better or worse; only the IC50 labels carry
that. Known biology: BRAF-mutant lines are *more* sensitive to BRAF inhibitors,
so a model that learned it should attribute a negative (sensitising)
contribution to the BRAF node on these pairs.

Untrained controls see the identical edges, so whatever sign they produce is
what wiring alone gives. Signs within one model are strongly correlated across
pairs (same weights, same BRAF embedding), so the per-pair binomial test is
reported for completeness but the honest comparison is model-level: where the
trained model's hit rate falls among the 50 controls' hit rates.

Run from the repository root:
    python experiments/15_interpretability_validation/direction_sign_test.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import binomtest, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import random_baseline as rb  # noqa: E402

RESULTS_CSV = rb.OUT_DIR / "direction_sign_test_results.csv"
SUMMARY_CSV = rb.OUT_DIR / "direction_sign_test_summary.csv"


def model_stats(braf_scores: np.ndarray, ln_ic50: np.ndarray) -> dict:
    """Hit rate (fraction sensitising) plus the Spearman between BRAF score and measured ln(IC50).

    A positive rho means pairs the lab measured as more sensitive also got a
    more negative BRAF attribution -- the scale-free version of "right sign".
    """
    rho = spearmanr(braf_scores, ln_ic50).statistic if np.std(braf_scores) > 0 else np.nan
    return {
        "hit_rate": float((braf_scores < 0).mean()),
        "mean_braf_score": float(braf_scores.mean()),
        "spearman_vs_ln_ic50": float(rho),
    }


def summarise(label: str, trained: dict, randoms: list[dict], n_pairs: int) -> dict:
    rand_hits = np.array([r["hit_rate"] for r in randoms])
    rand_rho = np.array([r["spearman_vs_ln_ic50"] for r in randoms])
    n_hits = round(trained["hit_rate"] * n_pairs)
    return {
        "subset": label,
        "n_pairs": n_pairs,
        "trained_hit_rate": trained["hit_rate"],
        "trained_binom_p_vs_0.5": binomtest(n_hits, n_pairs, 0.5).pvalue,
        "random_hit_rate_mean": rand_hits.mean(),
        "random_hit_rate_std": rand_hits.std(ddof=1),
        "random_models_all_sensitising": int((rand_hits == 1.0).sum()),
        "random_models_all_resistance": int((rand_hits == 0.0).sum()),
        # One-sided empirical p: share of controls at least as "right" as the trained model.
        "empirical_p_hit_rate": float((rand_hits >= trained["hit_rate"]).mean()),
        "trained_spearman": trained["spearman_vs_ln_ic50"],
        "random_spearman_mean": float(np.nanmean(rand_rho)),
        "random_spearman_std": float(np.nanstd(rand_rho, ddof=1)),
        "empirical_p_spearman": float((rand_rho >= trained["spearman_vs_ln_ic50"]).mean()),
    }


def main() -> int:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    graph_maps = rb.load_graph_and_maps(device)
    braf_idx = rb.protein_index_for_symbol(graph_maps, rb.BRAF_SYMBOL)
    pairs = rb.braf_case_study_pairs(graph_maps)
    cells, drugs = pairs["cell_idx"].tolist(), pairs["drug_idx"].tolist()
    print(f"[setup] device={device}  pairs={len(pairs)}  BRAF node={braf_idx}  controls={rb.N_RANDOM_MODELS}")

    trained = rb.load_trained_model(graph_maps, device)
    trained_braf = rb.score_pairs(trained, graph_maps, cells, drugs)[:, braf_idx].cpu().numpy()

    random_braf = np.empty((rb.N_RANDOM_MODELS, len(pairs)))
    start = time.time()
    for k, model in enumerate(rb.build_random_models(graph_maps, device)):
        random_braf[k] = rb.score_pairs(model, graph_maps, cells, drugs)[:, braf_idx].cpu().numpy()
        print(f"[control {k + 1:>2}/{rb.N_RANDOM_MODELS}] hit_rate={(random_braf[k] < 0).mean():.3f}  "
              f"({time.time() - start:.0f}s)")

    out = pairs[["sanger_model_id", "cell_line_name", "drug_id", "drug_name", "ln_ic50"]].copy()
    out["trained_braf_score"] = trained_braf
    out["trained_direction"] = np.where(trained_braf < 0, "sensitising", "resistance")
    out["random_braf_score_mean"] = random_braf.mean(axis=0)
    out["random_models_sensitising"] = (random_braf < 0).sum(axis=0)
    out.to_csv(RESULTS_CSV, index=False)

    ln_ic50 = pairs["ln_ic50"].to_numpy()
    subsets = [("all", np.ones(len(pairs), dtype=bool))] + [
        (name, (pairs["drug_id"] == drug_id).to_numpy()) for drug_id, name in rb.BRAF_INHIBITORS.items()
    ]
    rows = []
    for label, mask in subsets:
        t = model_stats(trained_braf[mask], ln_ic50[mask])
        r = [model_stats(random_braf[k, mask], ln_ic50[mask]) for k in range(rb.N_RANDOM_MODELS)]
        rows.append(summarise(label, t, r, int(mask.sum())))
    summary = pd.DataFrame(rows)
    summary.to_csv(SUMMARY_CSV, index=False)

    print("\n" + "=" * 100)
    print("Check 1: BRAF attribution direction, BRAF-mutant cell lines x BRAF inhibitors")
    print("=" * 100)
    with pd.option_context("display.width", 200, "display.max_columns", None, "display.float_format", "{:.4f}".format):
        print(summary.T.to_string(header=False))
    print(f"\nSaved {RESULTS_CSV.name}, {SUMMARY_CSV.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
