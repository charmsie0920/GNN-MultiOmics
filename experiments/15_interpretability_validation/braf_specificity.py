"""Follow-up to Check 1: is the BRAF-score / measured-sensitivity correlation specific and out-of-sample?

`direction_sign_test.py` found the trained model's BRAF attribution tracks
measured ln(IC50) across the 157 BRAF-mutant x BRAF-inhibitor pairs, and no
untrained control does. Two ways that could still be uninteresting:

1. Every protein's gradient x embedding scales with the prediction, so *any*
   protein would correlate with ln(IC50). Tested by ranking BRAF's correlation
   against every other protein with a non-constant score.
2. The pairs were in the training set, so the correlation is memorised.
   Tested by recomputing it on cell lines held out of training, using the
   split the checkpoint was trained on (data/raw/aligned_ic50_pairs.csv).

Run from the repository root:
    python experiments/15_interpretability_validation/braf_specificity.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import random_baseline as rb  # noqa: E402

SPLIT_CSV = rb.REPO_ROOT / "data" / "raw" / "aligned_ic50_pairs.csv"
RESULTS_CSV = rb.OUT_DIR / "braf_specificity_results.csv"
N_PERMUTATIONS = 10_000


def permutation_p(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> float:
    observed = spearmanr(x, y).statistic
    null = np.array([spearmanr(x, rng.permutation(y)).statistic for _ in range(N_PERMUTATIONS)])
    return float((np.sum(null >= observed) + 1) / (N_PERMUTATIONS + 1))


def main() -> int:
    rng = np.random.default_rng(0)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    graph_maps = rb.load_graph_and_maps(device)
    braf_idx = rb.protein_index_for_symbol(graph_maps, rb.BRAF_SYMBOL)
    pairs = rb.braf_case_study_pairs(graph_maps)

    split = pd.read_csv(SPLIT_CSV, usecols=["sanger_model_id", "split"], dtype=str).drop_duplicates()
    pairs = pairs.merge(split, on="sanger_model_id", how="left")
    print(f"[split] {pairs.drop_duplicates('sanger_model_id')['split'].value_counts().to_dict()} cell lines")

    model = rb.load_trained_model(graph_maps, device)
    cells, drugs = pairs["cell_idx"].tolist(), pairs["drug_idx"].tolist()
    scores = rb.score_pairs(model, graph_maps, cells, drugs).cpu().numpy()
    with torch.no_grad():
        predictions = model(
            graph_maps["x_dict"], graph_maps["edge_index_dict"],
            torch.tensor(cells, device=device), torch.tensor(drugs, device=device),
        ).cpu().numpy()

    ln_ic50 = pairs["ln_ic50"].to_numpy()
    braf = scores[:, braf_idx]

    varying = np.flatnonzero(scores.std(axis=0) > 0)
    rhos = np.array([spearmanr(scores[:, j], ln_ic50).statistic for j in varying])
    braf_rho = spearmanr(braf, ln_ic50).statistic
    braf_rank = int((rhos > braf_rho).sum()) + 1

    rows = []
    for label, mask in [
        ("all", np.ones(len(pairs), dtype=bool)),
        ("train cell lines", (pairs["split"] == "train").to_numpy()),
        ("held-out cell lines (val+test)", pairs["split"].isin(["val", "test"]).to_numpy()),
    ]:
        rows.append({
            "subset": label,
            "n_pairs": int(mask.sum()),
            "n_cell_lines": pairs.loc[mask, "sanger_model_id"].nunique(),
            "spearman_braf_score_vs_ln_ic50": spearmanr(braf[mask], ln_ic50[mask]).statistic,
            "permutation_p": permutation_p(braf[mask], ln_ic50[mask], rng),
            "spearman_prediction_vs_ln_ic50": spearmanr(predictions[mask], ln_ic50[mask]).statistic,
            "spearman_braf_score_vs_prediction": spearmanr(braf[mask], predictions[mask]).statistic,
            "braf_hit_rate": float((braf[mask] < 0).mean()),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print("BRAF specificity and out-of-sample check (trained model)")
    print("=" * 100)
    print(f"Proteins with a non-constant score across the {len(pairs)} pairs: {len(varying)}")
    print(f"BRAF Spearman vs ln(IC50) = {braf_rho:.4f} -> rank {braf_rank} of {len(varying)} "
          f"(top {100 * braf_rank / len(varying):.1f}%)")
    print(f"Other proteins: median rho {np.median(rhos):.4f}, 95th pct {np.percentile(rhos, 95):.4f}, "
          f"share with rho >= 0.5: {(rhos >= 0.5).mean():.1%}")
    with pd.option_context("display.width", 200, "display.max_columns", None, "display.float_format", "{:.4f}".format):
        print("\n" + summary.T.to_string(header=False))
    print(f"\nSaved {RESULTS_CSV.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
