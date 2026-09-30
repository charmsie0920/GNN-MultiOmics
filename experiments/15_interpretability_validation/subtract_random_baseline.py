"""Check 2 -- what is left of the gene ranking once untrained-model attribution is removed?

Today's gene panel ranks proteins by |gradient x embedding| from the trained
model. Part of that ranking is structural: proteins one hop from the drug or
cell line get large gradients under *any* weights. Subtracting what untrained
controls attribute isolates what training changed, and enrichment on that
residual shows whether the learned part still points at sensible biology.

Each model's score vector is first scaled to a share of its own total
|attribution|, because untrained and trained weights put gradients on
different scales and a raw subtraction would be dominated by whichever is
larger. Two residuals are reported:
- `residual`: trained share - mean control share (the plain subtraction).
- `zscore`:   (trained share - mean control share) / control std, i.e. how
  unusual each protein's trained score is relative to the controls.

Run from the repository root:
    python experiments/15_interpretability_validation/subtract_random_baseline.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import random_baseline as rb  # noqa: E402
from src.interpretation.enrichment import enrich  # noqa: E402
from src.interpretation.genes import annotate_genes, load_symbol_map, target_recovery  # noqa: E402
from src.interpretation.hetero_gnn_attribution import evidence_sets  # noqa: E402

RESULTS_CSV = rb.OUT_DIR / "subtract_random_baseline_results.csv"
TERMS_CSV = rb.OUT_DIR / "subtract_random_baseline_terms.csv"
TOP_K = 50
TERMS_SHOWN = 8
MAPK_KEYWORDS = ("MAPK", "RAF", "ERK", "MEK", "RAS")


def share(scores: np.ndarray) -> np.ndarray:
    return scores / np.abs(scores).sum(axis=-1, keepdims=True)


def representative_pairs(pairs: pd.DataFrame) -> pd.DataFrame:
    """The most sensitive BRAF-mutant line per BRAF inhibitor -- the clearest case of the biology."""
    return pairs.loc[pairs.groupby("drug_id")["ln_ic50"].idxmin()].reset_index(drop=True)


def main() -> int:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    graph_maps = rb.load_graph_and_maps(device)
    node_ids = graph_maps["protein_node_ids"]
    symbol_map = load_symbol_map(rb.REPO_ROOT / "data" / "processed" / "protein_symbol_map.csv")
    known_symbols = {s.upper() for s in symbol_map.values() if s}

    chosen = representative_pairs(rb.braf_case_study_pairs(graph_maps))
    cells, drugs = chosen["cell_idx"].tolist(), chosen["drug_idx"].tolist()
    print(f"[setup] pairs: {list(zip(chosen['cell_line_name'], chosen['drug_name'], chosen['ln_ic50'].round(2)))}")

    trained_raw = rb.score_pairs(rb.load_trained_model(graph_maps, device), graph_maps, cells, drugs).cpu().numpy()
    random_raw = np.stack(
        [rb.score_pairs(m, graph_maps, cells, drugs).cpu().numpy() for m in rb.build_random_models(graph_maps, device)]
    )  # [n_models, n_pairs, n_proteins]
    print(f"[scale] mean total |attribution| trained={np.abs(trained_raw).sum(-1).mean():.3g} "
          f"controls={np.abs(random_raw).sum(-1).mean():.3g}")

    trained = share(trained_raw)
    controls = share(random_raw)
    control_mean = controls.mean(axis=0)
    # Proteins outside every model's 2-hop radius have zero std; they carry no signal either way.
    control_std = np.where(controls.std(axis=0) > 0, controls.std(axis=0), np.inf)
    variants = {
        "raw": trained,
        "control_mean": control_mean,
        "residual": trained - control_mean,
        "zscore": (trained - control_mean) / control_std,
    }

    summary_rows, term_rows = [], []
    for p, pair in chosen.iterrows():
        drivers, targets = evidence_sets(graph_maps["edge_index_dict"], cells[p], drugs[p], node_ids)
        top_sets = {}
        print("\n" + "=" * 100)
        print(f"{pair['cell_line_name']} x {pair['drug_name']}  (ln IC50 {pair['ln_ic50']:.2f})")
        print("=" * 100)
        for name, matrix in variants.items():
            genes = annotate_genes(
                matrix[p].tolist(), node_ids, symbol_map, TOP_K, driver_proteins=drivers, target_proteins=targets
            )
            top_sets[name] = {g["gene_symbol"] for g in genes}
            recovery = target_recovery(genes, "BRAF", known_symbols)
            braf_row = next((g for g in genes if g["gene_symbol"] == "BRAF"), None)
            enrichment = enrich([g["gene_symbol"] for g in genes])
            terms = enrichment["terms"]
            significant = [t for t in terms if t["adjusted_p_value"] < 0.05]
            mapk_terms = [t for t in significant if any(k in t["term"].upper() for k in MAPK_KEYWORDS)]

            summary_rows.append({
                "cell_line": pair["cell_line_name"], "drug": pair["drug_name"], "ln_ic50": pair["ln_ic50"],
                "variant": name,
                "braf_rank": recovery["recovered"][0]["rank"] if recovery["recovered"] else None,
                "braf_direction": braf_row["direction"] if braf_row else None,
                "drug_targets_in_top": sum(g["is_drug_target"] for g in genes),
                "driver_mutations_in_top": sum(g["is_driver_mutation"] for g in genes),
                "overlap_with_raw_top": len(top_sets[name] & top_sets["raw"]),
                "enrichment_status": enrichment["status"],
                "significant_terms": len(significant),
                "significant_mapk_terms": len(mapk_terms),
                "top_term": terms[0]["term"] if terms else None,
                "top_term_adj_p": terms[0]["adjusted_p_value"] if terms else None,
                "top_genes": ";".join(g["gene_symbol"] for g in genes[:15]),
            })
            for t in terms[:TERMS_SHOWN]:
                term_rows.append({"cell_line": pair["cell_line_name"], "drug": pair["drug_name"], "variant": name,
                                  "term": t["term"], "library": t["library"],
                                  "adjusted_p_value": t["adjusted_p_value"], "overlap": t["overlap"]})

            r = summary_rows[-1]
            print(f"\n[{name}] BRAF rank={r['braf_rank']} ({r['braf_direction']})  "
                  f"targets/drivers in top{TOP_K}={r['drug_targets_in_top']}/{r['driver_mutations_in_top']}  "
                  f"overlap w/ raw={r['overlap_with_raw_top']}  sig terms={r['significant_terms']} "
                  f"(MAPK-ish {r['significant_mapk_terms']})  [{enrichment['status']}]")
            print(f"  top genes: {r['top_genes']}")
            for t in terms[:5]:
                print(f"  {t['adjusted_p_value']:.2e}  {t['term']}  ({t['library']})")

    pd.DataFrame(summary_rows).to_csv(RESULTS_CSV, index=False)
    pd.DataFrame(term_rows).to_csv(TERMS_CSV, index=False)
    print(f"\nSaved {RESULTS_CSV.name}, {TERMS_CSV.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
