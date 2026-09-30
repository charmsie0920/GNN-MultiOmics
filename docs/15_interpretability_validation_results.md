# Interpretability Validation — Learned Biology vs. Graph Wiring

**Scripts:** [`direction_sign_test.py`](../experiments/15_interpretability_validation/direction_sign_test.py), [`braf_specificity.py`](../experiments/15_interpretability_validation/braf_specificity.py), [`subtract_random_baseline.py`](../experiments/15_interpretability_validation/subtract_random_baseline.py) (shared setup in [`random_baseline.py`](../experiments/15_interpretability_validation/random_baseline.py))
**Raw results:** `experiments/15_interpretability_validation/*.csv`
**Model:** `models/checkpoints/best_hetero_gnn.pt` (`HeteroIC50GNN`), attribution = gradient × embedding (`src/interpretation/hetero_gnn_attribution.py`, unmodified)
**Date:** 2026-09-23 · **Runtime:** ~15 min on CPU

## What this covers

The gene panel's main validation so far has been *target recovery*: does a drug's GDSC target appear among its top-attributed genes? That test is weak, because the graph contains a `drug —targets→ protein` edge built from the same annotation. A drug's target sits one hop from the drug node, so it gets a large gradient under almost any weights. Recovering it may only mean the edge is there.

What the graph does **not** tell the model is *direction*. It knows a cell line carries a BRAF driver mutation and that PLX-4720 targets BRAF. It does not know whether that mutation makes the drug work better or worse, and only the IC50 labels carry that. So we test the model on direction and use **untrained controls** as the null. The controls are 50 `HeteroIC50GNN`s with the checkpoint's exact architecture and graph, initialised at seeds 1000–1049 and never trained. Whatever they show comes from the wiring and architecture alone.

**Case study:** every real GDSC2 measurement of a BRAF-driver-mutant cell line on a BRAF inhibitor. That is 53 cell lines × {PLX-4720, SB590885, Dabrafenib}, giving **157 pairs**. Known biology: these lines are *more* sensitive, so a model that learned the biology should give BRAF a negative (sensitising) attribution, most strongly where the measured sensitivity is highest.

## Check 1 — Direction of the BRAF attribution

| | All | PLX-4720 | SB590885 | Dabrafenib |
|---|---|---|---|---|
| Pairs | 157 | 53 | 53 | 51 |
| Trained: share of pairs where BRAF is sensitising | 0.707 | 0.604 | 0.774 | 0.745 |
| Untrained controls: mean ± SD of that share | 0.540 ± 0.361 | 0.539 ± 0.362 | 0.546 ± 0.367 | 0.535 ± 0.359 |
| Controls at least as "right" as trained (empirical p) | 0.44 | 0.46 | 0.42 | 0.42 |
| Trained: Spearman(BRAF score, measured ln IC50) | **0.647** | **0.648** | **0.714** | **0.664** |
| Controls: Spearman mean ± SD | −0.017 ± 0.190 | −0.017 ± 0.219 | −0.015 ± 0.213 | −0.014 ± 0.193 |
| Controls reaching the trained Spearman | **0 / 50** | 0 / 50 | 0 / 50 | 0 / 50 |

**The sign on its own does not separate trained from untrained.** The trained model calls BRAF sensitising in 71% of pairs. That would pass a naive binomial test (p < 1e-6), but that test is invalid here. Within one model the signs are strongly correlated across pairs, because the weights and the BRAF embedding are shared. The controls therefore land almost anywhere from 0% to 100% (8 of 50 are all-or-nothing), and 44% of them match or beat the trained model.

**Whether the BRAF score tracks measured sensitivity does separate them.** Across the 157 pairs, the more sensitive the lab measured a cell line, the more negative the trained model's BRAF attribution (ρ = 0.65). The controls average ρ ≈ 0, and none of the 50 reach it. The wiring is identical for all of these pairs (every one has the same mutation edge and the same target edge), so the wiring cannot produce this variation. It was learned.

### Follow-up: specific to BRAF, and out of sample?

(`braf_specificity.py`)

- **Specific.** Of the 5,672 proteins whose score varies across these pairs, BRAF's correlation with ln(IC50) ranks **4th (top 0.1%)**. The median protein is at ρ = 0.007 and only 0.4% of proteins reach ρ ≥ 0.5. So this is not a general effect where every attribution scales with the prediction.
- **Out of sample.** Using the checkpoint's own cell-line split (`data/raw/aligned_ic50_pairs.csv`):

| Subset | Cell lines | Pairs | ρ(BRAF score, ln IC50) | Permutation p | Sensitising share |
|---|---|---|---|---|---|
| Train | 38 | 113 | 0.682 | 1e-4 | 0.637 |
| **Held-out (val+test)** | **15** | **44** | **0.543** | **4e-4** | 0.886 |

The relationship holds on cell lines the model never trained on. It weakens somewhat (0.68 → 0.54), as expected when moving to held-out data.

One caveat: the BRAF score correlates ρ = 0.72 with the model's own prediction. That is consistent with BRAF being a main driver of these predictions. It also means the BRAF attribution is partly a readout of *how* sensitive the model thinks the line is, not an independent second signal.

## Check 2 — Subtracting the untrained baseline before enrichment

For the most sensitive BRAF-mutant line per drug, we compared four top-50 gene lists, each run through the existing `annotate_genes` → `target_recovery` → `enrich` pipeline. Before comparing, each model's scores were scaled to a share of its own total |attribution|, because untrained gradients are ~24× smaller.

- **raw**: today's panel (trained scores)
- **control_mean**: what the untrained controls alone rank highest
- **residual**: trained − control mean
- **zscore**: (trained − control mean) / control SD, i.e. how unusual each protein's trained score is

| Pair | Variant | BRAF rank | Overlap with raw top-50 | Top enriched term (adj. p) |
|---|---|---|---|---|
| SK-MEL-28 × PLX-4720 | raw | 1 | 50 | Disease / Oncogenic MAPK Signaling (7.8e-25) |
| | control_mean | 1 | 30 | Bladder cancer (2.0e-20); Melanoma, BRAF/RAF1 fusions in top 5 |
| | residual | 1 | 40 | Oncogenic MAPK Signaling (1.6e-24) |
| | zscore | 1 | 9 | Hepatitis B (2.7e-5); no MAPK term in top 5 |
| UACC-62 × SB590885 | raw | 1 | 50 | Oncogenic MAPK Signaling (1.5e-24) |
| | control_mean | 1 | 22 | Pathways in cancer (1.0e-22) |
| | residual | 1 | 31 | Oncogenic MAPK Signaling (2.1e-22) |
| | zscore | 1 | 12 | Signal Transduction (1.0e-11); PIP3/AKT |
| DU-4475 × Dabrafenib | raw | 1 | 50 | Oncogenic MAPK Signaling (9.0e-34) |
| | control_mean | 1 | 23 | Proteoglycans in cancer (1.0e-19); MAPK/BRAF in top 5 |
| | residual | 1 | 35 | BRAF/RAF1 Fusions (4.8e-31) |
| | zscore | 1 | 30 | BRAF/RAF1 Fusions (4.4e-24) |

Full term lists: `subtract_random_baseline_terms.csv`. (Every variant hits the 75-term cap, 25 per library × 3, so the count of significant terms is not informative.)

**Most of today's pathway panel is reproduced by untrained models.** The controls alone put BRAF at rank 1 and already produce strongly significant cancer/MAPK enrichment (top terms at adj. p ~1e-19 to 1e-22). They share 22–30 of the trained model's top 50 genes. So "top genes are enriched for MAPK signalling" and "BRAF is recovered" are **not** evidence of learning on their own. They mostly reflect the drug's 2-hop neighbourhood in STRING.

**After subtraction, the MAPK signal remains.** The plain residual keeps a MAPK / BRAF-fusion term as the #1 term for all three pairs (adj. p 1e-22 to 1e-31). The controls' #1 term is a generic cancer term in every case. The stricter z-score variant keeps it for Dabrafenib only. For PLX-4720 and SB590885 the z-score top genes shift to less canonical biology (AKT/PIP3, immune/metabolic terms at ~1e-5 to 1e-11). BRAF stays rank 1 in every variant.

## Interpretation

| Claim | Supported? |
|---|---|
| "The model recovers the drug's target" | Not on its own. Untrained models also rank BRAF first. |
| "Top genes are enriched for the right pathway" | Not on its own. Untrained models already get MAPK enrichment. |
| "The model learned the *direction* of the BRAF effect" | **Yes, in the graded sense.** BRAF attribution tracks measured sensitivity (ρ 0.65; 0/50 controls), ranks BRAF 4th of 5,672 proteins, and holds on held-out cell lines (ρ 0.54, p 4e-4). |
| "The model gets the sign right more often than chance" | Not demonstrable. Untrained sign rates vary too widely for 71% to stand out. |
| "Pathways the model *added* beyond wiring are biologically sensible" | Partly. MAPK survives plain subtraction on all three pairs, but survives the stricter z-score on only one of three. |

**Recommended wording for the report:** the model's BRAF attribution rises and falls with measured drug sensitivity across BRAF-mutant cell lines, including held-out ones. Untrained models on the same graph never reproduce this. Target recovery and pathway enrichment on their own should be presented as reflecting graph structure, not as evidence of learning.

**Limits.** One gene–drug-class case study, 3 pairs for Check 2, one attribution method. The strongest remaining test is option 3: remove the target edges for a set of drugs, retrain, and check whether the target is still recovered.
