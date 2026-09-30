# Target-Edge Ablation — Does the Model Need to Be Told a Drug's Target?

**Scripts:** [`build_ablated_graph.py`](../experiments/16_target_edge_ablation/build_ablated_graph.py), [`train_ablated.py`](../experiments/16_target_edge_ablation/train_ablated.py), [`evaluate_ablation.py`](../experiments/16_target_edge_ablation/evaluate_ablation.py)
**Raw results:** `experiments/16_target_edge_ablation/rmse_comparison.csv`, `attribution_comparison.csv`, `train_history_{original,ablated}.csv`
**Companion:** [`15_interpretability_validation_results.md`](./15_interpretability_validation_results.md)
**Date:** 2026-09-23 · **Runtime:** ~6 min training (2 runs, GTX 1650) + ~5 min evaluation

## What this covers

[Experiment 15](./15_interpretability_validation_results.md) found that the model's BRAF attribution tracks measured drug sensitivity (ρ = 0.65), which untrained models never do. But the graph wires `drug —targets→ BRAF` directly. The question here is the strongest version of the test: **if the model is never told the target, does it still work it out from the IC50 labels?**

We removed the `targets` edge from 7 probe drugs:

| Target | Drugs | Driver-mutant case-study pairs |
|---|---|---|
| BRAF | PLX-4720, SB590885, Dabrafenib | 157 (53 cell lines) |
| EGFR | Gefitinib, Erlotinib, AZD3759, Osimertinib | 48 (12 cell lines) |

That is 7 of 683 target edges. Nothing else changed. Drugs connect to proteins only through `targets`, so each probe drug is now **fully cut off from the protein graph**: its embedding comes from its Morgan fingerprint alone. The probe drugs' IC50 labels **stay in training**, so this is not leave-drugs-out. The model still sees how these drugs behave; it just isn't told what they bind.

Four models, all scored on the same test split (`data/raw/aligned_ic50_pairs.csv`, 16,840 test rows, cell-line-grouped):

| Model | Weights | Graph |
|---|---|---|
| `production` | `models/checkpoints/best_hetero_gnn.pt` | original |
| `zero_shot` | same production weights | **ablated** (edge vanishes, no chance to adapt) |
| `control_seed42` | retrained, seed 42 | original |
| `ablated_seed42` | retrained, seed 42 | **ablated** |

**The clean comparison is `control_seed42` vs `ablated_seed42`**: same seed, same data, same recipe as `05_train.py`, differing only in the 7 edges. The production checkpoint can't be the control, because `05_train.py` never sets a seed.

## Read this first: this cannot improve the model

Removing edges removes information. The best possible outcome is "no cost" (the edge was redundant). There is no mechanism by which deleting correct annotations makes a better model. This experiment measures **how much the model depends on target curation**, which matters for drugs that don't have one.

## Result 1 — RMSE

| Test subset | Rows | production | zero_shot | control_seed42 | ablated_seed42 |
|---|---|---|---|---|---|
| All test rows | 16,840 | 1.3035 | 1.3175 | **1.3002** | 1.3121 |
| Non-probe drugs | 16,284 | 1.3052 | 1.3052 | 1.3048 | 1.3159 |
| **Probe drugs** | 556 | 1.2523 | **1.6387** | **1.1589** | **1.1976** |
| BRAF probes | 236 | 1.1848 | 1.4703 | 1.0794 | 1.1362 |
| EGFR probes | 320 | 1.2998 | 1.7526 | 1.2142 | 1.2409 |

Training: control stopped at epoch 50 (best 35), ablated at epoch 40 (best 25). Best validation RMSE was 1.2744 vs 1.2754. Convergence was normal in both.

**Paired bootstrap over test cell lines** (2,000 resamples, ablated − control):

| Subset | Δ RMSE | 95% CI | P(Δ ≤ 0) |
|---|---|---|---|
| Probe drugs | **+0.039** | −0.009 to +0.092 | 0.064 |
| Non-probe drugs (edges untouched) | +0.011 | −0.009 to +0.031 | 0.145 |

What this shows:

1. **The trained model leans heavily on the edge.** Pull it out of the production model without retraining (`zero_shot`) and probe-drug RMSE jumps from 1.25 to **1.64**, with R² falling from 0.50 to 0.15. The non-probe rows don't move at all (1.3052 → 1.3052), which confirms the ablation touched only what it should.
2. **Retraining recovers almost all of it.** Trained without the edge from the start, probe RMSE is 1.198, only +0.039 worse than its same-seed control. That gap is suggestive but **not statistically established** (the CI crosses zero, p ≈ 0.06). For comparison, the untouched non-probe drugs moved by +0.011 between the same two runs. That is the noise you get just from a slightly different graph steering training differently. So the fingerprint plus the cell-line mutation path carries most of what the target edge provided.
3. **Overall RMSE is unaffected** (1.3002 vs 1.3121, inside the ±0.013–0.029 seed band in [13_seed_variance](./13_seed_variance_results.md)), as expected when only 3% of test rows involve a probe drug.

(Both seed-42 retrains happen to beat production on probe drugs: 1.16 and 1.20 vs 1.25. That is seed-to-seed variation in the production checkpoint, not a finding.)

## Result 2 — Does the attribution still track sensitivity without the edge?

Same statistic as Experiment 15 Check 1: Spearman correlation between the target protein's gradient × embedding score and the **measured** ln(IC50), across driver-mutant cell lines on that gene's inhibitors. The null is **50 untrained models on the ablated graph**, i.e. what wiring alone produces once the edge is gone.

**BRAF** (157 pairs; 15 of the 53 cell lines are held out from training):

| Model | ρ, all pairs | ρ, held-out lines only | Untrained models ≥ this ρ | Mean \|BRAF score\| |
|---|---|---|---|---|
| production | 0.647 | 0.543 | 0 / 50 | 0.510 |
| control_seed42 | 0.598 | 0.512 | 0 / 50 | 0.788 |
| zero_shot | 0.359 | 0.289 | 0 / 50 | 0.088 |
| **ablated_seed42** | **0.442** | **0.092** | 0 / 50 | 0.293 |
| untrained ×50 (ablated graph) | −0.011 ± 0.174 | −0.005 ± 0.236 | — | 0.004 |

**EGFR** (48 pairs; 6 of 12 cell lines held out): no model clearly beats the untrained floor. Production is at ρ = 0.29 (4 of 50 untrained models match it), and the ablated retrain is at 0.33 (3 of 50). With 12 cell lines this is underpowered, and **inconclusive in either direction**.

What this shows:

1. **The seed-42 control reproduces Experiment 15** (ρ 0.60 all / 0.51 held-out vs production's 0.65 / 0.54). The original finding is not a one-off from a lucky checkpoint.
2. **Without the edge, the model still learns *something* BRAF-specific on the lines it trained on.** ρ = 0.44, and none of 50 untrained models reach it. Because the drug is cut off from every protein, this can only come from the final layers learning to combine "this cell line has a BRAF mutation" (via the mutation edge) with "this drug's fingerprint looks like a RAF inhibitor". That combination was learned from IC50 labels, not wired in.
3. **But it does not generalise.** On held-out cell lines, the ablated model's ρ collapses to **0.09**, well within the untrained spread (±0.24). With the edge, held-out ρ is 0.51–0.54. So the out-of-sample BRAF finding from Experiment 15 **depends on the target edge**. Without it, the model fits the pattern on training lines but doesn't carry it to new ones.
4. **"BRAF is ranked #1" is uninformative again.** BRAF stays at median rank 1 in every model, and untrained models on the ablated graph still put it at median rank 3. BRAF is directly mutation-linked to every one of these cell lines, so it gets a large gradient under any weights.

## Verdict

| Question | Answer |
|---|---|
| Does removing target edges improve RMSE? | No, and it can't by design. |
| Does removing them hurt RMSE? | Barely, once retrained. +0.039 on probe drugs (not significant, p ≈ 0.06); overall unchanged. Pulled from an already-trained model without retraining, it hurts a lot (+0.39). |
| Can the model infer a drug's target from IC50 labels alone? | **Partly, on training cell lines only.** BRAF attribution still tracks sensitivity in-sample (ρ 0.44 vs untrained ~0), but not on held-out lines (ρ 0.09). |
| Is Experiment 15's held-out result learned independently of the wiring? | **No.** It needs the target edge. With the edge, the model learns *how much* BRAF matters from data, which untrained models can't do. Without the edge, that learning doesn't generalise. |

**Recommended wording for the report:** the model uses curated drug–target edges to learn which cell lines will respond, and this learned relationship generalises to unseen cell lines. It is not wired in, because untrained models on the same graph never show it. However, the model does not independently rediscover targets from response data: with the target edge removed it fits the pattern on training lines but loses it on held-out ones. Prediction accuracy is robust to a missing target annotation (probe-drug RMSE +0.04 after retraining), so the edge matters more for interpretability than for accuracy.

**Practical implication.** Missing target annotations cost little accuracy once the model is trained without them. So collecting more target annotations (e.g. ChEMBL, DrugBank) is unlikely to be a meaningful RMSE lever. It would mainly make gene-level explanations for those drugs trustworthy. The inverse is the operational risk: a model trained *with* an edge breaks badly if that edge later disappears (zero-shot +0.39). Graph edits should always be followed by a retrain.

**Limits.** One seed per arm. Two gene families, and EGFR is underpowered. Gradient × embedding is the only attribution method used. Both seed-42 checkpoints live under `experiments/16_target_edge_ablation/checkpoints/`; the production checkpoint and graph are untouched.
