# Ablation Study — Base + One Component at a Time

**Script:** [`src/final_model/run_ablation.py`](./final_model/run_ablation.py)
**Model:** [`src/models/mographdrp_aligned.py`](./models/mographdrp_aligned.py) (frozen Phase 0 baseline, imported by the script)
**Raw results:** `src/final_model/results/ablation_results.csv`
**Companions:** [`results.md`](../docs/results.md), [`13_seed_variance_results.md`](../docs/13_seed_variance_results.md), [`09_split_protocol_comparison.md`](../docs/09_split_protocol_comparison.md)
**Date:** 2026-10-03 · **Runtime:** 6.5 h (24 runs, 23,523 s) · **Hardware:** one laptop GTX 1650 (4 GB), torch 2.13.0+cu130
**Status:** Phase 2 filled (components a-e, 3 seeds, grouped split). Component f is not built; sections that depend on it, on the random split or on the final-model choice are still `TBD`.

> Anything marked *(report)* is meant to be lifted into the final report.

---

## 1. What this study answers

*(report)* One or two sentences. Suggested:

> Starting from the simplest pipeline that functions, which individual
> components improve prediction of ln(IC50) on unseen cell lines, by how much,
> and does combining them give the best model?

Requirement from the lecturer: **Base + a**, **Base + b**, **Base + c**, ...,
where each letter is one component, and Base + everything is the final model.

## 2. What is held fixed across every row

These must be identical for every rung, otherwise a delta is not attributable.

| Item | Setting |
|---|---|
| Pairs | 111,799 (cell line, drug) pairs; 531 cell lines (371 / 80 / 80 train / val / test; 532 have omics, one has no pair after the fingerprint filter); 240 drug IDs |
| Split | grouped by cell line, 70/15/15, `random_state=42` (77,668 / 17,291 / 16,840) |
| Omics preprocessing | per-modality scaling + PCA to 128 |
| Optimiser | Adam, lr 1e-4, weight decay **0** (Phase 0 frozen base) |
| Batch size / dropout | 128 / 0.4 |
| Epochs / early stopping | 200 epochs, **no early stopping** (Phase 0 frozen base) |
| LR schedule | **constant** (Phase 0 frozen base) |
| Checkpoint | best validation RMSE, in ln(IC50) units |
| Seeds | 42, 43, 44 (3 per row) |
| Hardware | all 24 runs on one GTX 1650, one process, commit `66205f8` |
| Metrics | RMSE, MAE, R², PCC, SCC; AUC and F1 at the shared median threshold |

Early-stopping and LR-scheduler patience are both 200, so neither fires within
the 200-epoch budget. The model is the frozen `MoGraphDRPAligned` (BAN
bilinear head with 3 heads, (512, 128) predictor, gated-sum drug fusion).
Aligned to these settings on 2026-10-03; no ladder run predates that.

**What the seeds vary.** A seed changes weight initialisation and batch order
only. The grouped split is fixed at `random_state=42`, so every run uses the
same 80 test cell lines. The std below is training noise; variance from the
choice of test cell lines is not measured.

Reference floors on the same test fold:

| Floor | Test RMSE |
|---|---|
| Global training mean | 2.7690 |
| **Per-drug training mean** (no omics, no model) | **1.4889** |

**Noise band:** ~0.03 RMSE, from [13_seed_variance](../docs/13_seed_variance_results.md).
This study's own seeds give a std of 0.0105 for `base` and 0.004-0.022 for the
other rungs (`base+bilinear` excepted, 0.49). Three seeds are too few to
replace the 0.03 figure, so it is kept. A delta smaller than the band is
reported as "no measurable effect".

## 3. The base

*(report)* Describe it as a pipeline, then justify each choice as "simplest".

| Stage | Base choice | Why this is the simplest functioning option |
|---|---|---|
| Omics inputs | GE + Mut_CNV | TBD (two inputs are the minimum for fusion to apply; closest to the benchmark's inputs) |
| Cell encoder | one 2-layer branch per omics, concatenated | TBD |
| Drug encoder | Morgan fingerprint (radius 2, 2048 bits) → MLP | TBD |
| Cell–drug interaction | concatenation | TBD |
| Predictor | MLP 512 → 128 → 1 (Linear → ReLU → BatchNorm per layer) | TBD |
| Targets | raw ln(IC50) | TBD |

Parameters: 3,184,257 (from the results CSV). Diagram: TBD (`docs/figures/`).

## 4. The components

One row per component. "Switch" is the single field it changes in `BASE`.

| ID | Component | Switch | Replaces / adds | Origin | Hypothesis (why it should help) |
|---|---|---|---|---|---|
| a | Cross-attention fusion | `fusion=cross_attention` | replaces concat of omics branches | this project | TBD |
| b | Molecular-graph GCN | `drug_mode=both` | adds a GCN over the drug's atoms beside the fingerprint | MoGraphDRP §2.2.3 | TBD |
| c | Bilinear head | `head=bilinear` | replaces concat of cell and drug vectors | MoGraphDRP §2.3 | TBD |
| d | Proteomics | `modalities=+Proteomics` | adds a third omics branch | this project | TBD |
| e | Per-drug target standardisation | `standardize_targets=True` | trains on per-drug z-scores, inverted before scoring | this project | TBD |
| f | PPI GNN + pair-specific attention | TBD | TBD | this project (novelty component) | TBD |

Known caveats to state honestly in the report:

- **a:** each modality is one vector, so attention weights are always 1.0. As
  built, this rung measures the extra projection layers, not attention.
  Decided 2026-10-03: reported with this caveat, not fixed.
- **e:** not usable under leave-drugs-out (unseen drugs have no statistics).
- **f:** covers only pairs whose drug has a known target: 71.75% of pairs
  (80,213 of 111,799; test 12,119 of 16,840). 169 of the 240 drugs have a
  target and 71 have none. (The "130 drugs with no target" figure counts all
  498 drugs in the graph, most of which are not in these pairs.)

## 5. Main result — grouped split

*(report)* This is the central table. Mean ± std over 3 seeds (42, 43, 44),
test fold, grouped split.

| Config | Components | RMSE | Δ vs base | Beyond noise? | PCC | R² | Gain over per-drug mean | Params | Fit (s) |
|---|---|---|---|---|---|---|---|---|---|
| per-drug mean | — | 1.4889 | — | — | — | — | 0 | 0 | — |
| `base` | none | 1.3232 ± 0.0105 | 0 | — | 0.8797 | 0.7706 | +0.1657 | 3,184,257 | 542 |
| `base+cross_attention` | a | 1.4144 ± 0.0040 | +0.0912 | yes, worse | 0.8608 | 0.7379 | +0.0744 | 3,382,273 | 814 |
| `base+mol_graph` | b | 1.3319 ± 0.0223 | +0.0087 | no | 0.8806 | 0.7676 | +0.1570 | 3,509,997 | 1,159 |
| `base+bilinear` | c | 1.8343 ± 0.4946 | +0.5111 | yes, worse | 0.8545 | 0.5378 | −0.3455 | 3,202,822 | 648 |
| `base+proteomics` | d | 1.2930 ± 0.0169 | −0.0302 | borderline | 0.8871 | 0.7810 | +0.1959 | 3,316,225 | 582 |
| `base+std_targets` | e | 1.3108 ± 0.0195 | −0.0124 | no | 0.8866 | 0.7749 | +0.1781 | 3,184,257 | 525 |
| `base+pair` | f | not built | | | | | | | |
| `aligned` | b + c | 1.3416 ± 0.0088 | +0.0184 | no | 0.8751 | 0.7642 | +0.1473 | 3,528,562 | 1,254 |
| `full` | a-e | 1.3394 ± 0.0081 | +0.0162 | no | 0.8789 | 0.7650 | +0.1495 | 4,270,706 | 2,309 |

Negative Δ means better than base. Fit is the mean training time per run.

**Headline:** no component clearly beats `base`. Proteomics is the only one
that is better on every seed, and its mean gain (0.0302) sits on the edge of
the ~0.03 band. `full` is not better than `base` (+0.0162, worse on all three
seeds) and is 0.046 worse than `base+proteomics`.

Per-seed test RMSE:

| Config | Seed 42 | Seed 43 | Seed 44 |
|---|---|---|---|
| `base` | 1.3141 | 1.3347 | 1.3208 |
| `base+cross_attention` | 1.4114 | 1.4130 | 1.4189 |
| `base+mol_graph` | 1.3185 | 1.3197 | 1.3577 |
| `base+bilinear` | 2.4003 | 1.6182 | 1.4846 |
| `base+proteomics` | 1.2791 | 1.3118 | 1.2881 |
| `base+std_targets` | 1.2937 | 1.3066 | 1.3320 |
| `aligned` | 1.3378 | 1.3353 | 1.3517 |
| `full` | 1.3329 | 1.3485 | 1.3369 |

**Is each delta real?** No significance test was run. The verdicts rest on the
mean delta against the ~0.03 band and on how many seeds agree in sign. A
paired bootstrap over test cell lines, as in
[16](../docs/16_target_edge_ablation_results.md), is still to do.

| Component | Δ RMSE | Seeds where it beat base | Per-seed Δ (42 / 43 / 44) | Verdict |
|---|---|---|---|---|
| a | +0.0912 | 0 / 3 | +0.097 / +0.078 / +0.098 | hurts — **tested-and-rejected** |
| b | +0.0087 | 1 / 3 | +0.004 / −0.015 / +0.037 | no effect — **tested-and-rejected** |
| c | +0.5111 | 0 / 3 | +1.086 / +0.284 / +0.164 | hurts, unstable — **tested-and-rejected** as configured |
| d | −0.0302 | 3 / 3 | −0.035 / −0.023 / −0.033 | borderline: consistent in sign, at the band edge |
| e | −0.0124 | 2 / 3 | −0.020 / −0.028 / +0.011 | no effect — **tested-and-rejected** |
| f | not built | | | |

**Validation and test disagree on the ranking.** Mean validation RMSE:

| Config | Val RMSE | Δ val vs base | Δ test vs base |
|---|---|---|---|
| `base+std_targets` | 1.2580 ± 0.0097 | −0.0368 | −0.0124 |
| `aligned` | 1.2703 ± 0.0037 | −0.0245 | +0.0184 |
| `base+proteomics` | 1.2808 ± 0.0066 | −0.0140 | −0.0302 |
| `base+mol_graph` | 1.2883 ± 0.0123 | −0.0065 | +0.0087 |
| `full` | 1.2904 ± 0.0072 | −0.0044 | +0.0162 |
| `base` | 1.2948 ± 0.0064 | 0 | 0 |
| `base+cross_attention` | 1.3240 ± 0.0101 | +0.0292 | +0.0912 |
| `base+bilinear` | 1.7901 ± 0.4697 | +0.4953 | +0.5111 |

`aligned` is second-best on validation and worse than base on test;
`std_targets` gains 0.037 on validation and 0.012 on test. Validation and test
are two different fixed sets of 80 cell lines, so differences of this size
between rungs are within what the choice of cell lines can move. This limits
how firmly any rung can be ranked from one split.

## 6. Per-component findings

*(report)* One short block per component.

### a. Cross-attention fusion — tested-and-rejected
- **Result:** 1.4144 ± 0.0040, Δ +0.0912. Worse than base on all 3 seeds, three
  times the noise band.
- **Interpretation:** this does not test attention. Each modality is one
  vector, so the attention weights are always 1.0, and the rung *replaces* the
  per-omics branches with projection layers. The result says that replacement
  is worse than the branches; it does not show that attention fails.
- **Why (evidence):** best epoch was 20-42 against 72-122 for base, and
  validation RMSE (1.3240) is 0.09 better than test (1.4144), the largest
  val-test gap of any stable rung. Both point to earlier overfitting. No
  further diagnosis was run.
- **Cost:** +198,016 parameters, 1.5× training time.

### b. Molecular-graph GCN — tested-and-rejected
- **Result:** 1.3319 ± 0.0223, Δ +0.0087. Inside the noise band; better than
  base on 1 seed of 3.
- **Interpretation:** adding the atom-level GCN beside the Morgan fingerprint
  gives no measurable gain on unseen cell lines.
- **Why (evidence):** not diagnosed. Consistent with the test being
  cell-line-grouped: every drug is seen in training, so a richer drug encoder
  has little to add. That reading is untested here.
- **Cost:** +325,740 parameters, 2.1× training time.
- Compare with [11_molecular_graph_results](../docs/11_molecular_graph_results.md), where the single-run ranking reversed across seeds. The same sign flip appears here (−0.015 on seed 43, +0.037 on seed 44).

### c. Bilinear head — tested-and-rejected as configured
- **Result:** 1.8343 ± 0.4946, Δ +0.5111. Per seed 2.4003 / 1.6182 / 1.4846:
  worse than base on every seed, worse than the per-drug-mean floor on two,
  and level with it on the third.
- **Interpretation:** on a fingerprint-only drug encoder the BAN head does not
  train to a usable model within 200 epochs. This is a statement about this
  configuration and budget, not about bilinear interaction in general:
  `aligned` (the same head plus the molecular graph) reaches 1.3416 ± 0.0088.
- **Why (evidence):**
  - Best epoch was 187 / 199 / 200, so validation RMSE was still falling when
    training stopped. The rung is not converged.
  - Validation RMSE was erratic throughout (mostly 3.0-4.1 on seed 42).
  - Ranking is intact while absolute error is not: PCC 0.85 on every seed
    against MAE 2.08 / 1.27 / 1.14. That pattern is an offset or scale error
    in the predictions, not missing signal.
  - Cause not diagnosed. The head and the predictor both contain BatchNorm,
    which behaves differently in training and evaluation; that is a
    hypothesis, not a finding.
- **Cost:** +18,565 parameters, 1.2× training time.
- The head was ported from the MoGraphDRP repo in Phase 0 (2026-10-01) and
  validated only under the random split (0.978), paired with the molecular
  graph. It was not re-checked against their code for this study.

### d. Proteomics — borderline
- **Result:** 1.2930 ± 0.0169, Δ −0.0302. Better than base on all 3 seeds
  (−0.035 / −0.023 / −0.033). The mean gain equals the ~0.03 band.
- **Interpretation:** the only component with a consistent gain. It is not a
  clear pass: one seed (43) is inside the band, and on validation the gain is
  smaller (−0.014).
- **Why (evidence):** not diagnosed beyond the numbers above.
- **Cost:** +131,968 parameters, 1.07× training time.
- Earlier flat-model result to reconcile: Mut_CNV hurt (E32 vs E04 in [results.md](../docs/results.md)). Not tested here: there is no rung without Mut_CNV.

### e. Per-drug target standardisation — tested-and-rejected
- **Result:** 1.3108 ± 0.0195, Δ −0.0124. Inside the noise band; better than
  base on 2 seeds of 3, worse on seed 44 (+0.011).
- **Interpretation:** no measurable effect on test RMSE for this model. It has
  the best validation RMSE of any rung (1.2580, −0.037 vs base), but that gain
  does not carry to the test cell lines.
- **Why (evidence):** best epoch was 14-36, among the earliest of the single
  components, so it fits fast and then overfits. Not diagnosed further.
- **Cost:** no extra parameters, same training time (0.97×).
- Earlier single-seed result on HeteroIC50GNN: 1.3449 → 1.2917 (−0.053). That gain does not reproduce on this model over 3 seeds.

### f. PPI GNN + pair-specific attention
- **Result / Interpretation / Why / Cost:** TBD
- **GNN depth ablation** (does message passing matter?):

| Message-passing layers | RMSE (all pairs) | RMSE (pairs with a known target) |
|---|---|---|
| 0 (attention on raw embeddings, no GNN) | TBD | TBD |
| 1 | TBD | TBD |
| 2 | TBD | TBD |

- **Coverage:** TBD % of test pairs have a target; report both subsets.

## 7. Do the components combine?

*(report)* Components can overlap or interfere, so the sum of single gains
rarely equals the gain of the full model.

| Quantity | RMSE change vs base |
|---|---|
| Sum of the individual deltas (a + b + c + d + e) | +0.5685 |
| Same sum without c (bilinear), whose single-rung result is unstable | +0.0574 |
| `full` measured | +0.0162 |
| Difference, `full` − sum of all five | −0.5523 |

- The components are **not additive**. The sum is dominated by bilinear's
  +0.51, which does not appear inside `full`: with the molecular graph
  present, the head trains normally (`aligned` +0.0184, `full` +0.0162).
- Does `full` beat the best single component beyond noise? **No.** `full`
  (1.3394) is 0.046 worse than `base+proteomics` (1.2930) and is not better
  than `base` on any seed.
- `full` overfits early: best epoch 47 / 21 / 10.
- Is any component harmful inside `full`? Not measured. The leave-one-out
  check (`full` minus one) was not run:

| Config | RMSE | Δ vs full |
|---|---|---|
| `full` | 1.3394 ± 0.0081 | 0 |
| `full − a` | not run | |
| `full − b` | not run | |
| ... | | |

## 8. Final model

*(report)*

- **Chosen configuration:** not chosen yet. This is a team decision; the
  evidence from Phase 2 is:
  - best on test: `base+proteomics`, 1.2930 ± 0.0169 (PCC 0.8871, R² 0.7810);
  - best on validation: `base+std_targets`, 1.2580 ± 0.0097 (test 1.3108);
  - `full` (base + a-e) is not better than `base` and should not be presented
    as the final model on these numbers.
- **Components dropped and why:** a (+0.0912, hurts), b (+0.0087, no effect),
  c (+0.5111, hurts and unstable as configured), e (−0.0124, no effect). d
  (−0.0302) is borderline.
- **Final test RMSE / PCC / R²:** TBD once the configuration is chosen.
- **Gain over the per-drug mean floor:** every stable rung beats 1.4889; the
  largest gain is `base+proteomics` at +0.1959. `base` itself is +0.1657.

Selection must be made on **validation** RMSE; the test numbers above are
reported after the choice, not used to make it. Not done yet. Note that
validation and test rank the rungs differently (section 5), so a choice made
on validation would pick `std_targets`, whose test gain is inside the noise
band.

## 9. Sanity check against the benchmark — random split

Not a ranking. It checks that the reproduction (`aligned`) is faithful.

| Config | Random-split RMSE | Reference |
|---|---|---|
| `aligned` (Phase 0, GE+Proteomics) | **0.978 ± 0.002** (3 seeds) | MoGraphDRP without XGBoost: 0.9497 |
| `base` | TBD | — |
| `full` | TBD | MoGraphDRP with XGBoost: 0.6622 |

**Phase 0 result (2026-10-01).** Source:
`experiments/14_benchmark_alignment/benchmark_alignment_results.csv`.

| Seed | Test RMSE | PCC | R² | Best epoch |
|---|---|---|---|---|
| 42 | 0.9756 | 0.9395 | 0.8819 | 195 |
| 43 | 0.9799 | 0.9386 | 0.8808 | 196 |
| 44 | 0.9798 | 0.9389 | 0.8808 | 200 |
| **mean ± std** | **0.9784 ± 0.0025** | 0.9390 | 0.8812 | |

- Settings: BAN head (3 heads), (512, 128) predictor, gated-sum drug fusion,
  weight decay 0, no early stopping, 200 epochs at a constant LR, 3,528,562
  parameters. Mean-only floor on this split: 2.838. The per-drug-mean floor
  (1.4889) is for the grouped split and does not apply here.
- The gap to 0.9497 is 0.029, at the edge of the ~0.03 noise band. Remaining
  known differences are all data: 2 PCA omics vs their 4 gene-filtered ones,
  Morgan only vs 3 fingerprints, our atom features, different pairs. Whether
  they explain the gap has not been tested.
- Seeds vary initialisation only; the random split is fixed.
- Best epoch 195-200: the model is still improving when the 200-epoch budget ends.
- **Not the same config as the ladder's `aligned`.** This run used
  GE+Proteomics; the ladder's `aligned` (= `base+mol_graph+bilinear`) uses
  GE+Mut_CNV. The two are not directly comparable.

Protocol gap for the final model (random − grouped): TBD. Compare with the
0.347 measured in [09](../docs/09_split_protocol_comparison.md).

## 10. Comparison with the other models on the same pairs

Same 111,799 pairs, same grouped split. TBD seeds each.

| Model | RMSE | PCC | Notes |
|---|---|---|---|
| Per-drug mean | 1.4889 | — | floor |
| RF | TBD | TBD | |
| MLP | TBD | TBD | |
| CrossAttention (E10) | 1.3029 ± 0.0126 | TBD | 5 seeds, [13](../docs/13_seed_variance_results.md) |
| HeteroIC50GNN (E11) | TBD | TBD | single run 1.3301 |
| `aligned` (MoGraphDRP-style) | 1.3416 ± 0.0088 | 0.8751 | 3 seeds, GE+Mut_CNV, this study |
| `base` (this ladder) | 1.3232 ± 0.0105 | 0.8797 | 3 seeds, this study |
| **Final model** | TBD | TBD | |

## 11. Further checks (fill if run)

- **Leave-drugs-out** (unseen compounds, no standardisation): TBD vs the MLP's 1.8516 in [10](../docs/10_leave_drugs_out_results.md).
- **Interpretability:** does attention on BRAF track measured sensitivity against untrained models? TBD, see [15](../docs/15_interpretability_validation_results.md).

## 12. Limitations to state

Tick the ones that apply and add numbers.

- [x] Seeds: 3 per row. The split itself is fixed at seed 42, so variance from the choice of test cell lines is not measured.
- [x] No hyperparameter tuning per rung; all use the benchmark's Table 1 settings.
- [x] Fixed 200-epoch budget with no early stopping. `base+bilinear` was still improving at epoch 200, so its result is a lower bound on what the head can do, not its converged value. Every other rung peaked before the end (best epoch 10-180).
- [x] Validation and test rank the rungs differently (section 5); both are single fixed sets of 80 cell lines.
- [x] No significance test; verdicts use the ~0.03 band and per-seed sign.
- [x] Base uses PCA features, not gene-level features; differs from the benchmark's COSMIC filtering.
- [ ] Only 532 of GDSC2's 969 cell lines (tri-omics-complete subset).
- [x] Component a is degenerate as built (see §4).
- [x] Component f covers 71.75% of pairs (169 of 240 drugs have a known target).
- [x] The benchmark comparison is against a reimplementation, not the authors' code on their data.

## 13. Figures to produce

| Figure | Content | Status |
|---|---|---|
| Ladder bar chart | Δ RMSE vs base per component, error bars = std over seeds, noise band shaded | TBD |
| Architecture diagram | base in grey, each component highlighted | TBD |
| Additivity plot | sum of single gains vs `full` | TBD |
| Protocol comparison | grouped vs random for base / aligned / final | TBD |

## 14. Report wording (draft once results are in)

*(report)* Draft from the Phase 2 numbers. Revise once component f and the
final configuration exist.

> The base model reaches RMSE 1.3232 ± 0.0105 on unseen cell lines (3 seeds),
> against 1.4889 for a per-drug mean. Of the five components tested
> individually, none improved on the base clearly beyond seed-to-seed
> variation (~0.03): adding proteomics gave the only consistent gain
> (−0.030, better on all three seeds), per-drug target standardisation and the
> molecular-graph encoder had no measurable effect, and cross-attention fusion
> (+0.091) and the bilinear head on its own (+0.511, unstable) made the model
> worse. Combining all five gives 1.3394 ± 0.0081, which is not better than
> the base.

## 15. Reproduction

```
python "src/final_model/run_ablation.py" --list
python "src/final_model/run_ablation.py" --protocols grouped --seeds 42 43 44
python "src/final_model/run_ablation.py" --configs aligned base full --protocols random
```

The second command produced the 24 rows in this file (run 2026-10-03, 15:33 to
22:05, as one process on a GTX 1650; launch details in `CLAUDE.md`). The third
(random split for `base` and `full`) has not been run.

Commit hash of the code that produced the CSV: `66205f8`.
