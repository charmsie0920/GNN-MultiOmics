# Ablation Study — Base + One Component at a Time

**Script:** [`src/final_model/run_ablation.py`](./final_model/run_ablation.py)
**Model:** [`src/models/mographdrp_aligned.py`](./models/mographdrp_aligned.py) (frozen Phase 0 baseline, imported by the script)
**Raw results:** `src/final_model/results/ablation_results.csv`
**Companions:** [`results.md`](../docs/results.md), [`13_seed_variance_results.md`](../docs/13_seed_variance_results.md), [`09_split_protocol_comparison.md`](../docs/09_split_protocol_comparison.md)
**Date:** 2026-10-03 · **Runtime:** 6.5 h (24 runs, 23,523 s) · **Hardware:** one laptop GTX 1650 (4 GB), torch 2.13.0+cu130
**Status:** Phase 2 filled (components a-e, 3 seeds, grouped split). Phase 3 filled 2026-10-04 (feature gate for component f, 10 runs over 5 seeds, 98 min, commit `cc303ee`): gate not passed, see section 6f. The attention form of f was built and rejected on 2026-10-05 (section 6f). Graph decision filled 2026-10-12 (15 runs, 5 seeds, 5.8 h, commit `2e425a9`): all three graph-bearing configurations are eligible by the rule, see the update in section 5. The final model was chosen on 2026-10-12 and is named `final` (section 8). Filled 2026-10-13: `final` under the random split (laptop B, section 9) and the RF and MLP rows of the comparison table (section 10). Leave-drugs-out filled 2026-10-13 (laptop B, 10 runs, 3.7 h, section 11). Every run of the plan is now in this file. Interpretability filled 2026-10-13 in a minimal form (section 11): a sensitivity test of the direct-hit signal on `final`, three pre-written rules, all passed. Still open: the `TBD` cells of sections 3, 4 and 13 (justifications, hypotheses, figures).

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
| Seeds | 42, 43, 44 (3 per row); the Phase 3 gate and Phase 4 stage A rows use 42-46 (5 per arm). Since 2026-10-07 `base`, `base+proteomics`, `base+std_targets` and `base+proteomics+std_targets` also have 42-46 in the ladder CSV, and since 2026-10-12 so do the three graph-decision configurations (section 5) |
| Hardware | all 24 runs on one GTX 1650 ("laptop A"), one process, commit `66205f8`. Every ladder row in this file is from laptop A unless a table says "laptop B". Laptop B (a teammate's machine, since 2026-10-11) ran the cross-machine anchor, `HeteroIC50GNN`, the random split and leave-drugs-out; see "Second machine" below |
| Metrics | RMSE, MAE, R², PCC, SCC; AUC and F1 at the shared median threshold |

Early-stopping and LR-scheduler patience are both 200, so neither fires within
the 200-epoch budget. The model is the frozen `MoGraphDRPAligned` (BAN
bilinear head with 3 heads, (512, 128) predictor, gated-sum drug fusion).
Aligned to these settings on 2026-10-03; no ladder run predates that.

**What the seeds vary.** A seed changes weight initialisation and batch order
only. The grouped split is fixed at `random_state=42`, so every run uses the
same 80 test cell lines. The std below is training noise; variance from the
choice of test cell lines is not measured.

**A seed does not pin a run that uses the molecular graph (found
2026-10-13).** On laptop A, three 1-epoch repeats of `final` with seed 42
gave test RMSE 1.4190 / 1.4370 / 1.4338, while two repeats of
`base+proteomics+std_targets` gave 1.3344 both times. So configurations
without `mol_graph` reproduce exactly on one machine (as checked for `base`
on 2026-10-04 and 2026-10-07), and configurations with it (`base+mol_graph`,
`aligned`, `full`, the graph-decision rows, `final`) do not. The cause was
not traced; the GCN's sum pooling on the GPU is the likely source. For those
rows a "seed" is one draw, the per-seed pairing with `base` is a pairing of
draws rather than of identical initial conditions, and a re-run lands within
seed-to-seed noise instead of on the same number. No conclusion in this file
depends on exact reproduction of such a row.

**Second machine (laptop B), anchor run 2026-10-11, results committed in `dcbc304`.** GPU
results are not bit-identical across hardware, so no delta in this file is
taken between a laptop A row and a laptop B row. To measure how far the two
machines differ, laptop B ran `base` under the grouped split on seeds 42-46
(`results/laptop_b_anchor_results.csv`):

| Machine | Validation RMSE | Test RMSE | Per-seed test (42 / 43 / 44 / 45 / 46) |
|---|---|---|---|
| Laptop A | 1.2961 ± 0.0059 | 1.3252 ± 0.0119 | 1.3141 / 1.3347 / 1.3208 / 1.3406 / 1.3156 |
| Laptop B | 1.2928 ± 0.0151 | 1.3290 ± 0.0169 | 1.3220 / 1.3272 / 1.3053 / 1.3468 / 1.3435 |
| B − A | −0.0033 | +0.0038 | +0.008 / −0.008 / −0.016 / +0.006 / +0.028 |

- **The means agree to within 0.004**, far inside the ~0.03 band. Same data
  on both: 111,799 pairs, per-drug mean floor 1.4889.
- **The same seed does not give the same run** on the two machines: per-seed
  differences reach 0.028, as large as seed-to-seed noise. A seed-paired
  comparison is only valid within one machine.
- Laptop B is slower on this model, not faster: 871 s per `base` run against
  546 s (0.63×).

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
| f | Pair-specific component. Two tested forms: hand-built pair features (no effect), and attention of target proteins over mutated proteins (Phase 4 stage A, rejected as configured). Skipped, decided 2026-10-06: PPI message passing under that attention (stage B) | `pair_module=features` or `pair_module=attention` (alternatives, never combined) | adds 5 per-pair features, or a 64-wide attention vector, to the predictor input | this project (novelty component) | the base sees the cell line and the drug separately; nothing tells it whether this cell line's mutated proteins are this drug's targets |

Known caveats to state honestly in the report:

- **a:** each modality is one vector, so attention weights are always 1.0. As
  built, this rung measures the extra projection layers, not attention.
  Decided 2026-10-03: reported with this caveat, not fixed.
- **e:** not usable under leave-drugs-out (unseen drugs have no statistics).
- **f:** covers only pairs whose drug has a known target: 71.75% of pairs
  (80,213 of 111,799; test 12,119 of 16,840). 169 of the 240 drugs have a
  target and 71 have none. (The "130 drugs with no target" figure counts all
  498 drugs in the graph, most of which are not in these pairs.)
- **f, attention form:** the graph has no protein features (`protein.x` is
  all zeros), so each protein is an embedding learned from the response loss
  alone. Only 634 of the 16,214 proteins are a target or a mutation of any
  pair (148 targets, 536 mutated, 50 both). 15 of the 251 proteins mutated in
  test cell lines are mutated in no training cell line, so their embeddings
  are untrained at test time.

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
| `base+pair_features` | f (features only) | 1.3216 ± 0.0161 (5 seeds) | −0.0035 | no | 0.8804 | 0.7712 | +0.1673 | 3,186,817 | 592 |
| `base+pair_attention` | f (attention, no message passing) | 1.3861 ± 0.0142 (5 seeds) | +0.0609 | yes, worse | 0.8677 | 0.7483 | +0.1028 | 4,271,489 | 1,041 |
| `aligned` | b + c | 1.3416 ± 0.0088 | +0.0184 | no | 0.8751 | 0.7642 | +0.1473 | 3,528,562 | 1,254 |
| `full` | a-e | 1.3394 ± 0.0081 | +0.0162 | no | 0.8789 | 0.7650 | +0.1495 | 4,270,706 | 2,309 |

Negative Δ means better than base. Fit is the mean training time per run.
`base+pair_features` comes from the Phase 3 gate run
(`pair_gate_results.csv`) with 5 seeds (42-46), and `base+pair_attention`
from the Phase 4 stage A run in the same CSV and seeds. Their Δ is against
the gate's own 5-seed `base`, 1.3252 ± 0.0119; on seeds 42-44 that re-run
matched the `base` row above exactly. Subset results are in section 6f.

### Update 2026-10-07 — five seeds for the candidates, and their combination

Seeds 45 and 46 were added for `base`, `base+proteomics` and
`base+std_targets`, and the combination `base+proteomics+std_targets` was run
on seeds 42-46 (11 runs, same machine and CSV, commits `1c00b00` and
`70c3d13`, whose code is identical). Mean ± std over 5 seeds, grouped split.
The interval is a t-interval over the five per-seed differences (training
noise; the split is fixed):

| Config | Fold | RMSE | Δ vs base | Seeds better | 95% CI, seeds | Per-seed Δ (42 / 43 / 44 / 45 / 46) |
|---|---|---|---|---|---|---|
| `base` | val | 1.2961 ± 0.0059 | 0 | | | |
| `base` | test | 1.3252 ± 0.0119 | 0 | | | |
| `base+proteomics` | val | 1.2796 ± 0.0060 | −0.0166 | 5 / 5 | −0.023 to −0.010 | −0.020 / −0.014 / −0.009 / −0.020 / −0.021 |
| `base+proteomics` | test | 1.2920 ± 0.0125 | −0.0331 | 5 / 5 | −0.043 to −0.023 | −0.035 / −0.023 / −0.033 / −0.045 / −0.030 |
| `base+std_targets` | val | 1.2599 ± 0.0074 | −0.0362 | 5 / 5 | −0.042 to −0.030 | −0.032 / −0.037 / −0.042 / −0.040 / −0.030 |
| `base+std_targets` | test | 1.3086 ± 0.0144 | −0.0166 | 4 / 5 | −0.041 to +0.008 | −0.020 / −0.028 / +0.011 / −0.040 / −0.006 |
| **`base+proteomics+std_targets`** | val | **1.2454 ± 0.0037** | **−0.0508** | 5 / 5 | −0.058 to −0.044 | −0.055 / −0.052 / −0.043 / −0.057 / −0.048 |
| **`base+proteomics+std_targets`** | test | **1.2797 ± 0.0192** | **−0.0455** | 5 / 5 | −0.060 to −0.031 | −0.064 / −0.047 / −0.032 / −0.041 / −0.044 |

Test fold, other metrics (mean over 5 seeds):

| Config | PCC | R² | MAE | Gain over per-drug mean (1.4889) | Params | Fit (s) |
|---|---|---|---|---|---|---|
| `base` | 0.8793 | 0.7699 | 0.9979 | +0.1637 | 3,184,257 | 546 |
| `base+proteomics` | 0.8864 | 0.7813 | 0.9683 | +0.1969 | 3,316,225 | 598 |
| `base+std_targets` | 0.8859 | 0.7757 | 0.9763 | +0.1803 | 3,184,257 | 525 |
| `base+proteomics+std_targets` | 0.8908 | 0.7854 | 0.9560 | +0.2092 | 3,316,225 | 589 |

- **`base+proteomics+std_targets` is the first configuration that beats
  `base` beyond the ~0.03 band on both folds**, on every seed, with a seed
  interval that excludes 0 on both.
- **The combination against its parts** (paired by seed):

  | Compared with | Fold | Δ | Seeds better | 95% CI, seeds |
  |---|---|---|---|---|
  | `base+proteomics` | val | −0.0342 | 5 / 5 | −0.039 to −0.029 |
  | `base+proteomics` | test | −0.0124 | 3 / 5 | −0.030 to +0.006 |
  | `base+std_targets` | val | −0.0145 | 5 / 5 | −0.024 to −0.005 |
  | `base+std_targets` | test | −0.0289 | 5 / 5 | −0.052 to −0.006 |

  It is clearly better than either part on validation. On test it is clearly
  better than standardisation alone, but **not established as better than
  proteomics alone** (−0.012, 3 of 5 seeds, interval includes 0).
- **Standardisation's own test gain is still not established** (−0.017, 4 of
  5, interval includes 0), although on validation it is the larger of the two
  effects (−0.036). What it adds on test, it adds on top of proteomics only
  on some seeds.
- Best epoch: `base` 122 / 111 / 72 / 22 / 67; proteomics 162 / 180 / 102 /
  106 / 170; standardisation 14 / 20 / 36 / 38 / 96; the combination 43 / 28
  / 13 / 189 / 45.
- **Limits:** the intervals cover training noise only. The split is one fixed
  set of 80 validation and 80 test cell lines, and no cell-line bootstrap was
  run (the Phase 2 runs for seeds 42-44 saved no predictions). The
  combination was chosen after seeing the 3-seed validation results of its
  parts, with the keep rule written before any new seed existed (below,
  section 8).

### Update 2026-10-12 — graph decision: the candidate plus the built graph components

The project objectives need a graph component in the final model, and the
candidate has none (section 8). Three configurations add the built graph
components to it, on seeds 42-46: 15 runs, 2026-10-11 21:51 to 2026-10-12
03:40 (20,964 s, 5.8 h), same machine and CSV as every row above, commit
`2e425a9` (no model, data or ladder code changed since `70c3d13`). Grouped
split, floor 1.4889. Below, `C` is the candidate `base+proteomics+std_targets`.

**Rule, written on 2026-10-11 before the runs:** a graph-bearing
configuration is eligible if its validation RMSE is no more than 0.03 worse
than the candidate's (1.2454, so 1.2754 or better). Among eligible ones the
team picks the one that meets the project objective; if more than one does,
the one with the lower validation RMSE. Test is reported after and is not
used to choose.

| Config | Val RMSE | Δ val vs `C` | Seeds better | Eligible (≤ 1.2754) | Test RMSE | Δ test vs `C` | Seeds better | Δ test vs `base` |
|---|---|---|---|---|---|---|---|---|
| `C` (candidate, no graph) | 1.2454 ± 0.0037 | 0 | | | 1.2797 ± 0.0192 | 0 | | −0.0455 |
| `C+mol_graph` | 1.2559 ± 0.0039 | +0.0106 | 0 / 5 | **yes** | 1.2731 ± 0.0182 | −0.0066 | 4 / 5 | −0.0521 |
| `C+mol_graph+pair_features` | 1.2441 ± 0.0102 | −0.0013 | 3 / 5 | **yes** | 1.2737 ± 0.0144 | −0.0059 | 3 / 5 | −0.0514 |
| `C+mol_graph+pair_attention` | 1.2745 ± 0.0159 | +0.0291 | 0 / 5 | **yes, by 0.0009** | 1.3225 ± 0.0285 | +0.0428 | 0 / 5 | −0.0027 |

Intervals on the Δ against `C`. "Seeds" is a t-interval over the five
per-seed differences; "cell lines" is a paired bootstrap over the fold's 80
cell lines (`pair_gate_report.py --tag graph_decision`):

| Config | Fold | Δ vs `C` | 95% CI, seeds | 95% CI, cell lines | Per-seed Δ (42 / 43 / 44 / 45 / 46) |
|---|---|---|---|---|---|
| `C+mol_graph` | val | +0.0106 | +0.004 to +0.017 | +0.001 to +0.020 | +0.017 / +0.009 / +0.012 / +0.012 / +0.003 |
| `C+mol_graph` | test | −0.0066 | −0.019 to +0.006 | −0.018 to +0.005 | −0.009 / −0.005 / −0.016 / −0.013 / +0.010 |
| `C+mol_graph+pair_features` | val | −0.0013 | −0.010 to +0.007 | −0.011 to +0.009 | −0.008 / +0.009 / −0.006 / +0.003 / −0.003 |
| `C+mol_graph+pair_features` | test | −0.0059 | −0.026 to +0.015 | −0.017 to +0.006 | +0.006 / −0.019 / −0.024 / −0.009 / +0.015 |
| `C+mol_graph+pair_attention` | val | +0.0291 | +0.012 to +0.046 | +0.006 to +0.052 | +0.026 / +0.047 / +0.036 / +0.027 / +0.009 |
| `C+mol_graph+pair_attention` | test | +0.0428 | +0.014 to +0.072 | +0.021 to +0.064 | +0.028 / +0.026 / +0.039 / +0.037 / +0.083 |

Test fold, other metrics (mean over 5 seeds):

| Config | PCC | R² | MAE | Gain over per-drug mean (1.4889) | Params | Fit (s) | Best epoch (42-46) |
|---|---|---|---|---|---|---|---|
| `C` | 0.8908 | 0.7854 | 0.9560 | +0.2092 | 3,316,225 | 589 | 43 / 28 / 13 / 189 / 45 |
| `C+mol_graph` | 0.8919 | 0.7877 | 0.9502 | +0.2158 | 3,641,965 | 1,213 | 32 / 19 / 19 / 40 / 42 |
| `C+mol_graph+pair_features` | 0.8923 | 0.7874 | 0.9529 | +0.2152 | 3,644,525 | 1,221 | 36 / 22 / 29 / 18 / 33 |
| `C+mol_graph+pair_attention` | 0.8823 | 0.7708 | 0.9872 | +0.1664 | 4,729,197 | 1,753 | 12 / 8 / 6 / 8 / 17 |

- **All three are eligible by the rule.** `C+mol_graph` and
  `C+mol_graph+pair_features` pass with room. `C+mol_graph+pair_attention`
  passes by 0.0009 (1.2745 against 1.2754), which is far inside its own seed
  std (0.0159): two of its five seeds are above the line. Its eligibility is
  not a stable result.
- **By the rule's tie-break, the lower validation RMSE is
  `C+mol_graph+pair_features` (1.2441).** Which eligible configurations meet
  the project objective is the team's decision and has not been made; no
  configuration is named `final` yet.
- **The molecular graph still has no measurable effect.** On the candidate it
  is 0.011 worse on validation (all 5 seeds, both intervals exclude 0) and
  0.007 better on test (4 of 5 seeds, both intervals include 0). Both are
  inside the ~0.03 band, and the signs disagree between the folds, as on
  `base` (+0.009 test, −0.007 validation, 3 seeds). It doubles training time
  (1,213 s against 589 s).
- **Pair features on top of the molecular graph** (one change, against
  `C+mol_graph`): validation −0.0119 on all pairs (5 of 5 seeds; cell-line
  interval −0.023 to −0.001, seed interval −0.024 to +0.000); test +0.0007
  (2 of 5, cell-line interval −0.011 to +0.011). No effect on all-pairs test
  RMSE, as in Phase 3. The direct-hit subset is below.
- **Pair attention is worse on this foundation too.** Against `C+mol_graph`:
  validation +0.0186 (0 of 5), test +0.0494 (0 of 5, cell-line interval
  +0.025 to +0.075). On `base` it was +0.061 on test (section 6f). It peaks
  at epoch 6-17, earlier than on `base` (13-24). With it, the model is no
  better than `base` on test (−0.003, 3 of 5 seeds, seed interval −0.038 to
  +0.033): it removes the candidate's whole test gain. The attention check
  (`pair_attention_check.py`) was not run on these checkpoints.
- **Against `base`,** both of the first two keep the candidate's gain: test
  −0.0521 (`C+mol_graph`, seed interval −0.069 to −0.035) and −0.0514
  (`C+mol_graph+pair_features`, −0.068 to −0.034), 5 of 5 seeds each.

**Direct-hit pairs** (a target of the drug is itself mutated: 213 validation
pairs, 211 test; per-drug-mean floor 1.8246 / 1.9609). Signed error in
brackets, positive = predicted ln(IC50) too high:

| Config | Validation RMSE | Test RMSE |
|---|---|---|
| `C` | 1.5221 ± 0.0322 [+0.430] | 1.5788 ± 0.0764 [+0.370] |
| `C+mol_graph` | 1.5311 ± 0.0615 [+0.421] | 1.5826 ± 0.0595 [+0.382] |
| `C+mol_graph+pair_features` | 1.3056 ± 0.0952 [+0.113] | 1.3842 ± 0.0561 [−0.067] |
| `C+mol_graph+pair_attention` | 1.5296 ± 0.1171 [+0.305] | 1.6629 ± 0.0457 [+0.252] |

| `C+mol_graph+pair_features` − `C+mol_graph`, direct hits | Δ | Seeds better | 95% CI, cell lines | 95% CI, seeds |
|---|---|---|---|---|
| validation | −0.2254 | 5 / 5 | −0.364 to −0.077 | −0.339 to −0.112 |
| test | −0.1984 | 5 / 5 | −0.417 to +0.054 | −0.305 to −0.092 |

- **The candidate has the same direct-hit bias as `base`** (about +0.4 too
  high on both folds); proteomics, standardisation and the molecular graph do
  not touch it. The pair features remove most of it on every seed.
- **The direct-hit RMSE gain is now better on 5 of 5 seeds on both folds**
  (Phase 3, on `base`: 4 of 5). The validation intervals exclude 0. On test
  the seed interval excludes 0 and the cell-line interval does not. By the
  yardstick Phase 3 fixed for this subset (every seed on validation, and a
  test bootstrap interval that excludes 0) the first half is met and the
  second is not. That rule was written for `base+pair_features`; here it is
  a reading after the fact, not a pre-registered test.
- Nothing outside direct hits: 1 hop, 2+ hops and the no-target control all
  have intervals that include 0 on test.
- These subsets are 1.26% of pairs, so they cannot move all-pairs RMSE by
  more than about 0.005.

Source: `results/ablation_results_subsets_graph_decision.csv` and
`_paired_graph_decision.csv` (against `C`), and the same two files with the
suffix `_vs_mol_graph` (against `C+mol_graph`).

The tables below this point are the original Phase 2 results on 3 seeds, kept
as run.

**Headline (Phase 2, 3 seeds):** no component clearly beats `base`. Proteomics is the only one
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
| `base+pair_features` | 1.3015 | 1.3398 | 1.3230 |

Gate rows only, seeds 45 / 46: `base` 1.3406 / 1.3156, `base+pair_features`
1.3342 / 1.3098.

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
| f (features) | −0.0035 | 3 / 5 | −0.013 / +0.005 / +0.002 / −0.007 / −0.006 (seeds 42-46) | no effect on all pairs — **tested-and-rejected**; direct-hit subset in 6f |

**Validation and test disagree on the ranking.** Mean validation RMSE:

| Config | Val RMSE | Δ val vs base | Δ test vs base |
|---|---|---|---|
| `base+std_targets` | 1.2580 ± 0.0097 | −0.0368 | −0.0124 |
| `aligned` | 1.2703 ± 0.0037 | −0.0245 | +0.0184 |
| `base+proteomics` | 1.2808 ± 0.0066 | −0.0140 | −0.0302 |
| `base+mol_graph` | 1.2883 ± 0.0123 | −0.0065 | +0.0087 |
| `full` | 1.2904 ± 0.0072 | −0.0044 | +0.0162 |
| `base+pair_features` (5 seeds, vs the gate's 5-seed `base` at 1.2961 ± 0.0059) | 1.2901 ± 0.0052 | −0.0060 | −0.0035 |
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
- **On the candidate (2026-10-12, 5 seeds):** `base+proteomics+std_targets+mol_graph`
  against `base+proteomics+std_targets` is +0.0106 on validation (worse on 5
  of 5 seeds) and −0.0066 on test (better on 4 of 5, interval includes 0).
  Still no measurable effect, at 2.1× the training time. It is eligible for
  the final model under the graph-decision rule (section 5 update); if it is
  included, it is there for the project's scope, not for accuracy.

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

### d. Proteomics — kept: consistent on both folds, just past the band on test
- **Result (Phase 2, 3 seeds):** 1.2930 ± 0.0169, Δ −0.0302. Better than base on all 3 seeds
  (−0.035 / −0.023 / −0.033). The mean gain equals the ~0.03 band.
- **Interpretation:** the only component with a consistent gain. It is not a
  clear pass: one seed (43) is inside the band, and on validation the gain is
  smaller (−0.014).
- **Five seeds (2026-10-07; seeds 45 and 46 added, run from commit
  `1c00b00`, same machine and CSV):**

  | Fold | `base` | `base+proteomics` | Δ | Seeds better | 95% CI, seeds | Per-seed Δ (42 / 43 / 44 / 45 / 46) |
  |---|---|---|---|---|---|---|
  | validation | 1.2961 ± 0.0059 | 1.2796 ± 0.0060 | −0.0166 | 5 / 5 | −0.023 to −0.010 | −0.020 / −0.014 / −0.009 / −0.020 / −0.021 |
  | test | 1.3252 ± 0.0119 | 1.2920 ± 0.0125 | −0.0331 | 5 / 5 | −0.043 to −0.023 | −0.035 / −0.023 / −0.033 / −0.045 / −0.030 |

  Better than `base` on every seed on both folds; the seed interval excludes
  0 on both. On test the mean gain is just past the ~0.03 band (PCC 0.8864
  against 0.8793, R² 0.7813 against 0.7699, gain over the per-drug mean floor
  +0.1969 against +0.1637). On validation it is about half that and inside
  the band. The interval covers training noise only: the split is fixed, and
  no cell-line bootstrap was possible because the Phase 2 runs (seeds 42-44)
  saved no predictions. `base` on seeds 45 and 46 reproduced the Phase 3 gate
  rows exactly.
- **Why (evidence):** not diagnosed beyond the numbers above.
- **Cost:** +131,968 parameters, 1.07× training time.
- Earlier flat-model result to reconcile: Mut_CNV hurt (E32 vs E04 in [results.md](../docs/results.md)). Not tested here: there is no rung without Mut_CNV.

### e. Per-drug target standardisation — kept by the validation rule; no established effect on test alone
- **Five seeds (2026-10-07):** validation 1.2599 ± 0.0074, Δ −0.0362, better
  on 5 of 5 seeds (seed interval −0.042 to −0.030), beyond the band. Test
  1.3086 ± 0.0144, Δ −0.0166, better on 4 of 5 (interval −0.041 to +0.008),
  inside the band. The verdict changed from the 3-seed "tested-and-rejected"
  because the choice of components is made on validation (section 8); on
  test alone its effect is still not established. Combined with proteomics
  it gives the best configuration so far (section 5 update).
- **Result (Phase 2, 3 seeds):** 1.3108 ± 0.0195, Δ −0.0124. Inside the noise band; better than
  base on 2 seeds of 3, worse on seed 44 (+0.011).
- **Interpretation:** no measurable effect on test RMSE for this model. It has
  the best validation RMSE of any rung (1.2580, −0.037 vs base), but that gain
  does not carry to the test cell lines.
- **Why (evidence):** best epoch was 14-36, among the earliest of the single
  components, so it fits fast and then overfits. Not diagnosed further.
- **Cost:** no extra parameters, same training time (0.97×).
- Earlier single-seed result on HeteroIC50GNN: 1.3449 → 1.2917 (−0.053). That gain does not reproduce on this model over 3 seeds.

### f. Pair-specific component — features: no effect (Phase 3); attention: tested-and-rejected as configured (Phase 4 stage A); message passing: skipped (stage B)

**Gate, 2026-10-04: not passed.** `base+pair_features` appends five hand-built
features to the predictor input (drug has a target, cell line has a mutation,
a target is itself mutated, 1 / (1 + min PPI hops), log count of mutated
proteins within one hop of a target). It is not a GNN. Run beside a re-run of
`base` in the gate's own CSV, commit `cc303ee`: seeds 42-44 first, then seeds
45-46 added the same day to settle the direct-hit result, so **5 seeds per
arm**. The `base` re-run reproduced the Phase 2 rows exactly on seeds 42-44
(same validation RMSE, test RMSE and best epoch).

- **Result (all pairs):** 1.3216 ± 0.0161 against the gate's 5-seed `base` at
  1.3252 ± 0.0119: Δ −0.0035, better on 3 seeds of 5. No effect —
  **tested-and-rejected** on the headline metric.
- **Coverage:** 71.75% of pairs have a drug with a known target. A target is
  directly mutated in only 1.26% of pairs (test: 211 of 16,840), so all-pairs
  RMSE could not have moved by more than about 0.005 even with perfect use of
  that fact. The gate was therefore read on subsets, with the rule fixed
  before the first run: direct-hit RMSE better on every seed on validation,
  and a test bootstrap interval that excludes 0.

Test fold, mean ± std over 5 seeds. Two intervals for Δ: a paired bootstrap
over test cell lines (2,000 resamples, seeds averaged), and a t-interval over
the five per-seed differences (training noise):

| Subset | Pairs | Per-drug mean | `base` | `base+pair_features` | Δ | Seeds better | 95% CI, cell lines | 95% CI, seeds |
|---|---|---|---|---|---|---|---|---|
| all | 16,840 | 1.4889 | 1.3252 ± 0.0119 | 1.3216 ± 0.0161 | −0.0035 | 3 / 5 | −0.016 to +0.009 | −0.012 to +0.005 |
| drug has a target | 12,119 | 1.4494 | 1.2632 ± 0.0163 | 1.2568 ± 0.0170 | −0.0064 | 4 / 5 | −0.022 to +0.009 | −0.020 to +0.007 |
| drug has no target (control) | 4,721 | 1.5857 | 1.4722 ± 0.0098 | 1.4751 ± 0.0167 | +0.0028 | 3 / 5 | −0.011 to +0.018 | −0.014 to +0.020 |
| **target directly mutated** | 211 | 1.9609 | 1.7157 ± 0.1212 | 1.5771 ± 0.1528 | −0.1386 | 4 / 5 | −0.304 to +0.067 | −0.374 to +0.097 |
| nearest mutation 1 hop away | 8,688 | 1.4699 | 1.2924 ± 0.0149 | 1.2899 ± 0.0164 | −0.0025 | 3 / 5 | −0.020 to +0.014 | −0.010 to +0.005 |
| nearest mutation 2+ hops away | 3,051 | 1.3578 | 1.1455 ± 0.0147 | 1.1421 ± 0.0120 | −0.0034 | 3 / 5 | −0.017 to +0.010 | −0.019 to +0.013 |

Direct-hit pairs, RMSE per seed (42 / 43 / 44 / 45 / 46):

| Fold | Pairs | `base` | `base+pair_features` | Δ per seed | Mean Δ | 95% CI, cell lines | 95% CI, seeds |
|---|---|---|---|---|---|---|---|
| test | 211 | 1.619 / 1.676 / 1.716 / 1.923 / 1.645 | 1.496 / 1.838 / 1.539 / 1.565 / 1.446 | −0.123 / +0.163 / −0.177 / −0.357 / −0.199 | −0.139 | −0.304 to +0.067 | −0.374 to +0.097 |
| val | 213 | 1.427 / 1.542 / 1.758 / 1.894 / 1.446 | 1.276 / 1.591 / 1.366 / 1.362 / 1.293 | −0.151 / +0.049 / −0.392 / −0.532 / −0.153 | −0.236 | −0.357 to −0.094 | −0.518 to +0.047 |

Mean signed error (prediction − measured) on direct-hit pairs, per seed:

| Fold | Config | Seed 42 | Seed 43 | Seed 44 | Seed 45 | Seed 46 | Mean |
|---|---|---|---|---|---|---|---|
| test | `base` | +0.275 | +0.367 | +0.527 | +0.609 | +0.244 | +0.404 |
| test | `base+pair_features` | +0.056 | +0.123 | −0.050 | −0.002 | +0.016 | +0.028 |
| val | `base` | +0.330 | +0.431 | +0.540 | +0.645 | +0.327 | +0.454 |
| val | `base+pair_features` | +0.125 | +0.112 | +0.090 | +0.074 | +0.082 | +0.097 |

- **Interpretation:**
  - `base` is biased on direct-hit pairs: it predicts ln(IC50) about 0.4 too
    high (too resistant) on every seed and both folds. The Mut_CNV PCA input
    does not give it this.
  - The features remove most of that bias on all 5 seeds and both folds.
  - The RMSE gain on those pairs is **likely but not established**. It is
    better on 4 seeds of 5 on both folds; seed 43 is worse on both. Mean Δ is
    −0.139 on test and −0.236 on validation. On validation the cell-line
    interval excludes 0; on test both intervals include 0, and on validation
    the seed interval does too. The gate rule is not met.
  - With 3 seeds the test Δ was −0.046 (2 of 3 seeds); the two added seeds
    both improved (−0.357, −0.199).
  - Nothing at one hop or more (−0.0025, −0.0034), no different from the
    no-target control (+0.0028), where the features are constant. The
    has-target Δ (−0.0064) is inside both intervals. Hand-built PPI proximity
    adds nothing; whatever signal there is sits in "a target is mutated".
- **Limits of this test:** 5 seeds; 211 test pairs; one fixed split. The hop
  features are gene-agnostic counts: they rule out proximity mattering in
  general, not a specific mutated neighbour mattering for a specific drug,
  which only a learned module could test. The direct-hit effect itself
  differs by drug (from −1.8 to +0.4 z on the training rows), which a single
  flag cannot express.
  Validation and test subsets were read in the same report. Seeds 45-46 were
  added after seeing the 3-seed result, to settle a 2-of-3 split, not planned
  in advance.
- **Cost:** +2,560 parameters, 1.03× training time (592 s vs 576 s).
- **Leakage:** mutation and target edges are inputs, not labels. A held-out
  cell line's mutation edges describe it in the same way its expression
  profile does, so their presence in the graph is not leakage.
- **Decision (2026-10-05):** Phase 4 is built in two stages. First
  attention between target and mutated proteins with no message passing,
  judged against these features on direct-hit pairs. Then 1 and 2 PPI
  message-passing layers as the depth ablation below. Plan in `CLAUDE.md`.
- **Decision (2026-10-06):** the first stage failed (below), and the second
  stage (message passing) is skipped for now. It was not built or run.

**Stage A, `base+pair_attention`: run 2026-10-05, failed its pass rule —
tested-and-rejected as configured.** 5 seeds (42-46), commit `9c7ba4f`, in
the gate CSV beside the Phase 3 rows, which were reused, not re-run.

- **Result (all pairs):** 1.3861 ± 0.0142 against `base` at 1.3252 ± 0.0119:
  Δ **+0.0609**, worse on 5 seeds of 5, beyond the ~0.03 band (cell-line
  interval +0.026 to +0.095, seed interval +0.032 to +0.089). Validation:
  1.3481 ± 0.0269 against 1.2961 ± 0.0059 (+0.0519, worse on 5 of 5).
  Gain over the per-drug mean floor (1.4889): +0.1028, against +0.1637 for
  `base`.
- **Pass rule, both conditions failed:**
  - Direct-hit pairs, validation, against `base+pair_features`: Δ +0.131,
    better on 1 seed of 5, cell-line interval +0.050 to +0.212. Worse, not
    better.
  - All pairs against `base`: +0.061, outside the band.

Test fold, mean ± std over 5 seeds, against `base`. Intervals as in the
Phase 3 table:

| Subset | Pairs | Per-drug mean | `base` | `base+pair_attention` | Δ | Seeds better | 95% CI, cell lines | 95% CI, seeds |
|---|---|---|---|---|---|---|---|---|
| all | 16,840 | 1.4889 | 1.3252 ± 0.0119 | 1.3861 ± 0.0142 | +0.0609 | 0 / 5 | +0.026 to +0.095 | +0.032 to +0.089 |
| drug has a target | 12,119 | 1.4494 | 1.2632 ± 0.0163 | 1.3327 ± 0.0214 | +0.0695 | 0 / 5 | +0.035 to +0.104 | +0.027 to +0.112 |
| drug has no target (control) | 4,721 | 1.5857 | 1.4722 ± 0.0098 | 1.5143 ± 0.0120 | +0.0421 | 0 / 5 | −0.001 to +0.088 | +0.034 to +0.050 |
| **target directly mutated** | 211 | 1.9609 | 1.7157 ± 0.1212 | 1.6164 ± 0.1369 | −0.0993 | 3 / 5 | −0.306 to +0.104 | −0.336 to +0.137 |
| nearest mutation 1 hop away | 8,688 | 1.4699 | 1.2924 ± 0.0149 | 1.3630 ± 0.0343 | +0.0706 | 0 / 5 | +0.037 to +0.102 | +0.013 to +0.129 |
| nearest mutation 2+ hops away | 3,051 | 1.3578 | 1.1455 ± 0.0147 | 1.1997 ± 0.0078 | +0.0541 | 0 / 5 | +0.015 to +0.092 | +0.032 to +0.076 |

Direct-hit pairs, RMSE per seed (42 / 43 / 44 / 45 / 46), attention against
the features:

| Fold | Pairs | `base+pair_features` | `base+pair_attention` | Δ per seed | Mean Δ | 95% CI, cell lines | 95% CI, seeds |
|---|---|---|---|---|---|---|---|
| val (decides) | 213 | 1.276 / 1.591 / 1.366 / 1.362 / 1.293 | 1.562 / 1.711 / 1.291 / 1.551 / 1.428 | +0.286 / +0.121 / −0.076 / +0.189 / +0.135 | +0.131 | +0.050 to +0.212 | −0.034 to +0.296 |
| test | 211 | 1.496 / 1.838 / 1.539 / 1.565 / 1.446 | 1.763 / 1.704 / 1.417 / 1.650 / 1.547 | +0.267 / −0.134 / −0.122 / +0.085 / +0.100 | +0.039 | −0.119 to +0.236 | −0.170 to +0.249 |

Mean signed error (prediction − measured) on direct-hit pairs, mean over 5
seeds (per-seed values in `pair_gate_results_subsets_attention_vs_base.csv`):

| Fold | `base` | `base+pair_features` | `base+pair_attention` |
|---|---|---|---|
| test | +0.404 | +0.028 | +0.160 |
| val | +0.454 | +0.097 | +0.294 |

Did the attention find the mutated target? For each target that is itself
mutated in the cell line, the attention it gives to the key that is the same
protein, mean over the 4 heads (`pair_attention_check.py`; trained = the 5
saved checkpoints, untrained = 20 random initialisations):

| Fold | Chance (1 / keys) | Untrained | Trained | Matching key ranked first, untrained | Trained |
|---|---|---|---|---|---|
| train (989 pairs) | 0.127 | 0.127 ± 0.010 | 0.145 ± 0.015 | 0.123 ± 0.076 | 0.214 ± 0.068 |
| val (213) | 0.126 | 0.125 ± 0.012 | 0.147 ± 0.017 | 0.110 ± 0.091 | 0.238 ± 0.095 |
| test (211) | 0.148 | 0.148 ± 0.013 | 0.167 ± 0.008 | 0.146 ± 0.090 | 0.208 ± 0.055 |

- **Interpretation:**
  - The rung makes the model worse everywhere except on direct-hit pairs,
    including on pairs whose drug has no target (+0.042, 0 of 5 seeds
    better), where the query is one constant token. The loss is therefore
    not about target-mutation matching.
  - It overfits early: best epoch 24 / 13 / 24 / 15 / 15, against
    122 / 111 / 72 / 22 / 67 for `base`. Validation RMSE varies more between
    seeds (std 0.027 against 0.006).
  - On direct-hit pairs it is between `base` and the features: it removes
    part of the base's bias (+0.40 to +0.16 on test) where the features
    remove nearly all of it (+0.03). Against `base` its RMSE there is −0.099
    on test (3 of 5 seeds, both intervals include 0) and −0.105 on
    validation (3 of 5; cell-line interval −0.186 to −0.005, seed interval
    includes 0).
  - **The attention did not learn the match.** A mutated target gives the
    matching key 0.145-0.167 of its attention, about 0.02 above chance and
    above an untrained module, and ranks it first on 21-24% of occasions.
    The checkpoint is taken at the best validation epoch, which came after
    13-24 epochs.
- **What this does and does not show:** this module, trained this way, is
  rejected. It is not evidence about pair-specific attention as an idea:
  the attention never found the mutated target, so the idea was not put to
  the test. The five hand-built features remain the better form of
  component f on every subset.
- **Why it did not learn, and why it hurts** (diagnosed 2026-10-05 from the
  5 checkpoints and the saved predictions; ad-hoc analysis, no committed
  script yet):
  - The attention stayed close to uniform: its entropy is 0.963 of the
    maximum, against 0.987 for untrained modules.
  - The embeddings barely moved from their random start. Mean row norm is
    7.96 for the proteins in use and 7.97 for rows that never receive a
    gradient. They start at unit variance, the learning rate is 1e-4, and
    the checkpoint is taken after 13-24 epochs.
  - The signal that would teach the match is 1.27% of training pairs, about
    1.6 pairs per batch of 128.
  - So the pair vector is close to an average of fixed random vectors for
    the cell line's mutated proteins. Over all 531 × 240 (cell line, drug)
    combinations, 54% of its variance depends on the cell line alone, 31% on
    the drug alone and 16% on the pair (untrained: 84% / 5% / 11%). It acts
    as a random fingerprint of the cell line, not as a pair-specific signal.
  - The predictor fits that fingerprint, which cannot transfer to held-out
    cell lines: 61 of the 80 test cell lines are worse than under `base`,
    more so for cell lines with few mutations (Spearman −0.26 between the
    loss and the mutation count, p 0.018). Validation then peaks early and
    the checkpoint is taken before the attention has learned anything.
  - That the predictor is fitting training cell lines is inferred from the
    points above; training-fold predictions were not saved, so it was not
    measured directly.
- **Limits of this test:** 5 seeds, one fixed split, 211-213 direct-hit
  pairs; no tuning of the embedding width, dropout or weight decay (attention
  dropout 0, weight decay 0, as the frozen base); validation and test subsets
  were read in the same report, with the rule read on validation.

What the rung is:

- **What it adds:** each protein is a 64-wide learned embedding. A drug's
  target proteins (up to 7) are the queries and the cell line's mutated
  proteins (up to 108, median 6) the keys and values of one 4-head attention
  layer. The attended targets are mean-pooled to one 64-wide vector per pair,
  concatenated with the base's interaction vector before the predictor. No
  message passing, so it is not a GNN.
- **Empty sets:** a drug with no target (71 of 240) queries with a learned
  "no target" token. A learned "no mutation" token is a key for every cell
  line, not only the 6 with no mutation: attention weights sum to 1, so a
  target with no relevant mutation needs somewhere to attend.
- **What it can add over the features:** which protein is mutated. Drug
  specificity alone is not new: the predictor already sees the direct-hit
  flag beside the drug vector.
- **Not built in:** nothing tells the module that a target and a mutated
  protein are the same protein. Both read the same embedding table, so it can
  learn that, from 989 direct-hit training pairs. It did not (table above).
- **Direct hits are concentrated:** 85 of the 240 drugs have any; drug `1931`
  has 372 of the 1,413; the median such drug has 6 in the training fold.
- **Cost:** +1,087,232 parameters (4,271,489 total): 1,037,824 in the
  embedding table, 16,640 in the attention layer, 32,768 in the widened
  predictor. At most 636 embedding rows (40,704 values) can receive a
  gradient; the rest of the table is unused (stage B, which would have used
  it, was skipped). Training time
  1.81× `base` (1,041 s against 576 s per run, 5 runs in 87 min on the
  GTX 1650); peak GPU memory 227 MiB.
- **Pass rule, fixed before the run:** on direct-hit pairs, better than
  `base+pair_features` on at least 4 of 5 seeds on validation, with the
  validation cell-line interval excluding 0; and no worse than `base` on all
  pairs beyond the ~0.03 band. Test is reported after, not used to decide.
- **Checks passed before the run:** 1-epoch CPU `base` (val 2.688133 / test
  2.747617) and `base+pair_features` (2.554799 / 2.649105) are unchanged by
  the new code; the Phase 3 feature matrix is byte-identical; the padded
  target and mutation sets reproduce the Phase 3 direct-hit flag on all
  111,799 rows; the report regenerates the two Phase 3 CSVs exactly.

- **GNN depth ablation** (does message passing matter?): the 0-layer row is
  Stage A. **The 1- and 2-layer rows (Stage B) were skipped, decided
  2026-10-06.** Not built, not run. Reasons: message passing would sit under
  an attention that did not learn, so a null result would say little about
  the PPI graph; and the 10 runs were estimated at about 8.5 h of GPU time.
  It is reopened only if the supervisor requires a measured message-passing
  rung.
- **What this leaves unanswered:** no trained rung in this ladder uses the
  STRING network with message passing, so the depth ablation is a single
  row. The evidence on the PPI graph is indirect: the Phase 3 hop features
  (no effect at one hop or more, 5 seeds) and `HeteroIC50GNN` as a compared
  baseline (section 10). This does not show that message passing fails; it
  was not tested here.

| Message-passing layers | RMSE (all pairs) | RMSE (pairs with a known target) |
|---|---|---|
| 0 (attention on raw embeddings, no GNN) | 1.3861 ± 0.0142 | 1.3327 ± 0.0214 |
| 1 | not run (skipped) | not run (skipped) |
| 2 | not run (skipped) | not run (skipped) |

For reference on the same rows: `base` 1.3252 ± 0.0119 / 1.2632 ± 0.0163,
`base+pair_features` 1.3216 ± 0.0161 / 1.2568 ± 0.0170.

Source: `src/final_model/results/pair_gate_results.csv`. Phase 3:
`pair_gate_results_subsets.csv`, `pair_gate_results_paired.csv`. Stage A:
the same two names with `_attention_vs_base` and `_attention_vs_features`,
and `pair_gate_results_attention_check.csv`. Produced by
`src/final_model/pair_gate_report.py` and
`src/final_model/pair_attention_check.py` (commands in section 15).

**Both forms on top of the candidate (2026-10-12, 5 seeds).** Added to
`base+proteomics+std_targets+mol_graph`, the two forms repeat their results
on `base`. The features do nothing to all-pairs test RMSE (+0.0007) and cut
direct-hit RMSE by 0.20-0.23, now on 5 of 5 seeds on both folds, with the
test cell-line interval still including 0. The attention is worse on all
pairs (+0.049 on test, 0 of 5 seeds) and on test direct hits (+0.080, 0 of
5; validation −0.002). Tables in
the section 5 update of 2026-10-12.

## 7. Do the components combine?

*(report)* Components can overlap or interfere, so the sum of single gains
rarely equals the gain of the full model.

| Quantity | RMSE change vs base |
|---|---|
| Sum of the individual deltas (a + b + c + d + e) | +0.5685 |
| Same sum without c (bilinear), whose single-rung result is unstable | +0.0574 |
| `full` measured | +0.0162 |
| Difference, `full` − sum of all five | −0.5523 |

The two kept components, 5 seeds (2026-10-07):

| Quantity | Validation | Test |
|---|---|---|
| d alone (`base+proteomics`) | −0.0166 | −0.0331 |
| e alone (`base+std_targets`) | −0.0362 | −0.0166 |
| Sum of the two | −0.0528 | −0.0497 |
| `base+proteomics+std_targets` measured | −0.0508 | −0.0455 |
| Measured as a share of the sum | 96% | 92% |

- **Proteomics and standardisation are close to additive** on both folds:
  one changes the inputs, the other the target scale. This is the opposite of
  `full`, where the same two sit beside three components that hurt or do
  nothing, and the result is worse than `base`.

The original five-component analysis (Phase 2, 3 seeds):

- The components are **not additive**. The sum is dominated by bilinear's
  +0.51, which does not appear inside `full`: with the molecular graph
  present, the head trains normally (`aligned` +0.0184, `full` +0.0162).
- Does `full` beat the best single component beyond noise? **No.** `full`
  (1.3394) is 0.046 worse than `base+proteomics` (1.2930) and is not better
  than `base` on any seed.
- `full` overfits early: best epoch 47 / 21 / 10.
- `full` here means base + a-e. It is pinned to those five in the script, so
  the pair component (f) is not part of it.
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

**Chosen 2026-10-12: `final` =
`base+proteomics+std_targets+mol_graph+pair_features`.** It is named `final`
in `run_ablation.py`. Grouped split, 5 seeds (42-46), laptop A:

| | Validation RMSE | Test RMSE | PCC | R² | MAE | Gain over per-drug mean (1.4889) |
|---|---|---|---|---|---|---|
| `final` | 1.2441 ± 0.0102 | **1.2737 ± 0.0144** | 0.8923 | 0.7874 | 0.9529 | +0.2152 |
| `base` | 1.2961 ± 0.0059 | 1.3252 ± 0.0119 | 0.8793 | 0.7699 | 0.9979 | +0.1637 |
| Δ, paired by seed | −0.0521, 5 of 5 | −0.0514, 5 of 5 | | | | |

Per-seed test RMSE: 1.2567 / 1.2694 / 1.2653 / 1.2906 / 1.2867. Best epoch
36 / 22 / 29 / 18 / 33. 3,644,525 parameters, about 20 min per run.

- **How it was chosen.** Two steps, both on validation RMSE with test
  reported after. First the keep rule of 2026-10-06 gave
  `base+proteomics+std_targets` (below). Then, because the project requires
  a graph component, the graph-decision rule of 2026-10-11 was applied to
  that candidate plus the built graph components (section 5, update of
  2026-10-12): a config was eligible if its validation RMSE was no more than
  0.03 worse than the candidate's. All three were eligible; this one has the
  lowest validation RMSE.
- **Where its gain comes from.** Proteomics and standardisation. The two
  graph-related parts do not move all-pairs RMSE beyond noise: the same
  model without them is 1.2797 ± 0.0192 on test (−0.006, 3 of 5 seeds in
  favour of `final`).
- **What the graph-related parts are.** `mol_graph` is a GCN over each
  drug's atoms and bonds, taken from MoGraphDRP; it is the model's only
  graph neural network. `pair_features` are five precomputed features from
  the STRING PPI, drug-target and driver-mutation edges; they are not
  message passing. They remove the base's bias on the 1.26% of pairs whose
  target is mutated, on 5 of 5 seeds.
- **What it is not.** Not a heterogeneous GNN, not a GAT, no cross-attention
  fusion. Those were in the proposal; the first is covered only by the
  compared `HeteroIC50GNN` baseline (section 10), and cross-attention fusion
  was tested and rejected (section 6a).
- **Under the random split** (laptop B, 3 seeds, section 9): 0.9622 ± 0.0058
  against `base` 1.0162 ± 0.0011 (−0.054, 3 of 3 seeds). Not the best
  configuration there: `base+std_targets` alone is 0.020 better on every
  seed.
- **On unseen drugs** (leave-drugs-out, laptop B, 5 seeds, section 11):
  run without standardisation, which cannot be used on unseen drugs, it
  reaches 2.0038 ± 0.0843 against `base` 2.1722 ± 0.0369 (−0.169, 4 of 5
  seeds) and a mean-only floor of 2.5500. The interval over the 36 held-out
  drugs includes 0, and it is not better than the earlier flat MLP's single
  run (1.8516). The model is built and selected for unseen cell lines; it
  makes no claim about new compounds.
- **What it does with the pair features** (section 11, 2026-10-13): with
  the direct hit removed at inference its RMSE on those pairs goes from 1.41
  to 1.64 on test, and its effect per drug tracks the measured one
  (Spearman +0.73 over 30 drugs; 0 of 50 untrained models reach it).
- **As served by the application (2026-10-26, `fbf49b3`):** the app averages
  the five `--save-weights` checkpoints of `final` (seeds 42-46,
  `results/final_weights_results.csv`) as an ensemble; nothing is retrained.
  Recomputed 2026-10-26 on laptop A from the app's own session code:
  validation RMSE 1.2118 (PCC 0.9072) against 1.2349-1.2594 for the single
  seeds, test RMSE 1.2377 (PCC 0.8977) against 1.2482-1.2803 (their mean
  1.2698). Each seed reproduces its recorded validation RMSE when the app
  loads it. The ensemble was not part of the selection rule: it is the
  chosen model averaged over its seeds, reported after the choice. The
  per-seed spread is the app's confidence score; it is not calibrated.
- **Nothing is left to run for it.**

The record of how the candidate was reached follows.

- **Candidate (2026-10-07): `base+proteomics+std_targets`.** It passes the
  rule below. Whether it is the final model is the team's decision, not yet
  made; it has no name of its own in `run_ablation.py` yet.
  - **Rule, written on 2026-10-06 before any of the new seeds existed:** a
    component is kept if its rung beats `base` on validation on at least 4 of
    5 seeds. The combined model is accepted if its validation RMSE is no
    worse than the best kept single rung. Test is reported afterwards, not
    used to decide.
  - **Applied:** proteomics is kept (5 of 5 on validation, −0.0166);
    standardisation is kept (5 of 5, −0.0362). The combination's validation
    RMSE is 1.2454 ± 0.0037, better than the best single rung (1.2599), so it
    is accepted.
  - **Test, reported after:** 1.2797 ± 0.0192 against `base` 1.3252 ± 0.0119
    (Δ −0.0455, 5 of 5 seeds, seed interval −0.060 to −0.031), PCC 0.8908,
    R² 0.7854, gain over the per-drug mean floor +0.2092 (`base` +0.1637).
  - **What it does not show:** on test it is not established as better than
    `base+proteomics` alone (−0.012, 3 of 5 seeds). The rule was written
    after seeing the 3-seed validation results of both parts. One fixed
    split; intervals cover training noise only.
  - **Constraint:** standardisation cannot be used under leave-drugs-out, so
    that run uses `base+proteomics`.
- **Earlier position (Phase 2, 3 seeds):** not chosen; the evidence then was:
  - best on test: `base+proteomics`, 1.2930 ± 0.0169 (PCC 0.8871, R² 0.7810);
  - best on validation: `base+std_targets`, 1.2580 ± 0.0097 (test 1.3108);
  - `full` (base + a-e) is not better than `base` and should not be presented
    as the final model on these numbers.
- **Components dropped and why:** a (+0.0912, hurts), b (+0.0087, no effect),
  c (+0.5111, hurts and unstable as configured), f as hand-built features
  (−0.0035, no effect on all pairs), f as attention without message passing
  (+0.0609, hurts as configured). **Kept:** d (proteomics) and e
  (standardisation), see above.
- **Not tested:** f with PPI message passing (stage B, skipped 2026-10-06).
  The final model therefore cannot claim a graph component from this ladder.
- **Status (2026-10-07):** a candidate exists and passes the rule; confirming
  it as final and naming it is the next step, and unblocks the full ablation
  (Phase 5 in `CLAUDE.md`). A further candidate not yet run is a growth-rate
  covariate (section 11).
- **Status (2026-10-11): the candidate is not the final model as it stands.**
  The project objectives require a graph component and the candidate has
  none. Not yet run: the candidate plus `mol_graph`, plus
  `mol_graph+pair_features`, plus `mol_graph+pair_attention`, 5 seeds each.
  Rule written before the runs: a graph-bearing configuration is eligible if
  its validation RMSE is no more than 0.03 worse than the candidate's
  (1.2454); among eligible ones the team picks the one that meets the
  objective; test is reported after. What is already known about the parts,
  on `base`: `mol_graph` +0.009 test and −0.007 validation (3 seeds, no
  measurable effect), `pair_features` −0.004 (5 seeds), `pair_attention`
  +0.061 (5 seeds). A graph component kept on these numbers is kept for the
  project's scope, not for accuracy, and is reported that way. Commands and
  the split of the remaining runs across two laptops:
  [`docs/18_final_runs_runbook.md`](../docs/18_final_runs_runbook.md).
- **Status (2026-10-12): the graph decision is run; the choice is not yet
  made.** All three graph-bearing configurations are eligible by the rule
  (validation 1.2754 or better): `+mol_graph` 1.2559, `+mol_graph+pair_features`
  1.2441, `+mol_graph+pair_attention` 1.2745. The third passes by 0.0009,
  inside its seed std of 0.0159, and is 0.043 worse than the candidate on
  test on every seed. Still the team's call: which of them meet the project
  objective. If both of the first two do, the rule's tie-break (lower
  validation RMSE) gives `base+proteomics+std_targets+mol_graph+pair_features`.
  Neither is established as better or worse than the candidate on test
  (−0.007 and −0.006, intervals include 0), so a graph component is included
  for the project's scope and reported as having no measurable effect on
  all-pairs RMSE. What the pair features do change is the 1.26% of pairs
  where a target is itself mutated (section 5 update). Note for the write-up:
  `mol_graph` is a GCN over the drug molecule; the pair features are 5
  numbers computed from the PPI, drug-target and mutation graph, with no
  trained message passing. No configuration is named `final` in
  `run_ablation.py` yet.
- **Final test RMSE / PCC / R²:** 1.2731 ± 0.0182 / 0.8919 / 0.7877 if
  `+mol_graph` is chosen; 1.2737 ± 0.0144 / 0.8923 / 0.7874 if
  `+mol_graph+pair_features` is (candidate without a graph: 1.2797 ± 0.0192 /
  0.8908 / 0.7854).
- **Gain over the per-drug mean floor:** every stable rung beats 1.4889. On 5
  seeds: `+mol_graph` +0.2158, `+mol_graph+pair_features` +0.2152, the
  candidate +0.2092, `base+proteomics` +0.1969, `+mol_graph+pair_attention`
  +0.1664, `base` +0.1637.

Selection is made on **validation** RMSE; test numbers are reported after the
choice, not used to make it. Validation and test rank the rungs differently
(section 5): standardisation is the larger effect on validation and
proteomics on test.

## 9. Sanity check against the benchmark — random split

Not a ranking. It checks that the reproduction (`aligned`) is faithful.

| Config | Random-split RMSE | Reference |
|---|---|---|
| `aligned` (Phase 0, GE+Proteomics) | **0.978 ± 0.002** (3 seeds) | MoGraphDRP without XGBoost: 0.9497 |
| `base` | 1.0162 ± 0.0011 (3 seeds, laptop B) | — |
| `final` | 0.9622 ± 0.0058 (3 seeds, laptop B) | MoGraphDRP with XGBoost: 0.6622. A reference, not a ranking: different pairs and data |
| `full` | not run under this split; `final` took its place | — |

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

**Ladder configs under the random split (laptop B, seeds 42-44,
`results/random_split_results.csv`).** The four single configs were run on
2026-10-11 (committed in `e34f390`) and `final` on 2026-10-12 (committed in
`c346fc9`), all on laptop B, 15 rows. A sanity check, never ranked against
the grouped rows. Floors on this split: global mean 2.8384, per-drug mean
1.5265.

| Config | Validation RMSE | Test RMSE | Δ vs `base` | Seeds better | PCC | R² | Best epoch |
|---|---|---|---|---|---|---|---|
| `base` | 1.0142 ± 0.0028 | 1.0162 ± 0.0011 | 0 | | 0.9338 | 0.8718 | 180 / 194 / 176 |
| `aligned` (GE+Mut_CNV) | 0.9855 ± 0.0041 | 0.9941 ± 0.0047 | −0.0222 | 3 / 3 | 0.9368 | 0.8773 | 195 / 198 / 196 |
| `base+proteomics` | 1.0105 ± 0.0038 | 1.0134 ± 0.0063 | −0.0028 | 2 / 3 | 0.9341 | 0.8725 | 139 / 197 / 177 |
| `base+std_targets` | 0.9343 ± 0.0024 | 0.9423 ± 0.0080 | −0.0739 | 3 / 3 | 0.9437 | 0.8898 | 140 / 158 / 123 |
| **`final`** | 0.9610 ± 0.0022 | **0.9622 ± 0.0058** | −0.0540 | 3 / 3 | 0.9412 | 0.8851 | 48 / 67 / 46 |

`final` per seed (42 / 43 / 44): test 0.9602 / 0.9578 / 0.9687, Δ vs `base`
−0.055 / −0.059 / −0.048. MAE 0.7145, 3,644,525 parameters, 2,649 s per run
on laptop B. Gain over the per-drug mean floor of this split: +0.5642
(`base` +0.5102, `base+std_targets` +0.5842).

- **`final` beats `base` under the random split too**, by about the same
  amount as under the grouped split: −0.054 here (laptop B, 3 of 3 seeds)
  against −0.051 there (laptop A, 5 of 5). Both are past the ~0.03 band.
- **It is not the best configuration under this split.**
  `base+std_targets` alone is better than `final` on every seed: by 0.020 on
  test (+0.024 / +0.018 / +0.018) and 0.027 on validation. That is inside the
  band, but consistent in sign. So under the random split the three parts
  `final` adds to standardisation (proteomics, the molecular graph, the pair
  features) make it slightly worse together, where under the grouped split
  they make it better (1.2737 against 1.3086 for `base+std_targets`).
- **`final` peaks early here**, at epoch 46-67, while the other four peak at
  123-198. Which added part causes that was not tested: only single
  components and `final` were run under this split.
- `final` was selected on the grouped split's validation fold. Nothing was
  selected on this split, and this result does not change the choice: the
  grouped split is the headline protocol.

- **The two kept components swap roles between the splits.** Under the
  random split, where every cell line is seen in training, standardisation is
  the large effect (−0.074, all 3 seeds) and proteomics does nothing
  (−0.003). Under the grouped split proteomics is the consistent test gain
  (−0.033) and standardisation's is not established (−0.017).
- `aligned` here is 0.9941; the Phase 0 figure above is 0.978. They are not
  the same run: different omics (GE+Mut_CNV here, GE+Proteomics there),
  different script and different machine.
- Protocol gap for `base` (grouped − random): 1.3290 − 1.0162 = 0.313, both
  on laptop B. Compare with the 0.347 measured in
  [09](../docs/09_split_protocol_comparison.md).
- Protocol gap for `final` (grouped − random): 1.2737 − 0.9622 = 0.312. This
  one crosses machines (grouped on laptop A, random on laptop B); the anchor
  in section 2 puts the machine difference at 0.004 on the mean. The gap is
  the same size as `base`'s, so `final` does not close any of the distance
  between seen and unseen cell lines.
- All 15 prediction files are in the repository. The one that was missing
  (`base+std_targets__random__seed43.npz`) was added in `edb1bf8`.

## 10. Comparison with the other models on the same pairs

Same 111,799 pairs, same grouped split, test fold. The number of seeds
differs by row and is given in the notes; a single run carries the ~0.03
band and is not a ranking.

| Model | RMSE | PCC | Notes |
|---|---|---|---|
| Per-drug mean | 1.4889 | — | floor |
| RF, GE+Mut_CNV+Proteomics, fingerprint | 1.3805 | 0.8692 | single run, laptop A, `experiments/06_full_matrix`. Same omics as `final`. Best RF row: GE+Proteomics, 1.3712 |
| MLP, GE+Mut_CNV+Proteomics, fingerprint | 1.5795 | 0.8357 | single run, laptop A, `experiments/06_full_matrix`. Same omics as `final`; worse than the per-drug mean |
| MLP, Proteomics, fingerprint (E04) | 1.3273 ± 0.0286 | 0.8802 | 5 seeds, laptop A, [13](../docs/13_seed_variance_results.md). The best MLP row; its recorded single run was 1.2843 |
| CrossAttention, GE+Proteomics, fingerprint (E10) | 1.3029 ± 0.0126 | 0.8856 | 5 seeds, laptop A, [13](../docs/13_seed_variance_results.md) |
| HeteroIC50GNN, raw targets | 1.3203 ± 0.0007 | 0.8815 | **laptop B**, 3 seeds (42-44), validation 1.2628. The only trained PPI GNN in the report |
| HeteroIC50GNN, per-drug standardised targets | 1.3092 ± 0.0147 | 0.8840 | **laptop B**, 3 seeds, validation 1.2462 |
| `base`, laptop B anchor | 1.3290 ± 0.0169 | 0.8785 | **laptop B**, 5 seeds; the reference for the two rows above |
| `aligned` (MoGraphDRP-style) | 1.3416 ± 0.0088 | 0.8751 | 3 seeds, GE+Mut_CNV, this study |
| `base` (this ladder) | 1.3232 ± 0.0105 | 0.8797 | 3 seeds, this study |
| **Final model** (`final`) | **1.2737 ± 0.0144** | 0.8923 | laptop A, 5 seeds. `base` on the same machine and seeds: 1.3252 ± 0.0119 |

- **RF, MLP and CrossAttention are the earlier experiments' results, not
  re-run here** (filled in 2026-10-13 from `experiments/06_full_matrix/` and
  `experiments/13_seed_variance/`). Same pairs, split and floor; each was
  trained with its own script and settings, not the ladder's. The RF and the
  three-omics MLP are one run each.
- **`final` is ahead of every compared model with seeds:** by 0.029 against
  CrossAttention (1.3029), 0.054 against the best MLP (1.3273) and 0.036 to
  0.047 against `HeteroIC50GNN`. These are unpaired differences of means, and
  the CrossAttention one is at the edge of the ~0.03 band.

**`HeteroIC50GNN` (laptop B, run 2026-10-11, results committed in `2e425a9`,
`results/hetero_gnn/`).** The older heterogeneous GNN over the PPI,
drug-target and mutation graph, 2,815,489 parameters, run with seeds for the
first time. Per-seed test RMSE (42 / 43 / 44): raw 1.3208 / 1.3207 / 1.3196;
standardised 1.3005 / 1.3009 / 1.3261.

- **It is about level with `base` on the same machine** (1.3290 ± 0.0169):
  −0.009 with raw targets, −0.020 with standardised targets. Both are inside
  the ~0.03 band, the comparison is unpaired, and the two models are trained
  differently (this one: AdamW, lr 1e-3, early stopping with patience 15,
  best epoch 6-32; the ladder: Adam, lr 1e-4, 200 epochs).
- **It is behind the graph-bearing candidates** (1.2731 and 1.2737 on laptop
  A, section 5) by about 0.04 to 0.05. That comparison crosses machines; the
  anchor above puts the machine difference at 0.004 on the mean.
- Its validation RMSE with standardised targets (1.2462) equals the
  candidate's (1.2454) while its test RMSE is 0.03 worse.
- Earlier single runs on laptop A, seed 42: 1.3449 raw and 1.2917
  standardised. The same seed gives 1.3208 and 1.3005 on laptop B.
- Seeds 45 and 46 were not run. A run takes about a minute on laptop B.

## 11. Further checks (fill if run)

### Where `base`'s remaining error is (2026-10-05)

Exploratory, from the saved predictions of the 5 `base` runs in the gate CSV.
Ad-hoc analysis, no committed script yet; nothing here was trained.

- **26% of `base`'s squared error on test (23% on validation) is one constant
  per cell line:** a held-out cell line is more sensitive or more resistant
  across all drugs than predicted. The offset has a standard deviation of
  0.68 ln(IC50) units and is the same across seeds (correlation 0.95), so it
  is systematic. Removing it exactly would take test RMSE from 1.3252 to
  1.1367. A per-drug constant is 6% of the error.
- For scale: a pair-specific component that only acts on direct-hit pairs
  (1.26% of pairs) cannot move all-pairs RMSE by more than about 0.005.
- `base` already predicts 46% (test) to 61% (validation) of the variance of
  a cell line's general sensitivity level.

Candidate explanations of that per-cell-line offset, each fitted on one set
of 80 held-out cell lines and applied to the other 80:

| Candidate | Result | Reading |
|---|---|---|
| Doubling time (GDSC `growth_rate_20220907.csv`, one number per cell line, known for 524 of 531) | Spearman +0.41 / +0.50 with the offset (val / test). One slope on log doubling time: test 1.3252 → 1.2962 (−0.029, 5 of 5 seeds), validation 1.2961 → 1.2790 (−0.017, 5 of 5) | slow-growing cell lines are more resistant than `base` predicts; `base`'s omics predict doubling time only weakly (correlation 0.36 on held-out lines) |
| Tissue, cancer type, growth properties | corrections made RMSE worse in both directions (+0.003 to +0.018) | do not transfer between sets of cell lines |
| Omics PCs (GE, Mut_CNV, proteomics; ridge) | −0.004 to +0.039 | the offset is not linearly predictable from the model's own inputs |
| Expression of the drug's own target genes, against `base`'s per-pair error (293 drug-target pairs, 20 draws of random genes as the null) | mean absolute Spearman 0.094 / 0.103 against 0.093 / 0.097 for random genes; per-drug ridge +0.001 / +0.004 | no signal beyond random genes |

- **Doubling time is the only candidate with out-of-sample support.** It is
  a post-hoc linear correction here, at the edge of the ~0.03 band; it only
  counts once it is run as a ladder rung. It is an assay covariate, not an
  omics layer, and a known confound of IC50 (cite the growth-rate inhibition
  literature if used). Not run, not decided.

### Leave-drugs-out — unseen compounds (laptop B, 2026-10-13)

Whole drugs are held out instead of cell lines: 168 / 36 / 36 drugs
(78,973 / 16,878 / 15,948 pairs), the same split as
[10](../docs/10_leave_drugs_out_results.md). Every cell line is seen in
training. Run on laptop B, 10 runs, 3.7 h, results committed in `5efc470`
(`results/ldo_results.csv`); the `ldo` protocol is commit `e4e3517`.
`final` cannot run here, because a held-out drug has no statistics to
standardise with, so the model tested is `final` without standardisation:
`base+proteomics+mol_graph+pair_features`. Seeds 42-46.

The per-drug mean floor does not exist under this split. References on the
test fold: global training mean 2.5500; each cell line's own training mean
2.3874.

| Config | Validation RMSE | Test RMSE | PCC | R² | MAE | Gain over the mean-only floor (2.5500) | Best epoch (42-46) |
|---|---|---|---|---|---|---|---|
| `base` | 2.2021 ± 0.0283 | 2.1722 ± 0.0369 | 0.5909 | 0.2603 | 1.7220 | +0.3777 | 181 / 166 / 140 / 182 / 114 |
| `final` without standardisation | 2.0715 ± 0.0268 | **2.0038 ± 0.0843** | 0.6835 | 0.3698 | 1.5440 | +0.5462 | 20 / 30 / 165 / 185 / 47 |

| Δ, `final` without standardisation − `base` | Δ | Seeds better | 95% CI, seeds | 95% CI, drugs | Per-seed Δ (42 / 43 / 44 / 45 / 46) |
|---|---|---|---|---|---|
| validation | −0.1306 | 5 / 5 | −0.187 to −0.075 | −0.315 to +0.023 | −0.161 / −0.105 / −0.107 / −0.194 / −0.086 |
| test | −0.1685 | 4 / 5 | −0.314 to −0.023 | −0.441 to +0.091 | −0.292 / −0.239 / −0.114 / +0.007 / −0.205 |

"Drugs" is a paired bootstrap over the fold's 36 held-out drugs, the unit
that is resampled when a new compound arrives. It and the per-drug figures
below are an ad-hoc analysis of the saved predictions, with no committed
script.

- **Both models beat the mean-only floor on unseen drugs**, so the drug
  encoders transfer something to compounds they have not seen. R² is 0.26 and
  0.37, against 0.77 and 0.79 on unseen cell lines: unseen drugs are the much
  harder setting.
- **The final model's parts help on average, but not established across
  drugs.** −0.169 on test is far larger than the ~0.03 band and the seed
  interval excludes 0, but the interval over drugs includes 0 on both folds.
  It is better on 20 of the 36 test drugs and 20 of the 36 validation drugs.
  The mean gain comes from a few drugs (drug `2011`: 3.06 to 1.05; `1020`:
  4.09 to 2.22) while others get worse (`1088`: 2.52 to 3.58). With 36 drugs
  per fold, which drugs are held out matters more than the seed.
- **Most of the error is the drug's overall level.** 76% of `base`'s squared
  test error and 70% of the other model's is one constant per drug: the
  model misjudges how potent an unseen compound is on average. With that
  removed, the two are level (1.064 and 1.093). So the gain is in guessing a
  new drug's mean response, not in ranking cell lines for it.
- **It cannot be attributed to one component.** Three things differ between
  the two configs (proteomics, the molecular graph, the pair features). The
  expectation is that the molecular graph matters most for unseen drugs
  (section 6b), but no single-component rung was run under this split.
- **Not better than the earlier flat models on this split.** The earlier
  single runs in [10](../docs/10_leave_drugs_out_results.md) were MLP 1.8516,
  RF 1.8719 and CrossAttention 1.9200 (fingerprint, GE+Proteomics, laptop A).
  This model is 2.0038 ± 0.0843, its best seed 1.9286. The comparison is a
  single run against a 5-seed mean, with different omics, training settings
  and machine, so it is not a ranking; but there is no evidence that the
  final model generalises to new compounds better than a flat MLP.
- **Validation does not select well here.** The validation drugs are 36
  different compounds. Across seeds, validation and test RMSE of the larger
  model are anti-correlated (−0.91): seeds 44 and 45 have the lowest
  validation RMSE, checkpoint late (epoch 165 and 185) and are the worst on
  test. Its test std (0.084) is more than twice `base`'s.
- By target annotation (test: 23 drugs with a known target, 13 without) the
  gain is similar in both groups (−0.186 and −0.135) and both drug intervals
  include 0. On the 23 test drugs with a target, both models are worse than
  the mean-only floor of that subset (2.21 and 2.03 against 2.00).
### Interpretability — what `final` does with the direct-hit signal (2026-10-13)

The test planned for Phase 6 read the pair attention's weights; that
attention is near uniform and is not in the final model. This is its minimal
replacement, on the one pair-specific input `final` has: the five pair
features. It is a sensitivity test on hand-built inputs. It is not a
gene-level attribution and not the SHAP module of the proposal. Script:
`final_model/pair_feature_check.py`; output
`results/final_weights_results_pair_feature_check*.csv`.

**Weights.** No `final` run had kept its weights, so `final` was trained
again on seeds 42-46 with `--save-weights` (laptop A, 2026-10-13 11:47 to
13:34, 6,425 s, `results/final_weights_results.csv`). These are new draws,
not the ladder's rows (section 2): validation 1.2442 ± 0.0092 and test
1.2698 ± 0.0126, against 1.2441 ± 0.0102 and 1.2737 ± 0.0144 in the ladder;
single seeds differ by up to 0.018. Against `base`: −0.055, 5 of 5 seeds. A
second set of five runs of `final` that agrees with the first.

**Method.** A direct-hit pair is predicted twice: with its real features,
and with the hit removed, meaning the same pair described as if its nearest
mutation were one hop from a target (`direct` 1 to 0, `proximity` 1.0 to
0.5; an input that 54,844 real pairs have). The *effect* of the hit is the
first prediction minus the second, in ln(IC50) units; negative means the hit
makes the model predict the pair more sensitive. The omics inputs are not
changed. The null is 50 untrained models of the same architecture, as in
[15](../docs/15_interpretability_validation_results.md). The three pass
rules were written in `CLAUDE.md` before the weights existed.

**1. Use — passed.** Direct-hit pairs of held-out cell lines:

| Fold | Pairs | RMSE, real features | RMSE, hit removed | Seeds where real is better | Signed error, real | Signed error, removed | Mean effect |
|---|---|---|---|---|---|---|---|
| validation | 213 | 1.3001 ± 0.040 | 1.6524 ± 0.040 | 5 / 5 | +0.161 | +0.539 | −0.378 |
| test | 211 | 1.4100 ± 0.053 | 1.6427 ± 0.062 | 5 / 5 | −0.016 | +0.461 | −0.477 |

Rule: lower RMSE with the real features on all 5 seeds on both folds. Met.
With the hit removed the trained model shows the same bias as models that
never had the features (about +0.4 to +0.5, sections 5 and 6f).

**2. Direction across drugs — passed.** For the 30 drugs with at least 3
direct-hit pairs on held-out cell lines, the model's mean effect against the
measured difference between hit and non-hit cell lines, both in units of the
drug's training std:

| Rows | Drugs | Spearman, trained (5 seeds) | Seeds above 0 | Untrained (50) | Controls reaching the trained mean |
|---|---|---|---|---|---|
| validation + test | 30 | **+0.725 ± 0.049** | 5 / 5 | +0.034 ± 0.312 | 0 / 50 |
| test only | 20 | +0.730 ± 0.047 | 5 / 5 | +0.008 ± 0.435 | 1 / 50 |

Rule: positive on all 5 seeds, at most 2 of 50 controls reaching the mean.
Met. Per drug (held-out cell lines, model column averaged over seeds;
negative = more sensitive when the target is mutated):

| Drug | Target | Direct-hit pairs | Measured, hit − non-hit (z) | Model's effect (z) |
|---|---|---|---|---|
| Gefitinib | EGFR | 6 | −2.81 | −0.51 |
| Osimertinib | EGFR | 6 | −2.80 | −0.52 |
| Dabrafenib | BRAF | 14 | −2.02 | −1.18 |
| PLX-4720 | BRAF | 15 | −1.79 | −1.03 |
| SB590885 | BRAF | 15 | −1.33 | −1.07 |
| Afatinib | EGFR, ERBB2 | 11 | −1.34 | −0.41 |
| Alpelisib | PI3Kalpha | 20 | −0.81 | −0.40 |
| Buparlisib | PI3K (pan) | 20 | −0.31 | −0.15 |
| AZD8186 | PI3Kalpha, PI3Kbeta | 20 | −0.02 | −0.19 |
| SGC-CBP30 | EP300, CBP | 13 | +0.25 | −0.18 |
| MIRA-1 | TP53 | 115 | +0.39 | −0.02 |

All 30 are in `results/final_weights_results_pair_feature_check_drugs.csv`.

- **The model's effect is ordered like the measurements**: largest for the
  BRAF and EGFR inhibitors, near zero for MIRA-1, whose TP53 hits are a
  quarter of all direct hits and are not sensitising.
- **It is smaller than the measured difference**: mean absolute size 0.36
  against 0.88. The two are not the same quantity. The measured column is
  the whole difference between hit and non-hit cell lines; the model's
  effect is what the pair features add on top of what it already reads from
  the omics. So a smaller effect is expected, and is not by itself an error.
- **It is almost always sensitising**: negative for 29 of the 30 drugs. The
  sign agrees with the measurement for 24. Of the 7 drugs whose hit cell
  lines are measured as more *resistant*, the model has the sign right for
  one. The features did not teach it resistance.

**3. BRAF case study — passed.** BRAF-mutant cell lines on PLX-4720,
SB590885 and Dabrafenib, the 157 pairs of experiment 15. Spearman
correlation across pairs between the model's effect and measured ln(IC50):

| Pairs | n | Trained (5 seeds) | Seeds above 0 | Within drug | Untrained (50) | Controls reaching |
|---|---|---|---|---|---|---|
| all | 157 | +0.695 ± 0.025 | 5 / 5 | +0.697 | −0.027 ± 0.193 | 0 / 50 |
| training cell lines | 113 | +0.716 ± 0.028 | 5 / 5 | +0.712 | −0.019 ± 0.212 | 0 / 50 |
| **held-out cell lines** | **44** | **+0.633 ± 0.026** | 5 / 5 | +0.658 | −0.049 ± 0.263 | 0 / 50 |

Rule: on the 44 held-out pairs, positive on at least 4 of 5 seeds and at
most 2 of 50 controls reaching the mean. Met. "Within drug" is the same
correlation in each drug's own z-units, so a difference in potency between
the three drugs cannot produce it.

- The hit is most sensitising where the cell line is measured as most
  sensitive, on cell lines the model never trained on. Experiment 15 found
  the same on the old graph model (0.65 overall, 0.54 held-out).
- Mean effect on held-out pairs: −1.69 ln(IC50) units, sensitising in 99% of
  pairs. Untrained models: +0.001 ± 0.005 and 43% ± 29%. As in experiment
  15, the share alone would not separate them from a trained model that is
  merely consistent; the correlation does.
- On these 44 pairs RMSE is 1.30 with the hit and 2.62 with it removed.
- **Caveat, the same as in experiment 15.** The effect correlates +0.47 with
  what the model predicts without the hit, so part of the +0.63 is
  inherited from the omics-based prediction: the model scales the hit by
  context it already reads elsewhere. And the hit mostly shifts the level
  rather than re-ranking cell lines: the rank correlation of the prediction
  with measured ln(IC50) is +0.80 with the hit and +0.78 without.

**What this shows and does not show.** The trained final model uses the
direct-hit signal, in the direction and in the drug order the measurements
support, and untrained models do not. It does not show a mechanism, it does
not explain any prediction outside direct hits (1.26% of pairs), and it says
nothing about the omics inputs, which are PCA components. The four other
pair features were not tested separately.

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
- [x] Component f was tested as hand-built features and as attention without message passing; PPI message passing (stage B) was skipped (decided 2026-10-06), so no trained rung in the ladder is a GNN and the PPI graph is tested only indirectly (Phase 3 hop features; `HeteroIC50GNN` as a compared baseline). The one subset with any signal (a target is directly mutated) is 1.26% of pairs, 211 on the test fold. The features' gain there is better on 4 of 5 seeds, but its test intervals include 0.
- [x] The attention form of f was rejected as configured, with no tuning (embedding width 64, no dropout on the attention, no weight decay). Its attention did not learn to find the mutated target (0.02 above chance), so the result does not test the idea itself. The graph has no protein features; 15 of the 251 proteins mutated in test cell lines appear in no training cell line.
- [x] The graph component of the final model is chosen under a tolerance rule (validation no more than 0.03 worse than the candidate), not because it improves accuracy: neither `mol_graph` nor the pair features changes all-pairs test RMSE beyond noise (2026-10-12). `pair_attention` meets the tolerance by 0.0009, inside its seed std.
- [x] `final` was selected under the grouped split only. Under the random split it beats `base` (−0.054, 3 of 3 seeds) but is 0.020 worse than `base+std_targets` alone on every seed (section 9). The random-split rows are from laptop B with 3 seeds.
- [x] Unseen drugs: the leave-drugs-out result rests on 36 test drugs. The gain over `base` (−0.169) has a seed interval that excludes 0 and a drug interval that does not; three components differ between the two configs, so it is not attributed; validation drugs do not select a good checkpoint for test drugs; and the model does not beat the earlier flat MLP's single run (1.8516 against 2.0038 ± 0.0843). Standardisation, one of the two components the grouped-split gain comes from, cannot be used there.
- [x] Interpretability is a sensitivity test on five hand-built inputs, on 1.26% of pairs. No gene-level attribution and no SHAP: the omics inputs are PCA components. The model's direct-hit effect is sensitising for 29 of 30 drugs, so it does not capture the drugs where a mutated target goes with resistance. The weights it was run on are a second set of five runs of `final`, not the ladder's rows.
- [x] The subset bootstrap resamples cell lines only; training noise is covered separately by a t-interval over 5 seeds.
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
> the base. Hand-built features relating each drug's targets to each cell
> line's mutations did not change overall RMSE (−0.004, 5 seeds). The base
> model over-predicted ln(IC50) by about 0.4 on the 1.3% of pairs where a
> target is itself mutated, and the features removed most of that bias on
> every seed. RMSE on those pairs fell from 1.72 to 1.58 (better on 4 of 5
> seeds), a difference whose 95% interval still includes zero. Replacing the
> features with learned attention between each drug's target proteins and
> each cell line's mutated proteins made the model worse overall
> (1.3861 ± 0.0142, +0.061, all 5 seeds) and no better than the features on
> the directly-mutated pairs. The trained attention gave a mutated target's
> own protein only 0.02 more weight than chance, so this result rejects the
> module as built rather than the idea. Message passing over the PPI network
> beneath that attention was not run.

Addition from the 5-seed runs of 2026-10-07 (replaces "none improved on the
base clearly" above once the final model is confirmed):

> With five seeds, adding proteomics improved on the base on every seed
> (test RMSE 1.2920 ± 0.0125, −0.033), and combining it with per-drug target
> standardisation gave 1.2797 ± 0.0192 (−0.046, every seed), against 1.4889
> for a per-drug mean. The two components were chosen on validation RMSE,
> where the combination was better than either alone; on test it was not
> clearly better than proteomics alone (−0.012, 3 of 5 seeds). Differences
> are over training seeds on one fixed split of cell lines.

Addition from the graph decision of 2026-10-12 (pick the sentence for the
configuration the team chooses):

> Adding a graph convolution over each drug's molecular graph to that model
> left accuracy unchanged within seed-to-seed variation (test RMSE
> 1.2731 ± 0.0182 against 1.2797 ± 0.0192; validation 0.011 worse). Adding
> five features computed from the protein-interaction, drug-target and
> mutation graph also left overall RMSE unchanged (1.2737 ± 0.0144), but on
> the 1.3% of pairs where a target of the drug is itself mutated it removed
> most of the model's bias and lowered RMSE on every seed (test 1.58 to 1.38;
> the interval over test cell lines still includes zero). Learned attention
> between target and mutated proteins made the model worse (1.3225 ± 0.0285,
> every seed).

Addition from the random split of 2026-10-12 (sanity check, second machine):

> Under a random pair split, where every test cell line is also seen in
> training, the final model reached 0.9622 ± 0.0058 against 1.0162 ± 0.0011
> for the base (3 seeds), a gain of the same size as on unseen cell lines.
> On that split per-drug target standardisation alone did slightly better
> (0.9423 ± 0.0080), and proteomics alone had no effect. Errors are about 0.31
> lower than on unseen cell lines for both the base and the final model, so
> the final model does not narrow that gap.

Addition from leave-drugs-out of 2026-10-13 (second machine):

> On compounds never seen in training (36 held-out drugs), the final model
> without target standardisation, which needs a drug's training responses,
> reached 2.0038 ± 0.0843 against 2.1722 ± 0.0369 for the base and 2.5500 for
> predicting the global mean (5 seeds). The improvement held on 4 of 5 seeds
> but not consistently across drugs (20 of 36; the interval over drugs
> includes zero), and about three quarters of the remaining error is a
> misjudged overall potency of the new compound. An earlier flat MLP reached
> 1.85 in a single run on the same split. The model is therefore presented as
> a predictor for new cell lines, not for new compounds.

Addition from the interpretability test of 2026-10-13:

> To test whether the model uses the biological prior it is given, each pair
> in which a drug's target is itself mutated was predicted twice, with its
> real pair features and with that signal removed. On held-out cell lines
> removing it raised RMSE on those pairs from 1.41 to 1.64 and brought back
> a bias of +0.46 ln(IC50) (all five seeds). Across the 30 drugs with enough
> such pairs, the size of the model's effect followed the measured difference
> between mutant and non-mutant cell lines (Spearman 0.73), and for BRAF
> inhibitors it was strongest in the BRAF-mutant cell lines measured as most
> sensitive (0.63 on held-out cell lines); none of 50 untrained models
> reached either value. The effect was sensitising for almost every drug, so
> the model does not capture targets whose mutation goes with resistance.
> This is a test of five hand-built inputs on 1.3% of pairs, not a
> gene-level explanation of the model.

## 15. Reproduction

```
python "src/final_model/run_ablation.py" --list
python "src/final_model/run_ablation.py" --protocols grouped --seeds 42 43 44
python "src/final_model/run_ablation.py" --configs aligned base full --protocols random
```

The second command produced the 24 rows in this file (run 2026-10-03, 15:33 to
22:05, as one process on a GTX 1650; launch details in `CLAUDE.md`). The third
was never run in that form: the random split was run on laptop B with a
different set of configs (below), and `full` has no random-split row.

Commit hash of the code that produced the CSV: `66205f8`.

Five-seed rows for the candidates (run 2026-10-06 23:44 to 2026-10-07 01:29,
same machine and CSV; first command from commit `1c00b00`, second from
`70c3d13`, which differs only in `CLAUDE.md`; rows already present were
skipped, so 4 and then 7 runs):

```
python src/final_model/run_ablation.py --configs base base+proteomics \
  --protocols grouped --seeds 42 43 44 45 46
python src/final_model/run_ablation.py \
  --configs base base+proteomics base+std_targets base+proteomics+std_targets \
  --protocols grouped --seeds 42 43 44 45 46
```

The ladder CSV now has 35 rows. These 11 runs also saved predictions, in
`results/ablation_results_predictions/`; the 24 Phase 2 runs did not.

Graph decision (run 2026-10-11 21:51 to 2026-10-12 03:40, same machine and
CSV, commit `2e425a9`; the 10 `base` and candidate rows were skipped, so 15
runs in 20,964 s):

```
python src/final_model/run_ablation.py \
  --configs base base+proteomics+std_targets \
            base+proteomics+std_targets+mol_graph \
            base+proteomics+std_targets+mol_graph+pair_features \
            base+proteomics+std_targets+mol_graph+pair_attention \
  --protocols grouped --seeds 42 43 44 45 46
C=base+proteomics+std_targets
python src/final_model/pair_gate_report.py \
  --results src/final_model/results/ablation_results.csv --tag graph_decision \
  --baseline $C --configs $C $C+mol_graph $C+mol_graph+pair_features $C+mol_graph+pair_attention
python src/final_model/pair_gate_report.py \
  --results src/final_model/results/ablation_results.csv --tag graph_decision_vs_mol_graph \
  --baseline $C+mol_graph --configs $C+mol_graph $C+mol_graph+pair_features $C+mol_graph+pair_attention
```

The ladder CSV now has 50 rows. The 15 runs saved predictions, and the 5
attention runs saved weights in `results/ablation_results_checkpoints/`
(18 MB each, 91 MB in all).

Laptop B blocks (a teammate's machine; each writes its own file, so no row
of laptop B is in the ladder CSV). Anchor and `HeteroIC50GNN` run 2026-10-11
(results committed in `dcbc304` and `2e425a9`); random split run 2026-10-11
for the four single configs (`e34f390`) and 2026-10-12 for `final`
(`c346fc9`, after pulling the commit that names it, `e1fd48f`):

```
python -u src/final_model/run_ablation.py --configs base \
  --protocols grouped --seeds 42 43 44 45 46 \
  --out src/final_model/results/laptop_b_anchor_results.csv
python -u experiments/07_gnn_ablation/hetero_ic50_gnn_matrix.py --seed 42 \
  --out src/final_model/results/hetero_gnn/raw_seed42.csv      # and seeds 43, 44
python -u experiments/07_gnn_ablation/hetero_ic50_gnn_matrix.py --seed 42 \
  --standardize-targets \
  --out src/final_model/results/hetero_gnn/std_seed42.csv      # and seeds 43, 44
python -u src/final_model/run_ablation.py \
  --configs base aligned base+proteomics base+std_targets final \
  --protocols random --seeds 42 43 44 \
  --out src/final_model/results/random_split_results.csv
```

Leave-drugs-out (laptop B, results committed 2026-10-13 in `5efc470`, 10
runs, 13,270 s; the protocol is commit `e4e3517`):

```
python -u src/final_model/run_ablation.py \
  --configs base base+proteomics+mol_graph+pair_features \
  --protocols ldo --seeds 42 43 44 45 46 \
  --out src/final_model/results/ldo_results.csv
```

Interpretability test (laptop A, 2026-10-13; the retrain 11:47 to 13:34, 5
runs in 6,425 s, from `5efc470` plus the `--save-weights` flag, which was not
yet committed; the check itself trains nothing):

```
python src/final_model/run_ablation.py --configs final --protocols grouped \
  --seeds 42 43 44 45 46 --save-weights \
  --out src/final_model/results/final_weights_results.csv
python src/final_model/pair_feature_check.py
```

The weights are in `results/final_weights_results_checkpoints/` (15 MB each,
70 MB in all). Runs with `mol_graph` are not bit-reproducible (section 2), so
a repeat of the first command gives runs within seed noise of these, not the
same weights.

The interval over drugs and the per-drug figures in section 11 were computed
ad hoc from `results/ldo_results_predictions/` with `bootstrap_delta` from
`pair_gate_report.py`, grouping by drug; `pair_gate_report.py` itself only
reads grouped-split rows.

Phase 3 gate (run 2026-10-04, same machine, commit `cc303ee`; seeds 42-44
19:42 to 20:40, seeds 45-46 20:57 to 21:36, 10 runs):

```
python src/data/pair_features.py
python src/final_model/run_ablation.py --configs base base+pair_features \
  --protocols grouped --seeds 42 43 44 45 46 --out src/final_model/results/pair_gate_results.csv
python src/final_model/pair_gate_report.py --configs base base+pair_features
```

`--configs` was added on 2026-10-05. The report's bootstrap draws depend on
which configs are reported together, so naming the two Phase 3 configs keeps
these files reproducible after other rungs are added to the same CSV
(checked: byte-identical output).

Phase 4 stage A (run 2026-10-05, same machine, commit `9c7ba4f`, 15:48 to
17:14, 5 runs; the 10 Phase 3 rows were skipped and reused):

```
python -m src.models.pair_graph_drp
python src/final_model/run_ablation.py --configs base base+pair_features base+pair_attention \
  --protocols grouped --seeds 42 43 44 45 46 --out src/final_model/results/pair_gate_results.csv
python src/final_model/pair_gate_report.py --tag attention_vs_features \
  --configs base+pair_features base+pair_attention --baseline base+pair_features
python src/final_model/pair_gate_report.py --tag attention_vs_base \
  --configs base base+pair_attention
python src/final_model/pair_attention_check.py
```

The last command reads the trained weights in
`results/pair_gate_results_checkpoints/` (17 MB per run). It was written
after the run was launched and is not part of commit `9c7ba4f`; it only
reads the checkpoints.
