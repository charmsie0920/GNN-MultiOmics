# Ablation Study — Base + One Component at a Time

**Script:** [`src/final_model/run_ablation.py`](./final_model/run_ablation.py)
**Model:** [`src/models/mographdrp_aligned.py`](./models/mographdrp_aligned.py) (frozen Phase 0 baseline, imported by the script)
**Raw results:** `src/final_model/results/ablation_results.csv`
**Companions:** [`results.md`](../docs/results.md), [`13_seed_variance_results.md`](../docs/13_seed_variance_results.md), [`09_split_protocol_comparison.md`](../docs/09_split_protocol_comparison.md)
**Date:** 2026-10-03 · **Runtime:** 6.5 h (24 runs, 23,523 s) · **Hardware:** one laptop GTX 1650 (4 GB), torch 2.13.0+cu130
**Status:** Phase 2 filled (components a-e, 3 seeds, grouped split). Phase 3 filled 2026-10-04 (feature gate for component f, 10 runs over 5 seeds, 98 min, commit `cc303ee`): gate not passed, see section 6f. The attention form of f is not built; sections that depend on it, on the random split or on the final-model choice are still `TBD`.

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
| Seeds | 42, 43, 44 (3 per row); the Phase 3 gate rows use 42-46 (5 per arm) |
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
| f | Pair-specific component. Tested form: hand-built pair features. Built, not run: attention of target proteins over mutated proteins (Phase 4 stage A). Planned: PPI message passing under that attention (stage B) | `pair_module=features` or `pair_module=attention` (alternatives, never combined) | adds 5 per-pair features, or a 64-wide attention vector, to the predictor input | this project (novelty component) | the base sees the cell line and the drug separately; nothing tells it whether this cell line's mutated proteins are this drug's targets |

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
| `aligned` | b + c | 1.3416 ± 0.0088 | +0.0184 | no | 0.8751 | 0.7642 | +0.1473 | 3,528,562 | 1,254 |
| `full` | a-e | 1.3394 ± 0.0081 | +0.0162 | no | 0.8789 | 0.7650 | +0.1495 | 4,270,706 | 2,309 |

Negative Δ means better than base. Fit is the mean training time per run.
`base+pair_features` comes from the Phase 3 gate run
(`pair_gate_results.csv`) with 5 seeds (42-46). Its Δ is against the gate's
own 5-seed `base`, 1.3252 ± 0.0119; on seeds 42-44 that re-run matched the
`base` row above exactly. Subset results are in section 6f.

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

### f. Pair-specific component — feature gate run (Phase 3); attention built, not run (Phase 4 stage A)

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

**Stage A, `base+pair_attention`: built 2026-10-05, not run.** No result yet;
this is what the rung is, so the numbers can be read against it when it runs.

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
  learn that, from 989 direct-hit training pairs. Whether it did is checked
  from the saved weights after the run.
- **Direct hits are concentrated:** 85 of the 240 drugs have any; drug `1931`
  has 372 of the 1,413; the median such drug has 6 in the training fold.
- **Cost:** +1,087,232 parameters (4,271,489 total): 1,037,824 in the
  embedding table, 16,640 in the attention layer, 32,768 in the widened
  predictor. At most 636 embedding rows (40,704 values) can receive a
  gradient; the rest of the table is unused until stage B. Timed over 3
  epochs on the GTX 1650: 5.4 s per epoch against 2.8 s for `base` (1.9×),
  about 18 min per 200-epoch run; peak GPU memory 227 MiB.
- **Pass rule, fixed before the run:** on direct-hit pairs, better than
  `base+pair_features` on at least 4 of 5 seeds on validation, with the
  validation cell-line interval excluding 0; and no worse than `base` on all
  pairs beyond the ~0.03 band. Test is reported after, not used to decide.
- **Checks passed before the run:** 1-epoch CPU `base` (val 2.688133 / test
  2.747617) and `base+pair_features` (2.554799 / 2.649105) are unchanged by
  the new code; the Phase 3 feature matrix is byte-identical; the padded
  target and mutation sets reproduce the Phase 3 direct-hit flag on all
  111,799 rows; the report regenerates the two Phase 3 CSVs exactly.

- **GNN depth ablation** (does message passing matter?): not run yet; the
  0-layer row is Stage A of Phase 4.

| Message-passing layers | RMSE (all pairs) | RMSE (pairs with a known target) |
|---|---|---|
| 0 (attention on raw embeddings, no GNN) | not run | not run |
| 1 | not run | not run |
| 2 | not run | not run |

Source: `src/final_model/results/pair_gate_results.csv`,
`pair_gate_results_subsets.csv`, `pair_gate_results_paired.csv`; produced by
`src/final_model/pair_gate_report.py --configs base base+pair_features`.

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

- **Chosen configuration:** not chosen yet. This is a team decision; the
  evidence from Phase 2 is:
  - best on test: `base+proteomics`, 1.2930 ± 0.0169 (PCC 0.8871, R² 0.7810);
  - best on validation: `base+std_targets`, 1.2580 ± 0.0097 (test 1.3108);
  - `full` (base + a-e) is not better than `base` and should not be presented
    as the final model on these numbers.
- **Components dropped and why:** a (+0.0912, hurts), b (+0.0087, no effect),
  c (+0.5111, hurts and unstable as configured), e (−0.0124, no effect), f as
  hand-built features (−0.0035, no effect on all pairs). d (−0.0302) is
  borderline.
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
- [x] Component f was tested only as hand-built features. Its one subset with any signal (a target is directly mutated) is 1.26% of pairs, 211 on the test fold. The gain there is better on 4 of 5 seeds, but its test intervals include 0.
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
> seeds), a difference whose 95% interval still includes zero.

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
