# Split Protocol Comparison — Why Our Numbers Aren't Comparable to MoGraphDRP's

**Script:** [`experiments/09_split_protocol/split_comparison.py`](../experiments/09_split_protocol/split_comparison.py)
**Raw results:** `experiments/09_split_protocol/split_comparison_results.csv`
**Companions:** [`results.md`](./results.md), [`08_ensemble_refinement_results.md`](./08_ensemble_refinement_results.md), [`06_cross_attention_ablation_results.md`](./06_cross_attention_ablation_results.md)
**Date:** 2026-08-24 · **Runtime:** 428.2s

## The question

MoGraphDRP reports **RMSE 0.6622**. Our best model reports **RMSE 1.2442**. Those
look like they say something about model quality — but they are measured under
different train/test protocols, so they measure **different tasks**.

This experiment runs our best model (cross-attention, GE+Proteomics, one-hot —
row E01 in [`results.md`](./results.md)) under **both** protocols with
everything else held fixed, to quantify how much of the apparent gap is
protocol rather than architecture.

| | MoGraphDRP (their §2.5) | This project |
|---|---|---|
| Split unit | Individual (cell line, drug) **pairs** | Whole **cell lines** |
| Ratio | Random 80/10/10 | 70/15/15 (`GroupShuffleSplit`) |
| Can a test cell line appear in training? | Yes | No (asserted in code) |
| Task measured | **Imputation** — fill a gap in a partly-screened row | **Generalization** — predict for an unseen cell line |

Neither is wrong. A random split is the correct evaluation for imputation, and
imputation is explicitly what the paper does (their §3.4: the model "was able to
impute 9,117 missing IC50 values"). Our proposal (§1, §4) instead describes a
researcher uploading a **new** cell line, which requires the grouped split.

## Measured leakage

The `leakage_report()` helper quantifies exactly what each protocol exposes:

| | Grouped by cell line | Random pairs |
|---|---|---|
| Test cell lines | 80 | 530 |
| ...of which also in training | **0** | **530 (all)** |
| Test rows whose cell line was seen in training | **0.0%** | **100.0%** |
| Mean training rows per test cell line | **0** | **208** |

Under the random split, *every* test prediction concerns a cell line the model
has already seen an average of **208 times**. Each of the 531 cell lines
contributes ~279 measurements sharing one omics profile, so random shuffling
scatters them across all three folds. The model does not need to generalize to a
new biological sample; it needs to interpolate within a profile it has largely
memorized.

## Results

| Protocol | RMSE | MAE | R² | PCC | SCC | AUC | F1 | Floor |
|---|---|---|---|---|---|---|---|---|
| Grouped by cell line | 1.2442 | 0.9336 | 0.7883 | 0.8897 | 0.8480 | 0.9191 | 0.8423 | 2.7097 |
| Random pairs | **0.8971** | 0.6668 | 0.8955 | 0.9464 | 0.9194 | 0.9532 | 0.8830 | 2.7749 |
| *MoGraphDRP (published)* | *0.6622* | *0.4884* | *0.9388* | *0.9689* | *0.9524* | *—* | *—* | *—* |

**Protocol effect: 1.2442 → 0.8971, a 0.3471 RMSE improvement from changing
nothing but the split.** That is 28% of our grouped-split error, produced
entirely by allowing cell-line leakage.

## Decomposing the apparent gap

```
Apparent gap  (our grouped 1.2442  vs their 0.6622)  =  0.5820
  of which protocol   (1.2442 -> 0.8971)             =  0.3471   (60%)
  of which everything else (0.8971 -> 0.6622)        =  0.2349   (40%)
```

**They remain genuinely ahead by 0.2349 RMSE on a like-for-like protocol** —
that part is real and should not be explained away. Plausible contributors,
all deliberately out of scope for our matrix:

- **Four omics types** (they add DNA methylation and GSVA pathway-activity
  scores) versus our three.
- **Three fingerprint types** (Morgan + PubChem + ESPF) fused by attention,
  versus our single Morgan fingerprint.
- **An atom-level molecular graph GCN** per drug — a representation we
  explicitly scoped out in favour of the heterogeneous PPI graph.
- **Bilinear attention** for drug–cell interaction, versus our concatenate-then-MLP head.
- COSMIC gene filtering, reducing each omics branch to ~600–700 cancer-relevant
  genes rather than PCA over the full feature set.

## The XGBoost leakage hypothesis: confirmed

[`08_ensemble_refinement_results.md`](./08_ensemble_refinement_results.md) recorded a
falsifiable prediction — that the paper's +19.7% refinement gain depends on
cell-line leakage, and that refinement should therefore start helping under a
random split. It does:

| Protocol | Base RMSE | + XGBoost | Δ |
|---|---|---|---|
| Grouped by cell line | 1.2442 | 1.2557 | **−0.92%** |
| Random pairs | 0.8971 | 0.8885 | **+0.96%** |

**The sign flips.** Identical refiner, identical hyperparameters, identical base
architecture — only the split differs. This confirms the mechanism: a residual
corrector can learn "this cell line reads 0.3 high" and apply it at test time
*only if it has seen that cell line*. With 100% of test rows carrying ~208 prior
examples of their own cell line, it can; with 0%, there is nothing to transfer
and the refiner only adds variance.

**Caveat on magnitude:** our +0.96% is far short of their +19.7%, so leakage
explains the *direction* but not the full size. Their interaction vector comes
from a multi-head bilinear attention module and likely carries more correctable
structure than our head's penultimate activation. Leakage is a necessary
condition for the gain, not a complete account of it.

## What to report

1. **Keep the grouped split as the primary protocol.** It matches the stated use
   case in the proposal and is the defensible number for the report.
2. **Never present 1.2442 and 0.6622 side by side without qualification.** They
   measure different tasks. Use the decomposition above if a direct comparison
   is needed.
3. Suggested wording:

   > Results are reported under a cell-line-grouped split, measuring
   > generalization to previously unseen cell lines. MoGraphDRP's published
   > RMSE of 0.6622 uses a random split in which every test cell line also
   > appears in training (~208 measurements each), evaluating imputation rather
   > than generalization. Run under that same protocol, our model achieves
   > RMSE 0.8971; the remaining 0.235 gap reflects genuine architectural
   > differences.

4. The random-split row is a **comparability bridge only** — it should not be
   cited as a headline result, and `random_pair_split()` carries a docstring
   saying so.

## Caveats

- One model, one configuration. The protocol effect is measured on E01 only;
  its magnitude may differ for other architectures (a graph model with learned
  protein embeddings could plausibly exploit leakage differently).
- The random split uses the paper's 80/10/10 ratio rather than our 70/15/15,
  so the training set is slightly larger (107,811 vs 93,507 pairs) — a small
  confound in the protocol comparison, retained deliberately to match their
  setup exactly.
- The random-split run needed 96 epochs to converge versus 24 for the grouped
  split, consistent with there being more exploitable structure to fit.
- A third protocol, **leave-drugs-out** (grouped by drug), remains unrun. It is
  the only setting that would test whether Morgan fingerprints earn their place,
  since structural encoding can generalize to unseen drugs where one-hot
  structurally cannot. See
  [`phase2_wide_matrices_and_split_plan.md`](./phase2_wide_matrices_and_split_plan.md).
