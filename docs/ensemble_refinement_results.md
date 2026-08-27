# XGBoost Ensemble Refinement — Negative Result

**Scripts:** [`experiments/Ensemble Refinement/ensemble_baseline.py`](../experiments/Ensemble%20Refinement/ensemble_baseline.py), [`src/models/ensemble_refinement.py`](../src/models/ensemble_refinement.py)
**Raw results:** `experiments/Ensemble Refinement/ensemble_results.csv`
**Companions:** [`cross_attention_ablation_results.md`](./cross_attention_ablation_results.md), [`gnn_ablation_results.md`](./gnn_ablation_results.md), [`results.md`](./results.md)
**Date:** 2026-08-24 · **Runtime:** 2 base models in 225.7s

## What this tests

MoGraphDRP's final architectural component (their §2.4). Rather than using the
neural model's output directly, its compressed interaction vector `f` and its
initial prediction `ŷ` are concatenated (`Z = f ⊕ ŷ`, their eq. 14) and handed
to a gradient-boosted tree ensemble that learns to correct residual error.

They report it as one of their largest single wins: **RMSE 0.8334 → 0.6744 on
their independent CCLE test set, a 19.7% reduction** (their Table 7), and
describe the base model as lacking "sufficient accuracy without this
component."

Applied here to the two best models from the matrix. XGBoost hyperparameters
match the paper exactly (`n_estimators=100, max_depth=6, learning_rate=0.05,
subsample=0.8`). Our analogue of their 128-dim `f_ij` is the penultimate
activation of the prediction head — the last learned representation before the
scalar readout — captured via a forward hook, also 128-dim.

The refiner is fit on the **training fold only** and applied to val/test,
inheriting the same cell-line-grouped split as every other experiment.

## Results

| Base model | Config | Metric | Without XGBoost | With XGBoost | Δ |
|---|---|---|---|---|---|
| **CrossAttention** | GE+Proteomics, onehot | RMSE | **1.2442** | 1.2557 | **−0.0115 (−0.92%)** |
| | | MAE | 0.9336 | 0.9341 | −0.0005 |
| | | R² | 0.7883 | 0.7843 | −0.0040 |
| | | PCC | 0.8897 | 0.8898 | +0.0001 |
| | | AUC | 0.9191 | 0.9188 | −0.0003 |
| | | F1 | 0.8423 | 0.8456 | +0.0033 |
| **GNN-GCN** | tri-omics, fingerprint, +mutation edges | RMSE | **1.3305** | 1.3394 | **−0.0089 (−0.67%)** |
| | | MAE | 0.9970 | 1.0044 | −0.0074 |
| | | R² | 0.7681 | 0.7650 | −0.0031 |
| | | PCC | 0.8802 | 0.8785 | −0.0017 |
| | | AUC | 0.9116 | 0.9099 | −0.0017 |
| | | F1 | 0.8308 | 0.8306 | −0.0002 |

**Refinement made both models slightly worse.** The effect is small (<1%) but
consistent in direction across both architectures and across nearly every
metric — the opposite sign to the paper's +19.7%.

## Why the discrepancy: the split, not the method

The most likely explanation is a **methodological difference in how train/test
are separated**, not an implementation error.

| | MoGraphDRP | This project |
|---|---|---|
| Split | Random 80/10/10, "maintaining the distribution of IC50 values" (their §2.5) | `GroupShuffleSplit` **grouped by cell line**, 70/15/15 |
| Can a test cell line appear in training? | **Yes** | **No** (asserted in code) |

Under a random split, the same cell line appears in both training and test with
different drugs. A residual corrector can then learn **cell-line-specific
offsets** — "this cell line's IC50s run 0.3 high" — and apply them at test time,
because it has seen that cell line before. That is a legitimate and sizeable
source of correctable residual, and it would plausibly account for a large
share of a 19.7% gain.

Under our grouped split, every test cell line is entirely unseen. There is no
cell-line-specific residual structure available to transfer, so the refiner has
nothing systematic left to correct and only adds variance.

The refiner's own feature importances support this reading:

```
CrossAttention:  predicted_ic50 0.711,  f6 0.016,  f69 0.013,  f12 0.013, ...
GNN-GCN:         predicted_ic50 0.464,  f65 0.067, f45 0.051,  f1 0.047, ...
```

`predicted_ic50` dominates in both — the refiner is largely passing the base
prediction through, with the remaining thinly-spread interaction features
contributing corrections that do not generalize to held-out cell lines. Note
the paper's own Fig 10 shows the same dominance pattern (`predicted_ic50`
second-highest importance), consistent with the mechanism differing only in
whether those residuals transfer.

## What this means

1. **Do not adopt the ensemble stage.** It costs an extra model, an extra
   dependency, and a more complex inference path, for a consistent small loss
   under leakage-free evaluation.
2. **Treat cross-dataset comparisons with the paper cautiously.** Their headline
   RMSE 0.6622 is measured under a random split; our numbers are measured under
   a cell-line-grouped split, which is a strictly harder task (predicting for
   cell lines never seen in training — the actual clinical use case described
   in the proposal's §1 and §4). The two figures are **not directly
   comparable**. Quantified in
   [`split_protocol_comparison.md`](./split_protocol_comparison.md): protocol
   alone accounts for 0.347 RMSE (60% of the apparent gap); the genuine
   architectural difference is 0.235.
3. **The confirmatory experiment was run, and the hypothesis held.** Re-running
   this refinement under a random (non-grouped) split flips its sign:

   | Protocol | Base RMSE | + XGBoost | Δ |
   |---|---|---|---|
   | Grouped by cell line | 1.2442 | 1.2557 | **−0.92%** |
   | Random pairs | 0.8971 | 0.8885 | **+0.96%** |

   Identical refiner and hyperparameters; only the split differs. Under the
   random split 100% of test rows belong to a cell line seen ~208 times in
   training, and the refiner starts helping. Under the grouped split there is
   no such structure and it only adds variance. Full analysis in
   [`split_protocol_comparison.md`](./split_protocol_comparison.md).

   Note the magnitude still falls short of the paper's +19.7%, so leakage is a
   *necessary* condition for the gain rather than a complete account of it —
   their bilinear-attention interaction vector likely carries more correctable
   structure than our head's penultimate activation.

## Caveats

- **The GNN base number differs slightly from its matrix run** (1.3305 here vs
  1.3513 in [`gnn_ablation_results.md`](./gnn_ablation_results.md)). Both are
  retrained from the same seed and split, so the gap reflects GPU
  non-determinism in the conv kernels. The before/after comparison in this
  document is internally consistent — both columns come from the *same* trained
  model — so the conclusion is unaffected.
- Only the two best models were refined, per plan. It is possible refinement
  helps weaker base models (more residual left to correct), but that would be
  an argument for refinement as a crutch, not as an improvement to the final
  architecture.
- Our interaction vector is the head's penultimate activation, which is an
  interpretation of the paper's `f_ij` (theirs comes out of a bilinear
  attention module we did not implement). A different choice of extraction
  point might behave differently, though the dominance of `predicted_ic50` in
  the importances suggests the specific vector matters less than the absence of
  transferable residual structure.
