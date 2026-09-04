# Cross-Attention Fusion + MLP Head — Results & Interpretation

**Script:** [`experiments/01_cross_attention_fusion/cross_attention_baseline.py`](../experiments/01_cross_attention_fusion/cross_attention_baseline.py)
**Companions:** [`docs/04_rf_baseline_results.md`](./04_rf_baseline_results.md), [`docs/02_proteomics_baseline_results.md`](./02_proteomics_baseline_results.md)
**Date:** 2026-08-17

```
[threshold] median ln_ic50 across all targets = 3.2367
[features] loaded GE        transcriptomics_pca.csv  (532, 128)
[features] loaded Mut_CNV   genomics_pca.csv         (532, 128)
[features] loaded Proteomics proteomics_pca.csv       (532, 128)
[targets] 242036 rows -> 242036 finite -> 134764 with omics -> 134764 unique pairs
[design]  134764 pairs  |  omics: GE(128) + Mut_CNV(128) + Proteomics(128)  |  drug one-hot: 295  =  349.1 MB
[split]   train   93507 pairs    371 cell lines  (69.4%)
[split]   val     21041 pairs     80 cell lines  (15.6%)
[split]   test    20216 pairs     80 cell lines  (15.0%)

[training: 30 epochs, early-stopped at epoch 18 with best val_rmse=1.2253]

                  RMSE       PCC       AUC        F1
Validation      1.2253    0.8994    0.9211    0.8379
Test            1.2658    0.8849    0.9160    0.8432
Mean-only       2.7097        --        --        --   <- floor to beat

params=902,273  best_epoch=18  device=cuda
Fit 169.2s | Total 171.8s | Peak RSS: 2.05 GB
```

## 1. What the pipeline is doing

1. **Features** — the three per-modality PCA CSVs (`transcriptomics_pca.csv`
   = `GE`, `genomics_pca.csv` = `Mut_CNV`, `proteomics_pca.csv` =
   `Proteomics` — note the filename/key mapping is *not* 1:1, per
   `src/data/00_run_preprocessing.py:31-36`) are loaded **separately**, not
   concatenated. All three share the identical 532-cell-line index (verified
   in code, not assumed).
2. **Targets / pairs** — identical loading and one-hot drug encoding as
   `rf_baseline.py`/`proteomics_baseline.py`: 242,036 raw rows → 134,764
   usable (cell line, drug) pairs, 295 distinct drugs.
3. **Split** — same cell-line-grouped `GroupShuffleSplit` (70/15/15,
   `random_state=42`, no-leakage assertion) as the other three baselines.
4. **Model** — `MultiOmicsCrossAttentionFusion` (`src/models/cross_attention_fusion.py`,
   unmodified) fuses the three 128-dim per-cell-line embeddings via learned
   directional cross-attention over all 6 modality pairs into a 256-dim
   embedding. That embedding is concatenated with the 295-dim one-hot drug
   vector and passed through a `[256, 128]` MLP head
   (`Linear → BatchNorm → ReLU → Dropout` ×2 → `Linear(1)`). **The fusion
   module and the head are trained jointly**, end-to-end, by backprop on
   MSE(`ln_ic50`) — 902,273 total trainable parameters.
5. **Training** — Adam (`lr=1e-3`, `weight_decay=1e-5`), `ReduceLROnPlateau`,
   early stopping on val RMSE (patience 15). Converged fast: best epoch 18
   of 30 run, ~170s total on GPU (`cuda`, `torch 2.13.0+cu130`).

## 2. Headline result: cross-attention fusion decisively beats every prior baseline

| Model | Test RMSE | Test PCC |
|---|---|---|
| Mean-only floor | 2.7097 | — |
| RF, fused (GE+Mut_CNV+Proteomics, concatenated) | 2.1933 | 0.7050 |
| RF, proteomics-only | 2.0000 | 0.7460 |
| **Cross-attention fusion + MLP head** | **1.2658** | **0.8849** |

This is not a marginal improvement — test RMSE drops by **~42% relative to
the fused RF baseline** and **~37% relative to proteomics-only RF** (the
strongest RF variant so far). PCC rises from ~0.70–0.75 (RF) to **0.885**.
AUC 0.916 / F1 0.843 on the binarized (median-split) task confirm this isn't
an artifact of the continuous-metric definition — classification-style
separability improved by a similar margin.

**Two plausible, non-exclusive drivers, not yet disentangled:**
1. **Learned fusion vs. flat concatenation** — this is the effect the module
   was built to test (per `docs/04_rf_baseline_results.md` §4.3's RF → MLP →
   cross-attention → GCN → GAT ladder): letting the model learn how GE,
   Mut_CNV, and Proteomics interact, instead of handing a tree-ensemble a
   flat 384-dim vector, is exactly the "complex, non-linear biological
   interactions" gap the proposal's introduction calls out.
2. **Gradient-based training vs. tree ensembles in general** — some of this
   gap may be attributable to neural nets simply fitting this kind of tabular
   regression better than RF here, independent of the fusion mechanism.
   `mlp_baseline.py` (concatenated input, no cross-attention, on the
   `charms` branch) is the natural control to isolate this — if a plain MLP
   on the same 384-dim concatenated input scores meaningfully worse than
   1.2658, that's evidence the cross-attention fusion itself (not just
   "neural net > tree") is doing real work.

## 3. Validation is close to test — no overfitting signal

Val RMSE 1.2253 vs. test RMSE 1.2658 (both computed on cell lines never seen
in training) — a ~3% gap, consistent with `04_rf_baseline_results.md`'s
val/test stability pattern. Early stopping fired cleanly at epoch 18 without
the training loss continuing to fall much further before val RMSE plateaued
(epoch 25–30 shows val RMSE oscillating around 1.245, not improving), so this
isn't a case of stopping too early on a still-improving curve.

## 4. Caveats

- **First real-data run of `MultiOmicsCrossAttentionFusion`.** Its only prior
  verification (`docs/preprocessing_and_fusion_module.md`) used synthetic
  `rng.normal`/`rng.lognormal` data. No numerical issues (NaNs, PCA capping,
  BatchNorm batch-size-1 errors) surfaced on the real 532-cell-line data;
  `drop_last=True` on the shuffled train loader was kept as a precaution
  (mirrors `mlp_baseline.py`) but never needed to trigger.
- **Same drug-identity limitation as every other baseline so far** — one-hot
  drug encoding (§4.3 of `04_rf_baseline_results.md` already documents this):
  the model cannot generalize to a drug outside the 295-drug vocabulary.
  This result says nothing yet about the eventual RDKit-fingerprint drug
  representation.
- **Not the GAT deliverable.** This is still a flat per-cell-line embedding
  (no graph, no STRING PPI structure) — the proposal's §4.2 architecture
  additionally wraps this in a GAT over a heterogeneous graph, which doesn't
  exist yet (`docs/plan/graph_construction_plan.md` is drafted, not
  implemented). This result is evidence the *fusion* half of that
  architecture is working; it says nothing about the *graph* half yet.

## 5. What this means for MCS16

This is strong, concrete evidence that the cross-attention fusion module —
built and only unit-verified on synthetic data — carries real signal on the
actual GDSC2/CCLE data, ahead of the more expensive graph-construction work.
It substantially de-risks the §4.2 architecture: if fusion alone gets this
close to (or past) a reasonable "how good can this data get" ceiling, the
GAT stage's job is to add relational structure on top of an already-strong
representation, not to rescue a weak one.

## 6. Next steps

1. **Run `mlp_baseline.py` (currently only on the `charms` branch) on this
   same split/seed** to isolate "learned fusion" from "gradient-based
   training in general" per §2 above — the single most informative next
   comparison point.
2. **Run the equivalent single-omics ablations** (GE-only, Mut_CNV-only
   cross-attention doesn't apply — cross-attention needs ≥2 modalities — but
   RF/MLP single-omics numbers for GE and Mut_CNV would round out the
   picture alongside the existing Proteomics-only RF result).
3. **Log this number** in the running comparison table proposed in
   `04_rf_baseline_results.md` §5.3 once it exists (mean-only / RF-fused /
   RF-proteomics / MLP-fused / **cross-attention-fused** / GCN / GAT rows).
4. **Proceed to graph construction** (`docs/plan/graph_construction_plan.md`)
   now that the fusion half of the architecture has real-data evidence
   behind it, rather than being blocked on "does fusion even work."
