# Proteomics-Only Baseline — Results & Interpretation

**Script:** [`experiments/02_proteomics/proteomics_baseline.py`](../experiments/02_proteomics/proteomics_baseline.py)
**Companion:** [`docs/04_rf_baseline_results.md`](./04_rf_baseline_results.md) (the fused GE+Mut_CNV+Proteomics baseline this run is compared against)
**Date:** 2026-08-17

```
[threshold] median ln_ic50 across all targets = 3.2367
[features] loaded proteomics_pca.csv (532, 128)
[targets] 242036 rows -> 242036 finite -> 134764 with omics -> 134764 unique pairs
[design]  134764 pairs x 423 features (128 omics + 295 drug one-hot) = 217.5 MB
[split]   train   93507 pairs    371 cell lines  (69.4%)
[split]   val     21041 pairs     80 cell lines  (15.6%)
[split]   test    20216 pairs     80 cell lines  (15.0%)

                  RMSE       PCC       AUC        F1
Validation      2.0698    0.7530    0.7723    0.6438
Test            2.0000    0.7460    0.7520    0.6435
Mean-only       2.7097        --        --        --   <- floor to beat
```

## 1. What the pipeline is doing

1. **Features** — a single omics modality (proteomics) loaded directly from
   `data/processed/proteomics_pca.csv` (already z-scored, median-imputed, and
   PCA-compressed to 128 dims by `OmicsPreprocessingPipeline`) — no early
   fusion, no other modality involved.
2. **Targets** — same GDSC2 `ln_ic50` loading/filtering as the fused baseline:
   242,036 raw rows → 134,764 usable (cell line, drug) pairs.
3. **Design matrix** — each pair gets the 128-dim proteomics vector for its
   cell line, concatenated with the same 295-dim one-hot drug identity vector
   used in the fused baseline → 423 features × 134,764 rows.
4. **Split** — identical cell-line-grouped `GroupShuffleSplit` methodology as
   `rf_baseline.py` (70/15/15, no cell-line leakage, asserted in code), but
   with `random_state=42` pinned (the fused baseline left this `None`), so
   this split is reproducible run-to-run.
5. **Model** — the same 100-tree Random Forest (`min_samples_leaf=5`,
   `max_features="sqrt"`) as the fused baseline, regressing `ln_ic50`.
6. **Metrics** — RMSE/PCC as before, plus AUC/F1 on the target binarized at
   the median `ln_ic50` (3.2367) across all targets — a pattern carried over
   from the team's Transcriptomics/Genomics single-omics script so all
   single-omics runs report on the same four metrics.

## 2. A real surprise: proteomics-only beats the fused baseline

| | RMSE (test) | PCC (test) |
|---|---|---|
| Fused (GE + Mut_CNV + Proteomics, 384-dim) | 2.1933 | 0.7050 |
| **Proteomics only (128-dim)** | **2.0000** | **0.7460** |
| Mean-only floor | 2.7097 | — |

Proteomics alone cuts RMSE by **~26%** vs. the mean-only floor (vs. ~19% for
the fused model) and correlates more strongly with true `ln_ic50` (PCC 0.746
vs. 0.705). This is not a coverage artifact: `proteomics_pca.csv` is indexed
on the **same 532 cell lines** as `fused_early.npy` (both are downstream of
`OmicsPreprocessingPipeline.align_by_cell_line()`, which intersects all three
modalities before fitting), and `[targets]` resolves to the identical 134,764
pairs in both runs. So this is a like-for-like comparison on the same rows —
adding GE and Mut_CNV features on top of proteomics made the Random Forest
**worse**, not better.

**Plausible explanation, not yet confirmed:** with `max_features="sqrt"`, a
tree at each split samples `sqrt(423) ≈ 21` candidate features from the
423-dim proteomics design matrix vs. `sqrt(679) ≈ 26` from the 679-dim fused
one — but the fused matrix's 384 omics columns dilute the pool with 256 extra
GE/Mut_CNV dimensions competing against the 295 drug one-hot columns and 128
proteomics columns for split selection, so any given split is statistically
less likely to land on the drug identity or a genuinely predictive proteomics
feature. This is consistent with the literature review's note (quoted in
`04_rf_baseline_results.md`) that proteomics is under-represented and
comparatively information-dense per feature in CCLE-scale data. Worth a
follow-up ablation (GE-only, Mut_CNV-only) to confirm proteomics specifically
is the strongest single modality rather than "any single modality beats the
concatenation" being a general RF/`max_features` artifact.

## 3. Reproducibility

Unlike `rf_baseline.py` (`RANDOM_STATE = None`), this script pins
`RANDOM_STATE = 42` for both the split and the Random Forest, so re-running
it produces byte-identical `RMSE`/`PCC`/`AUC`/`F1` on any machine (same
caveat as `04_rf_baseline_results.md` §3 about unpinned `scikit-learn` versions
in `requirements.txt` still applies).

## 4. What this means for MCS16

Per `04_rf_baseline_results.md` §4.3, the RF baseline (fused or single-omics) is
a **flat, non-graph, non-relational** placeholder — this result doesn't
change the target architecture (cross-attention fusion + GAT per §4.2), but
it is directly useful input for it:

- It's a second, independently-produced data point (alongside `shaahid`'s
  Genomics run) for the "GCN/GAT single-modality ablation" already planned in
  the proposal's §4.3 — the eventual GAT model should be checked against
  per-modality RF numbers like this one, not just the fused RF number.
- It suggests proteomics carries a disproportionate amount of the fused
  model's signal, which is worth flagging to the team before the cross-attention
  fusion module is trained on real data — if this holds up, the attention
  weights between modality pairs involving proteomics are a specific thing to
  sanity-check once real (non-synthetic) `MultiOmicsCrossAttentionFusion`
  training happens.

## 5. Next steps

1. Run the equivalent single-omics ablation for GE (transcriptomics) — the
   `shaahid` branch already covers Genomics and Proteomics in
   `experiments/03_transcriptomics_genomics/TransGen_Baseline.py`; a
   Transcriptomics-only number would complete the three-way comparison and
   test whether "single modality beats fused" is proteomics-specific.
2. Re-run with `max_features` swapped to a fixed count (rather than `"sqrt"`)
   across both fused and single-omics runs, to check whether the RMSE gap in
   §2 is a `max_features`-driven artifact of dimensionality rather than a
   genuine proteomics-signal effect.
3. Log this number in the running `docs/results.md` comparison table proposed
   in `04_rf_baseline_results.md` §5.3, once that table exists.
