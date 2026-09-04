# Random Forest Baseline — Results & Interpretation

**Script:** [`src/data_engineering/models/rf_baseline.py`](../src/data_engineering/models/rf_baseline.py)
**Run:** identical output reproduced independently on every team member's machine
**Date:** 2026-08-12
**Proposal reference:** [`MCS16 - Final Proposal Report.docx`](./MCS16%20-%20Final%20Proposal%20Report.docx)

```
[features] loaded cache fused_early.npy (532, 384)
[targets] 242036 rows -> 242036 finite -> 134764 with omics -> 134764 unique pairs
[design]  134764 pairs x 679 features (384 omics + 295 drug one-hot) = 349.1 MB
[split]   train   93507 pairs    371 cell lines  (69.4%)
[split]   val     21041 pairs     80 cell lines  (15.6%)
[split]   test    20216 pairs     80 cell lines  (15.0%)

                  RMSE       PCC
Validation      2.2556    0.7245
Test            2.1933    0.7050
Mean-only       2.7097        --   <- floor to beat
```

## 1. What the pipeline is doing

1. **Features** — three omics modalities (genomics, transcriptomics, proteomics) were each z-scored and PCA-compressed, then early-fused by concatenation (`load_and_fuse.py`) into a single `532 cell lines × 384 features` matrix, cached as `fused_early.npy`.
2. **Targets** — GDSC2 drug-response records (`ln_ic50`) are filtered to finite values, restricted to cell lines that have omics coverage, and de-duplicated per (cell line, drug) pair. 242,036 raw rows collapse to 134,764 usable pairs.
3. **Design matrix** — each pair gets the 384-dim omics vector for its cell line, concatenated with a 295-dim one-hot drug identity vector → 679 features × 134,764 rows (349 MB as float32).
4. **Split** — a `GroupShuffleSplit` groups by **cell line**, not by pair, so no cell line's omics profile appears in both train and test (train/val/test ≈ 70/15/15 by pairs, no group overlap — asserted in code).
5. **Model** — a 100-tree Random Forest (`min_samples_leaf=5`, `max_features="sqrt"`) regresses `ln_ic50` on the 679 features.

## 2. What the results mean

- **Mean-only floor (RMSE 2.7097):** the error from always predicting the training mean — a naive baseline with zero learned signal.
- **RF Test RMSE 2.1933 / PCC 0.7050:** the model cuts error by ~19% versus the naive floor and produces predictions that correlate strongly (r ≈ 0.70) with true `ln_ic50` on cell lines it never saw during training.
- **Val vs. test are close** (RMSE 2.26 vs 2.19, PCC 0.72 vs 0.71), which means the model isn't overfitting to the validation split — performance is stable across two independent held-out groups of cell lines.

**Conclusion:** the PCA-compressed, early-fused omics features carry real, learnable signal about drug response — this is not a spurious fit. The data engineering pipeline (alignment, PCA, fusion) upstream of this script is producing a genuinely usable feature set. A PCC around 0.70 for a tabular baseline on GDSC-scale IC50 prediction is a credible, non-trivial result, not a floor that's trivially easy to beat.

## 3. Why every machine produced identical numbers

This is expected, not a coincidence, and it's a good sign for reproducibility. Every source of randomness in the script is pinned:

| Source of randomness | How it's controlled |
|---|---|
| Train/val/test split | Both `GroupShuffleSplit` calls use `random_state=RANDOM_STATE` (42) |
| Random Forest bootstrap sampling & feature subsampling | `RandomForestRegressor(..., random_state=42)` |
| Parallel tree fitting (`n_jobs=-1`) | scikit-learn draws one RNG seed **per tree** from the master `random_state` *before* dispatching to worker threads/processes — thread scheduling never feeds back into which data/features a tree sees, so results are identical whether it runs on 4 cores or 32 |
| Input data | The script reads a pre-built cache (`fused_early.npy`, `fused_cell_lines.csv`, `gdsc2_response_master.csv`) rather than regenerating it — as long as everyone is running against the same cached files, the design matrix is byte-identical on every machine |
| No GPU / no stochastic layers | RF training is pure CPU with no dropout, augmentation, or non-deterministic reduction order — nothing hardware-dependent influences the *result* |

**What legitimately does vary by machine** (and isn't a concern): `Fit`/`Predict` wall-clock time and `Peak RSS` are hardware-dependent — different CPUs/core counts will show different timings even though `RMSE`/`PCC` stay identical.

**One fragility worth flagging:** [`requirements.txt`](../requirements.txt) currently only pins GUI/API packages (PySide6, fastapi, uvicorn) — it does not pin `numpy`/`pandas`/`scikit-learn` versions. The run happened to be identical across everyone's current environments, but a future scikit-learn release *could* change RF's internal tie-breaking or split-finding behavior and silently break this reproducibility. Recommend pinning `numpy`, `pandas`, and `scikit-learn` (exact or `~=` versions) once the baseline numbers are considered "final," so they stay citable in the report.

## 4. What this means for MCS16 & data meaning

### 4.1 This isn't a generic baseline — it's the literal metric named in the proposal

Section 3.1.1 (*Product Characteristics & Requirements*) of the Final Proposal Report defines the project's core success criterion in one sentence:

> "the platform must achieve a lower Root Mean Squared Error (RMSE) than **standard flat vector concatenation baselines**"

That is exactly what this script is: a flat-vector concatenation model (PCA'd omics concatenated with a one-hot drug vector, fed to a non-graph learner). Section 2.2.1 of the literature review describes this same category directly — *"Traditional machine learning approaches, including support vector regression and random forest... established reproducible baselines for DRP but treated molecular features as flat, independent inputs."* This run is that sentence, executed on our own data. **RMSE 2.1933 / PCC 0.7050 (test) is the number the eventual cross-attention + GAT model has to beat** to satisfy the proposal's own stated success metric — not an arbitrary internal sanity check.

### 4.2 It also confirms Sprint S1–S2 deliverables are real, not just planned

Section 4.1.2 (*Preprocessing Pipeline*) specifies: per-modality z-score normalization, PCA compression to a shared dimension, and a cell-line-stratified 70/15/15 train/val/test split "to prevent data leakage." `load_and_fuse.py` and `rf_baseline.py` implement precisely that spec — this is the Data Engineering track (Sprints S1–S2, per §3.3) already working end-to-end, not just designed on paper. Today (2026-08-12) sits inside that window, ahead of the **Data Freeze milestone on 24 August 2026** — this baseline is effectively the checkpoint that proves the pipeline will be ready in time for that freeze.

### 4.3 Where this baseline diverges from the target architecture — and why that matters

The final MCS16 model (§4.2) is not this script scaled up; it replaces two specific pieces this baseline uses as placeholders:

| Component | This RF baseline | Target MCS16 architecture (§4.1.2, §4.2.1, §4.2.2) |
|---|---|---|
| Drug representation | 295-dim one-hot identity vector (a drug is just an index) | 2048-bit RDKit Morgan fingerprint from SMILES — encodes actual chemical structure |
| Omics fusion | Simple concatenation (early fusion) | Cross-attention fusion — Q/K/V attention across all 6 modality pairs |
| Relational structure | None — flat feature vector per pair | Graph Attention Network over a heterogeneous graph incorporating STRING PPI edges |
| Generalization to new drugs | Impossible — a drug absent from the one-hot vocabulary has no representation | Possible in principle — fingerprints are computed from structure, not identity |

This matters for **data meaning**, not just modeling: the one-hot drug encoding here means the model can only ever say "which of these 295 known drugs does this look like," never "what does this molecule's structure imply." The PCC ≈ 0.70 we're seeing is entirely attributable to the omics side and drug identity — it contains zero chemical-structure signal. That's a meaningfully different (and easier) problem than what §4.3's planned ablation study is designed to test, so this number should be read as "the floor with drug-identity-only," not yet "the floor with drug-structure information," which is the actual baseline the ablation study in §4.3 calls for.

### 4.4 A real coverage gap this run surfaces

GDSC2 has drug-response records for **969** distinct cell lines, but `fused_early.npy` only has omics for **532** of them (≈55%) — the rest were dropped before this script ever ran, during modality alignment (a cell line needs genomics **and** transcriptomics **and** proteomics coverage to survive the fusion step). Section 2.2.5 of the literature review calls out proteomics under-representation as a field-wide gap; this pipeline is now hitting that same wall concretely — proteomics coverage in CCLE is well known to be far sparser than transcriptomics/CNV/mutation coverage, and this ~45% cell-line drop is consistent with that being the binding constraint here. The proposal's own conclusion (§5) already acknowledges the related limitation that predictions rely on cell line data rather than clinical data and exclude some modalities — this is the same class of constraint, now visible as a concrete number rather than a caveat in prose. It means every metric in this document (and every future GNN metric trained on the same fused matrix) describes performance on the ~55% of cell lines with complete tri-omics coverage, not the full GDSC2 cohort.

## 5. Next steps

1. **Decide the "official" concatenation baseline before Data Freeze (24 Aug 2026).** Re-run this script with the drug block swapped from one-hot to RDKit Morgan fingerprints (the same drug encoding the final model will use, per §4.1.2) so the number reported against §3.1.1's success criterion and used in §4.3's ablation table is on the same drug representation as the GAT model, not an easier one-hot version.
2. **Pin `numpy`, `pandas`, and `scikit-learn` in `requirements.txt`** now, while these numbers are being treated as citable in the report — closes the reproducibility fragility noted in Section 3 before more results accumulate on unpinned versions.
3. **Log this baseline number durably**, e.g. a running `docs/results.md` table (mean-only / RF-one-hot / RF-fingerprint / GCN / GAT rows), so the "floor to beat" from §3.1.1 has a permanent, auditable home instead of only living in a terminal transcript.
4. **Flag the 532-vs-969 cell-line coverage gap to the team** (Aditya, as Bioinformatics Data Engineer per §3.2.2) — decide whether to attempt recovering cell lines with partial omics (e.g. imputation) or formally document the ~55% coverage as a scope boundary in §3.1.2 before graph construction begins, since the GAT model will inherit the same cap unless it's addressed.
5. **Re-run the same `GroupShuffleSplit`/`random_state=42` methodology when evaluating GCN and GAT variants** in the §4.3 ablation study, so RF vs. GCN vs. GAT comparisons in the final report sit on an identical train/val/test partition rather than being confounded by different splits.
