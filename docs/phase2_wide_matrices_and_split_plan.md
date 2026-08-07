# Phase 2: Wide Matrices — What Was Done, and the Train/Val/Test Split Plan

## What was done

### 1. Fixed a memory inefficiency before running at real scale

`build_wide_matrices.py`'s row filter used `.astype(str) != ""` on columns
loaded as `dtype="category"`. That materializes a full object-dtype copy of
the column, which defeats the entire point of loading it as `category` in
the first place. With only ~9.6GB RAM available locally and 78.8M real
RNA-seq rows to process, this was flagged as a real OOM risk before running
(not theoretical — the run measured 2GB of swap usage even *after* the fix).
Replaced with `dropna()` + a direct comparison against the categorical,
which never leaves the compact categorical representation.

### 2. Ran `build_wide_matrices.py` against the real aligned output

```
.venv/bin/python src/data/build_wide_matrices.py
```

Exit code 0. Real results:

| Matrix | Shape (cell lines × features) | Size |
|---|---|---|
| `GE_wide.csv` | 1,898 × 41,143 | 418M |
| `Mut_CNV_wide.csv` | 1,570 × 21,582 (20,797 mutation + 785 CNV) | 136M |
| `Proteomics_wide.csv` | 948 × 8,453 | 51M |
| **Cell lines common to all 3** | **919** | — |

This is the first real validation of the RNA-seq pivot path specifically —
earlier validation only covered real CNV, real proteomics, and a 2M-row real
mutation subset (see `docs/preprocessing_and_fusion_module.md`).

### 3. Found (and mapped a fix for) a second instance of the ID-space bug

The `depmap_id` (`ACH-*`) vs `sanger_model_id` (`SIDM*`) mismatch documented
in `docs/data_ingestion_and_alignment.md` — originally found in the
diagnostic `omics_availability_by_model.csv` columns — turns out to matter
for something load-bearing, not just diagnostics: **`gdsc2_response_master.csv`
uses `ACH-*` as `standard_model_id`, but the wide matrices built in this
phase are indexed by `SIDM*`** (inherited from `rnaseq_aligned.csv` /
`cnv_aligned.csv` / `mutations_aligned.csv`, which use the raw CCLE file's
own `model_id` column as-is).

Verified directly: checked the first 50,000 GDSC response rows against
`GE_wide.csv`'s index — **1 accidental overlap**. A direct join between GDSC
labels and the wide omics matrices silently fails almost entirely.

**Fix (planned, not yet applied):** `model_crosswalk.csv` has 1,771 rows
with both IDs present, enough to map `standard_model_id` (`ACH-*`) →
`sanger_model_id` (`SIDM*`) before joining GDSC labels to the wide matrices.
Verified this resolves it: joining through the crosswalk instead of directly
recovers **228,654 usable (cell line, drug) rows** with a valid `ln_ic50`,
across **913** of the 919 common cell lines (250.4 rows/cell line on
average) — this is the real, usable labeled dataset size for training.

## Train/val/test split — options evaluated

The proposal report (§4.1.2) specifies: *"the dataset is stratified by cell
line and split into training (70%), validation (15%), and test (15%) sets."*
That single sentence under-specifies a real decision, so here's what was
actually considered, using the real 913-cell-line / 228,654-row dataset
above (42 distinct cancer types present, smallest has 4 cell lines, largest
82 — full distribution below).

| # | Option | Leakage-free? | Tests generalization to | Verdict |
|---|---|---|---|---|
| 1 | Random split by (cell line, drug) **row** | ❌ No | Nothing new — same cell lines seen in train & test | **Reject.** A cell line's omics profile would appear in training (with one drug) and test (with another). The fusion module can learn cell-line-specific baseline sensitivity during training and "recognize" that cell line at test time regardless of the drug, inflating metrics without proving the model generalizes to a genuinely new patient. |
| 2 | Group split by **cell line only** (no cell line's rows split across sets) | ✅ Yes | Unseen cell lines | **Matches the proposal literally.** Solid default, but doesn't control for cancer-type imbalance across splits — by chance, a rare cancer type (e.g. the 4 Biliary Tract Carcinoma lines) could land entirely in train, meaning the model is *never* evaluated on it. |
| 3 | **Group split by cell line, stratified by `cancer_type`** | ✅ Yes | Unseen cell lines, with balanced representation per cancer type in every split | **Recommended.** Strictly improves on Option 2 at the same leakage guarantee — same no-leakage property, more statistically reliable test metric, and enables a per-cancer-type performance breakdown later (useful groundwork for the SHAP interpretability module). |
| 4 | Group split by **drug** (leave-drugs-out) | ✅ Yes (drug axis) | Unseen drugs, not unseen patients | **Reject as primary.** Answers a different clinical question (drug repurposing/screening) than this project's stated framing (predicting response for a new patient's cell line — §1's "personalized medicine" motivation, and §4's inference flow where a researcher uploads a *cell line*). Worth keeping as a **secondary ablation** later, alongside the GCN-vs-GAT and single-modality ablations already planned in §4.3. |
| 5 | Two-way disjoint (cell line **and** drug both held out) | ✅ Yes (both axes) | The strictest "cold start" case | **Reject as primary.** Shrinks usable test data to the intersection of two held-out sets, and isn't what any of the reviewed benchmarks (NIHGCN, MoGraphDRP, DeepMoDRP) use as their primary metric. Interesting future robustness check, not a blocker now. |

### Recommended implementation (Option 3)

Split happens at the **cell-line grain**, not the row grain — each of the
913 labeled cell lines has exactly one `cancer_type`, so this doesn't need
`GroupKFold`/`StratifiedGroupKFold` machinery; a plain two-stage
`sklearn.model_selection.train_test_split` on the list of unique cell lines,
stratified by `cancer_type`, achieves the same thing more simply:

```python
train_ids, temp_ids = train_test_split(
    cell_line_ids, test_size=0.30, stratify=cancer_type_per_cell_line, random_state=42
)
val_ids, test_ids = train_test_split(
    temp_ids, test_size=0.50, stratify=cancer_type_per_temp_cell_line, random_state=42
)
```

Every (cell line, drug) row then inherits its split from its cell line —
giving ~70/15/15 by cell line count (~639/137/137), which at 250 rows/cell
line average works out to roughly 160K/34K/34K labeled pairs per split.

**Known edge case to handle when implementing:** `train_test_split`'s
`stratify` requires at least 2 members per class per split. The smallest
cancer types (4–7 cell lines: Biliary Tract Carcinoma, Prostate Carcinoma,
Non-Cancerous, Chondrosarcoma, Hodgkin's Lymphoma, Mesothelioma,
T-Cell Non-Hodgkin's Lymphoma, Other Blood Cancers, Rhabdomyosarcoma) will
likely error out of a naive two-stage stratified split, since a 4-member
class split 70/15/15 can leave a side with fewer than 2 members. Resolution:
collapse cancer types below a small threshold (e.g. < 10 cell lines) into a
shared `"Other/Rare"` stratification bucket *for splitting purposes only* —
the real `cancer_type` label is untouched everywhere else (features, eval
breakdowns, etc.), it's purely a workaround for the stratifier's minimum-count
requirement.

### Where this fits in the pipeline

The split must happen **before** `OmicsPreprocessingPipeline.fit_transform()`,
not after: `fit_transform()` should only ever see `train_ids`' omics rows, so
the fitted scalers/PCA don't leak information from val/test cell lines.
`val_ids`/`test_ids` go through `.transform()` using the pipeline already
fitted (and saved) on `train_ids` alone.

## What's next

1. Apply the crosswalk-based ID fix and implement the split above (with the
   rare-class bucketing) as a real script/module — not built yet, this
   document is the plan, not the code.
2. Run `OmicsPreprocessingPipeline.fit_transform()` on `train_ids` only, then
   `.transform()` on `val_ids`/`test_ids`, against the real wide matrices —
   first real (non-synthetic) end-to-end preprocessing run.
3. From there: drug encoding (RDKit) and STRING PPI graph construction
   (Phase 4 in the earlier runbook), which are unblocked by this split since
   they don't depend on it.

## Cancer type distribution (for reference)

913 cell lines, 42 types, sorted ascending:

```
Biliary Tract Carcinoma: 4        Osteosarcoma: 11
Prostate Carcinoma: 5              Burkitt's Lymphoma: 12
Non-Cancerous: 6                   Cervical Carcinoma: 13
Chondrosarcoma: 6                  Squamous Cell Lung Carcinoma: 13
Hodgkin's Lymphoma: 6              Glioma: 14
Mesothelioma: 7                    Other Sarcomas: 15
T-Cell Non-Hodgkin's Lymphoma: 7   Plasma Cell Myeloma: 15
Other Blood Cancers: 7             Thyroid Gland Carcinoma: 15
Rhabdomyosarcoma: 7                Hepatocellular Carcinoma: 15
Esophageal Carcinoma: 9            Head and Neck Carcinoma: 16
Chronic Myelogenous Leukemia: 10   T-Lymphoblastic Leukemia: 17
Endometrial Carcinoma: 10          Bladder Carcinoma: 18
B-Lymphoblastic Leukemia: 19       Kidney Carcinoma: 32
Oral Cavity Carcinoma: 22          Glioblastoma: 34
Ewing's Sarcoma: 23                B-Cell Non-Hodgkin's Lymphoma: 35
Esophageal Squamous Cell Carcinoma: 26   Ovarian Carcinoma: 40
Acute Myeloid Leukemia: 27         Colorectal Carcinoma: 46
Gastric Carcinoma: 27              Breast Carcinoma: 49
Other Solid Cancers: 28            Melanoma: 51
Pancreatic Carcinoma: 29           Small Cell Lung Carcinoma: 55
Neuroblastoma: 30                  Non-Small Cell Lung Carcinoma: 82
```
