# Multi-Omics Preprocessing & Cross-Attention Fusion Module

This documents the code added to build **Part 1** of the MCS16 architecture: the
data preprocessing pipeline that reduces GE / Mut_CNV / Proteomics matrices to a
shared dimension, and the cross-attention module that fuses them into a single
cell-line embedding. It covers what was built, why, the design decisions made,
and what is intentionally still out of scope.

## Why this was built

An initial task brief (drafted with Gemini) proposed this module. Before
implementing it, it was checked against `docs/MCS16 - Final Proposal Report.docx`
and `PROJECT_CONTEXT.md` to confirm it actually matched the project's design —
several details in the brief either weren't specified anywhere, conflicted
between the two source docs, or contained a hard technical bug. Those were
resolved with the project owner before writing any code (see
[Decisions](#decisions-made-before-implementation) below).

## Files added

| File | Role |
|---|---|
| `src/data/omics_preprocessing.py` | `OmicsPreprocessingPipeline` — alignment, normalization, PCA compression |
| `src/models/cross_attention_fusion.py` | `PairwiseCrossAttention`, `MultiOmicsCrossAttentionFusion`, verification script |
| `requirements.txt` | Pins `torch`, `scikit-learn`, `pandas`, `numpy`, `joblib` |
| `src/__init__.py`, `src/data/__init__.py`, `src/models/__init__.py` | Make `src` an importable package (`python -m src.models.cross_attention_fusion`) |
| `.gitignore` | Added `models/` — fitted PCA/scaler artifacts shouldn't be committed |

Environment: `torch 2.13.0+cu130`, `scikit-learn 1.9.0`, `pandas 3.0.5`,
`numpy 2.5.1`, `joblib 1.5.3`, installed into the existing `.venv`. CUDA is
available in this environment.

## `OmicsPreprocessingPipeline` (`src/data/omics_preprocessing.py`)

Takes three wide DataFrames (cell line ID as index, features as columns) and
returns three PCA-compressed `(n_samples, d_target)` arrays, one per modality.

**Alignment** — `align_by_cell_line()` intersects the three DataFrames'
indices and reindexes all of them onto the same sorted set of cell line IDs,
so downstream tensors line up row-for-row across modalities.

**Per-modality preprocessing** (each modality gets its own `sklearn.Pipeline`,
fitted independently):

- **GE (transcriptomics):** `log1p` → `StandardScaler` → `PCA(d_target)`.
- **Mut_CNV (genomics):** `StandardScaler` → `PCA(d_target)` — mutations and
  CNV are treated as one combined matrix, matching the proposal report
  §4.1.2 rather than the binarize-mutations approach mentioned in
  `PROJECT_CONTEXT.md`; the two docs disagreed and this was the option chosen.
- **Proteomics:** median imputation (`SimpleImputer`) → `StandardScaler` →
  `PCA(d_target)`.

**PCA component safeguard** — `PCA` cannot produce more components than
`min(n_samples, n_features)`. The pipeline now caps `n_components` to that
bound automatically and prints a warning if it had to, instead of throwing.
This was a real bug in the original brief's own verification example (32
dummy samples, `d_target=128` — mathematically impossible), fixed here by
capping defensively *and* by raising the demo sample count (below).

**Persistence** — `save()`/`load()` serialize each modality's fitted
`Pipeline` (scaler + imputer + PCA weights) via `joblib` to `models/omics_pca/`,
so the exact training-time transform can be replayed on new cell lines at
inference without re-fitting (avoids train/inference skew).

## `PairwiseCrossAttention` / `MultiOmicsCrossAttentionFusion` (`src/models/cross_attention_fusion.py`)

- **`PairwiseCrossAttention`** — one directional multi-head attention block:
  `query_modality`, `key_value_modality` are each `(B, d_model)`. Internally
  reshaped to `(B, 1, d_model)` and passed through `nn.MultiheadAttention`,
  then residual-added and `LayerNorm`'d.
- **`MultiOmicsCrossAttentionFusion`** — builds one `PairwiseCrossAttention`
  per ordered pair, for all 6 directional pairs across
  `{GE, Mut_CNV, Proteomics}`. Each pair's output goes through a residual FFN
  (`Linear → GELU → Dropout → Linear`, then `LayerNorm`), all 6 are
  concatenated (`6 × d_model`), and projected through a final
  `Linear → GELU → Dropout(0.2) → LayerNorm` head to `(B, out_dim)`.

**Known nuance, not a bug:** since each modality is already pooled to a
single vector per sample before fusion, every attention call has exactly one
query and one key — softmax over a single element is always `1.0` by
construction. The representational power of this block comes from the
learned Q/K/V projections and the residual FFN stack, not from a non-trivial
attention weight distribution. Worth knowing if attention weights are later
used for interpretability (the proposal report's interpretability plan is
SHAP-based, not attention-based, so this likely doesn't block anything).

## Decisions made before implementation

The following were ambiguous or conflicting between the brief, the proposal
report, and `PROJECT_CONTEXT.md`, and were confirmed with the project owner
before writing code:

1. **Scope:** Build the full pipeline (preprocessing + fusion) in one pass,
   rather than splitting strictly along the AA/WH ownership lines in
   `PROJECT_CONTEXT.md`.
2. **Mutation encoding:** Combine mutations + CNV into one matrix, z-score,
   PCA — per the proposal report, not the binarize-then-separate-PCA
   approach noted in `PROJECT_CONTEXT.md`.
3. **Data scope:** Synthetic verification data only for this pass. Real CCLE/GDSC
   data was deliberately **not** wired up yet (see below).
4. **Setup:** Install `torch`/`scikit-learn`/`pandas`/`numpy` into `.venv`,
   and fix the PCA/sample-count bug by raising the synthetic demo to 200
   samples (`d_target=128` needs `n_samples > d_target`).

## Verified

Ran `python -m src.models.cross_attention_fusion` end to end:

```
Input shapes:
  GE: (200, 10000)
  Mut_CNV: (200, 5000)
  Proteomics: (200, 2000)

Post-PCA (intermediate pair) shapes:
  GE: (200, 128)
  Mut_CNV: (200, 128)
  Proteomics: (200, 128)

Final fused cell-line embedding shape: (200, 256)

MultiOmicsCrossAttentionFusion parameter count: 727,168 (trainable: 727,168)
Approx. parameter memory (fp32): 2.77 MB
```

Also separately verified the `save()`/`load()` round trip: fitted pipelines
persisted to `models/omics_pca/`, reloaded, and used to `transform()` a fresh
batch of 50 synthetic cell lines without re-fitting — confirms the
train/inference split works as intended.

727K parameters (~2.8MB fp32) is a small fraction of the 4–8GB VRAM budget —
comfortable headroom is left for the downstream GAT and MLP prediction head.

## Relationship to `src/data/ingest_and_align.py`

`ingest_and_align.py` and `omics_preprocessing.py` are sequential pipeline
stages, not alternatives — see the conversation for the full reasoning.
Short version:

- `ingest_and_align.py` turns **raw** GDSC/CCLE/STRING files into ID-crosswalked,
  **long**-format aligned CSVs (`data/processed/aligned/`). It solves ID
  matching across `model_id` / `BROAD_ID` / `CCLE_ID` / `COSMIC_ID`, which
  `omics_preprocessing.py` has no logic for at all.
- `omics_preprocessing.py` takes an already-aligned **wide** matrix
  (cell line × feature) and does normalization + PCA. It has no ID-matching
  logic and can't parse the raw `.zip`/`.xlsx` files.

**Recommendation: keep both.** `ingest_and_align.py` is not overkill — it
implements the ID-matching step explicitly called out in `PROJECT_CONTEXT.md`.
It could be simplified now that `pandas` is installed (its hand-rolled
XLSX/zip parsing predates having any dependencies available), but that's a
refactor, not a deletion.

## Fix: `ingest_and_align.py` was buffering unbounded row lists in memory

While preparing to run the script against the real CCLE/GDSC data (ahead of a
Colab run), a blocking bug surfaced: `load_omics_table` collected every parsed
row into a Python list before writing any of it to disk. The RNA-seq archive
alone unpacks to **78.8M rows** — buffering that as a list of dicts measured
out to **~69GB of RAM** (verified with `tracemalloc` on a 1M-row sample and
extrapolated), which would `MemoryError` on any local machine or Colab runtime
before producing a single output file.

Fixed by turning `load_omics_table` into `stream_omics_table`, a generator
that yields rows one at a time straight into `write_csv` (which already wrote
incrementally), with `counts_by_model` updated as a side effect during
iteration instead of a second pass over a materialized list. `main()` was
reordered so the rnaseq/cnv/mutation streamed writes happen first (populating
the counts), then `availability_rows`/`response_master_rows` are computed from
those counts afterward. Verified with `tracemalloc` that memory now stays flat
(under 1MB) regardless of row count, and cross-checked the full real run in
the background (see below).

## New: `src/data/build_wide_matrices.py` — the long → wide pivot + proteomics loader

Closes the gap noted below by bridging `ingest_and_align.py`'s long-format
output to the wide format `OmicsPreprocessingPipeline` expects, without
modifying `ingest_and_align.py`'s ID-matching logic at all:

- **`pivot_long_to_wide()`** — pivots `rnaseq_aligned.csv` (→ GE, on
  `rsem_tpm`) and `cnv_aligned.csv` (→ CNV, on `total_copy_number`) from one
  row per (cell line, gene) into one row per cell line. Duplicate
  (cell line, gene) rows are averaged (`aggfunc="mean"`) rather than erroring.
- **`pivot_mutation_presence()`** — mutations can't be pivoted the same way
  since a cell line can carry several distinct variants in the same gene
  (multiple rows per (cell line, gene)); collapsed into a binary
  presence/absence matrix instead (`aggfunc="max"` on a constant 1).
- **`combine_mutation_and_cnv()`** — column-wise union (`MUT_<gene>` /
  `CNV_<gene>` prefixes, outer join, missing filled with 0) into the single
  Mut_CNV matrix the proposal report's §4.1.2 preprocessing step expects.
- **`load_proteomics_wide()`** — proteomics needed no pivot at all.
  `Protein_matrix_averaged_20250211.tsv` (inside
  `data/raw/auxiliary/Proteomics_20250211.zip`) already ships wide (rows =
  Sanger `model_id`, columns = `uniprot_id`), just with an unusual 3-row
  header block (uniprot IDs / gene symbols / row-key labels) that needed a
  bespoke parse.

All three writers share `standard_model_id` (Sanger `SIDM*` IDs) as the join
key, matching the ID space `ingest_and_align.py` already standardizes on.

**Verified against real data:**
- CNV: full real pass (not a sample) → `(1507, 785)` wide matrix, plausible
  copy-number values (2–4 range, i.e. near-diploid to amplified).
- Proteomics: full real pass → `(948, 8453)` wide matrix, 38.6% missing
  (as expected for proteomics — this is exactly what the median-imputation
  step in `OmicsPreprocessingPipeline` exists to handle).
- Mutations: real 2M-row subset → binary `{0.0, 1.0}` presence matrix,
  confirmed no other values leak through.
- GE (RNA-seq, 78.8M rows): validated the streaming fix on a 500K-row real
  slice; the full real run was kicked off in the background against actual
  CCLE data as end-to-end confirmation before recommending this for Colab.

## Explicitly not done in this pass

- **Wiring `build_wide_matrices.py`'s output into `OmicsPreprocessingPipeline`
  end-to-end.** The wide matrices are now produced correctly, but
  `fit_transform()` hasn't yet been run against them (only against synthetic
  data so far) — that's the natural next check.
- **Downstream PyG integration.** The fused `(B, out_dim)` embedding is meant
  to bind into heterogeneous graph node features (per the task's stated
  purpose), but no PyTorch Geometric graph-construction code was touched —
  that's the GAT/graph module, out of scope here.
