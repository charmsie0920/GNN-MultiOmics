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

## Explicitly not done in this pass

- **Real data wiring.** `omics_preprocessing.py` has not been connected to
  `data/processed/aligned/*`. Two things are needed first:
  1. A long → wide pivot of `ingest_and_align.py`'s output (currently one row
     per gene per model; `OmicsPreprocessingPipeline` needs one row per model).
  2. Proteomics ingestion — `ingest_and_align.py` doesn't parse
     `data/raw/auxiliary/Proteomics_20250211.zip` yet.
- **Downstream PyG integration.** The fused `(B, out_dim)` embedding is meant
  to bind into heterogeneous graph node features (per the task's stated
  purpose), but no PyTorch Geometric graph-construction code was touched —
  that's the GAT/graph module, out of scope here.
