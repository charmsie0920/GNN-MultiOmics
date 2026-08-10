# Data Ingestion & Alignment — `ingest_and_align.py` Real-Data Run

This documents the real run of `src/data/ingest_and_align.py` against the
project's full CCLE/DepMap + GDSC + STRING data (`data/raw/`, 7.3GB), what
each output file contains, and one real bug the run surfaced.

## What the script does

Turns raw, heterogeneous CCLE/GDSC/STRING files into ID-crosswalked,
long-format CSVs (`data/processed/aligned/`) — one row per (cell line, gene)
for the omics tables, one row per (cell line, drug) for GDSC responses. It
does **not** normalize, pivot, or reduce dimensionality; that's
`build_wide_matrices.py` + `OmicsPreprocessingPipeline`'s job (see
`docs/preprocessing_and_fusion_module.md`).

**Why a crosswalk is needed at all:** the raw sources don't share one ID
system. CCLE/DepMap's own files use Sanger Cell Model Passport IDs
(`SIDM*`, e.g. `SIDM00776`) in their `model_id` column. GDSC's dose-response
file separately references `SANGER_MODEL_ID`. `model_list_20260724.csv` also
carries a `BROAD_ID` (DepMap-style `ACH-*`, e.g. `ACH-000205`), `CCLE_ID`,
and `COSMIC_ID` for the same cell line. `load_model_crosswalk()` builds a
lookup so every downstream row can be joined on a shared identity regardless
of which ID system the original file used.

## The run

```
.venv/bin/python src/data/ingest_and_align.py
```

Exit code 0, ~23 minutes locally (no GPU used — this is pure CPU/text
processing; see the "Colab" note at the bottom). `missing_gdsc_model_ids` in
the summary came back **empty** — every GDSC response row resolved through
the crosswalk successfully.

## Output files

All under `data/processed/aligned/` (gitignored — regenerate by re-running
the script, don't commit these):

| File | Rows | Size | Grain |
|---|---|---|---|
| `model_crosswalk.csv` | 2,266 | 224K | one row per cell line |
| `gdsc2_response_master.csv` | 242,036 | 64M | one row per (cell line, drug) |
| `rnaseq_aligned.csv` | 78,839,610 | 8.5G | one row per (cell line, gene) |
| `mutations_aligned.csv` | 14,370,318 | 2.3G | one row per (cell line, variant) |
| `cnv_aligned.csv` | 1,230,390 | 110M | one row per (cell line, gene) |
| `omics_availability_by_model.csv` | 2,266 | 80K | one row per cell line |
| `string_summary.json` | — | — | existence check only, not parsed |

### `model_crosswalk.csv`

The unified ID table. Columns: `sanger_model_id, depmap_id, cell_line_name,
broad_id, ccle_id, cosmic_id, rrid, tissue, cancer_type`.

```
sanger_model_id,depmap_id,cell_line_name,...,tissue,cancer_type
SIDM01774,ACH-000205,PK-59,...,Pancreas,Pancreatic Carcinoma
```

### `gdsc2_response_master.csv`

GDSC2 drug-response rows (the IC50 labels), enriched with the crosswalk and
tacked-on omics availability flags. Key columns: `standard_model_id,
drug_id, drug_name, putative_target, pathway_name, ln_ic50, auc`.
`ln_ic50` is the regression target.

```
gdsc2_response,ACH-001711,...,SIDM01132,PFSK-1,...,1003,Camptothecin,TOP1,...,ln_ic50=-1.463887,auc=0.93022,...
```

### `rnaseq_aligned.csv`, `cnv_aligned.csv`, `mutations_aligned.csv`

The three omics tables `build_wide_matrices.py` pivots. All keyed on
`standard_model_id`, which for these three files is the raw CCLE/DepMap
`model_id` value as-is (Sanger `SIDM*` format — see the bug note below for
why this matters).

- **RNA-seq**: value column used downstream is `rsem_tpm`. Confirmed real,
  raw-scale (not pre-log-transformed) — e.g. `rsem_tpm=3.5521` — so the
  `log1p` step in `OmicsPreprocessingPipeline` is correct and won't
  double-transform.
- **CNV**: value column is `total_copy_number` (e.g. `2.0` = neutral/diploid).
- **Mutations**: one row per called variant per gene per cell line (not one
  row per gene) — a gene with 3 distinct variants in one cell line produces
  3 rows. This is why `build_wide_matrices.py` collapses mutations into a
  binary presence/absence matrix rather than pivoting a single value column.

### `omics_availability_by_model.csv` / the `*_rows` columns — known bug, not yet fixed

**These numbers are wrong for ~86% of cell lines and shouldn't be trusted
yet.** Verified: only 313 of 2,266 cell lines show any nonzero
`rnaseq_rows`/`cnv_rows`/`mutation_rows`, even though the omics data is
actually present for far more of them.

**Root cause:** `summarize_availability()` builds its lookup key from
`model.depmap_id or model.sanger_model_id` — i.e. it prefers the `ACH-*`
DepMap ID when available (which is nearly always). But `stream_omics_table()`
populates `standard_model_id` directly from the raw file's `model_id` column,
which is Sanger `SIDM*` format, not `ACH-*`. Confirmed directly: checked a
500K-row RNA-seq sample against all 2,266 crosswalk records —
**`sanger_model_id` matched 1,203 of them; `depmap_id` matched 0.**
`build_response_master()` has the same bug (joins on `depmap_id`), which is
why `gdsc2_response_master.csv`'s `rnaseq_rows`/`cnv_rows`/`mutation_rows`
columns are also unreliable (e.g. the sample row above shows all zeros for a
cell line that does have omics data).

**This does not affect the actual pipeline.** `build_wide_matrices.py`
pivots `rnaseq_aligned.csv`/`cnv_aligned.csv`/`mutations_aligned.csv`
directly on their own internally-consistent `standard_model_id` column
(Sanger `SIDM*` throughout) — it never goes through `summarize_availability()`
or `build_response_master()`. Only the diagnostic availability CSV and the
three `*_rows` columns bolted onto `gdsc2_response_master.csv` are wrong;
if those are ever used to decide which cell lines have full omics coverage,
that decision would currently be wrong for most of the dataset.

**Fix, not yet applied:** key `summarize_availability()`/
`build_response_master()` on `sanger_model_id` instead of (or before)
`depmap_id`, since that's the ID space the omics tables actually populate.

## What's next

`build_wide_matrices.py` hasn't been run against this real output yet (only
validated so far against real CNV, real proteomics, and a 2M-row real
mutation subset — see `docs/preprocessing_and_fusion_module.md`). That's the
next step, now that real `rnaseq_aligned.csv`/`cnv_aligned.csv`/
`mutations_aligned.csv` exist to pivot.

## Note on Colab

This run used no GPU — it's pure CPU/text processing, ~23 minutes locally.
Running it on Colab instead buys nothing for this stage; Colab is worth
reserving for the actual GNN training later. If you do still want it on
Colab, see the file list and mount/copy instructions from the prior runbook
discussion in this conversation.
