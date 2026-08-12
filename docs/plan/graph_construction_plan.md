# Plan: `03_graph_construction.py` — HeteroData Graph Construction

## Context

The pipeline so far has three completed stages, all documented in `docs/`:
`ingest_and_align.py` (raw → ID-crosswalked long CSVs) → `build_wide_matrices.py`
(long → wide per-modality matrices) → `omics_preprocessing.py` /
`00_run_preprocessing.py` (wide → PCA-compressed CSVs). The user has now
uploaded the real PCA outputs (`data/processed/{transcriptomics,genomics,proteomics}_pca.csv`,
533 rows × 128 dims each, keyed on Sanger `standard_model_id`), unblocking
the next stage: **Phase 4, graph construction**, explicitly called out as
"unblocked but not yet started" in `docs/phase2_wide_matrices_and_split_plan.md`.

Investigation surfaced three real gaps between the user's original task
description and repo state, resolved with the user directly:

1. **PCA files** — now present (confirmed above). No special-casing needed;
   `03_graph_construction.py` reads them directly and fails loudly if missing.
2. **Drug SMILES** — no SMILES file exists anywhere in the repo. User wants
   it fetched live from the PubChem API, keyed off GDSC's own drug list
   (`data/raw/gdsc/screened_compounds_rel_8.5.csv`, 621 rows / 542 unique
   drug names), cached to disk, with a documented report of source/size/coverage
   so the user can decide whether to keep the raw cache.
3. **Protein ID space mismatch** — STRING PPI edges use Ensembl protein IDs
   (`ENSP*`); GDSC's `TARGET` column is free-text gene symbols (e.g. `EGFR`,
   or multi-target families like `PDGFR, KIT, VEGFR, FLT3, RET, CSF1R`).
   No alias file exists locally. User approved fetching
   `9606.protein.aliases.v12.0.txt.gz` from STRING (matching the existing
   `9606.protein.links.v12.0.txt.gz` version) to bridge gene symbol → `ENSP`.

This plan keeps the acquisition/caching work in small prerequisite scripts
(matching the existing `ingest_and_align.py` / `build_wide_matrices.py`
separation-of-concerns pattern) and keeps `03_graph_construction.py` itself
purely a **local-file → `HeteroData`** builder, per the user's explicit
request to scope this script to graph construction only.

## Deliverables

0. `docs/plan/graph_construction_plan.md` — this plan, committed to the repo as the first step of implementation (per user request, so it's reviewable/version-controlled alongside the other `docs/*.md` pipeline-stage writeups).
1. `src/data/fetch_drug_smiles.py` — one-time/idempotent PubChem SMILES fetch.
2. `src/data/fetch_string_aliases.py` — one-time STRING alias file download.
3. `src/data/03_graph_construction.py` — the main deliverable: builds and saves the `HeteroData` object.
4. `docs/graph_construction.md` — detailed documentation of sources, decisions, and stats (matching the existing `docs/*.md` style).
5. `requirements.txt` — add `rdkit` and `torch_geometric` (no `requests` needed; stdlib `urllib.request` matches `ingest_and_align.py`'s existing no-extra-deps style for the two fetch scripts).

## 1. `src/data/fetch_drug_smiles.py`

- Reads `data/raw/gdsc/screened_compounds_rel_8.5.csv` (`DRUG_ID, DRUG_NAME, TARGET`), dedupes on `DRUG_ID`.
- For each drug, queries PubChem PUG REST
  (`.../compound/name/{drug_name}/property/CanonicalSMILES/JSON`), with:
  - Rate limiting (PubChem's documented ~5 req/s ceiling).
  - A synonym fallback: if `DRUG_NAME` fails, retry with each entry in `SYNONYMS`.
  - Resumability: skip drugs already present in the output CSV so re-running after a partial failure doesn't re-hit the API.
- Output: `data/raw/pubchem/gdsc_drug_smiles.csv` (`drug_id, drug_name, canonical_smiles, pubchem_cid`).
- Final `print()`: rows attempted, rows resolved, rows unresolved (with drug names logged), and the output file's size on disk — directly answers the user's "let me know how big the raw data would be" ask.

## 2. `src/data/fetch_string_aliases.py`

- Downloads `9606.protein.aliases.v12.0.txt.gz` from STRING's public download endpoint (same v12.0 as the existing `data/raw/string/9606.protein.links.v12.0.txt.gz`) to `data/raw/string/`.
- Skips download if the file already exists (idempotent).
- Final `print()`: downloaded file size, so the user can decide whether to keep it (same reporting pattern as the SMILES fetch).

## 3. `src/data/03_graph_construction.py`

### Node construction

- **`cell_line`**: inner-join `transcriptomics_pca.csv` / `genomics_pca.csv` / `proteomics_pca.csv` on `standard_model_id` (defensively — log if any of the 533 rows fail to match across all three, rather than assuming). Concatenate the three 128-dim blocks → `x` of shape `(num_cell_lines, 384)`. Store the ordered `standard_model_id` list as `data['cell_line'].node_ids` for future traceability (e.g. wiring GDSC labels back on in a later stage) — cheap, standard practice, not scope creep.
- **`drug`**: load `data/raw/pubchem/gdsc_drug_smiles.csv`, parse each `canonical_smiles` with RDKit (`Chem.MolFromSmiles`), drop unparseable rows (log count), compute a 2048-bit Morgan fingerprint per molecule (`rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)`, the current non-deprecated RDKit API) → `x` of shape `(num_drugs, 2048)`, `float32`. Store `data['drug'].node_ids` (GDSC `drug_id`) similarly.
- **`protein`**: node universe = union of (a) every `ENSP` endpoint that survives the STRING PPI score filter below, and (b) every `ENSP` a resolved drug target maps to (some targets may have no high-confidence PPI edges but still need a node for the `targets` edge). Features: zero-initialized `(num_proteins, 128)` embedding — explicitly a placeholder per the task ("basic embeddings... secondary interactors"), matching the cell-line PCA dim loosely for consistency, to be learned by an `nn.Embedding` in the GNN stage (not built here).

### Edge construction

- **`(protein, interacts_with, protein)`**: stream `9606.protein.links.v12.0.txt.gz` in chunks (`pandas.read_csv(..., chunksize=1_000_000)`, never materializing all 13.7M rows at once — directly satisfies the memory constraint). Filter `combined_score >= 700` (STRING's standard "high confidence" cutoff; called out as a named constant so it's easy to tune later). Build the `ENSP → node index` mapping incrementally as new proteins are encountered, appending to plain Python int lists (cheap — no large string-tuple lists retained), converting to a `(2, num_edges)` `int64` tensor only at the end. Verify during implementation whether STRING already lists both directions (A→B and B→A) for each interaction; if not, mirror edges so the symmetric relation is represented in both directions.
- **`(drug, targets, protein)`**: for each drug's `TARGET` column (comma-separated gene symbols; skip empty/`"not defined"`), map each symbol to `ENSP` via a `gene_symbol → ENSP` dict built from `9606.protein.aliases.v12.0.txt.gz` (loaded once, filtered to only the symbols actually referenced — memory-bounded by the small target vocabulary, not the full alias file). Family-level names that don't exact-match any alias (e.g. `PDGFR`, `VEGFR` as opposed to `PDGFRA`/`FLT1`) are logged as unresolved rather than guessed — no fuzzy/family expansion, matching "don't over-engineer." Add any newly-referenced target `ENSP` as a protein node if not already present from the PPI step.

### Memory-efficiency techniques used (to call out explicitly in the script's docstring/comments)

- Chunked STRING PPI parsing with a confidence-score filter (no full 13.7M-row DataFrame in memory).
- Incremental dict-based node indexing instead of building and deduplicating a full edge string list.
- Alias file filtered to only the referenced gene-symbol vocabulary, not loaded wholesale.
- All ID→index bookkeeping done in plain Python dicts/lists (cheap); only final tensors go to `torch`.

### Output

- `print()` of full graph metadata: per node type `num_nodes` and feature dim, per edge type `num_edges`; then total tensor memory (sum of `tensor.numel() * tensor.element_size()` across every node/edge store) reported in MB.
- `torch.save(data, "data/processed/hetero_graph.pt")`.
- Light integrity check before saving: assert every edge type's `edge_index.max() < num_nodes` for its endpoint type (catches an indexing bug immediately rather than silently producing a broken graph).

## 4. `docs/graph_construction.md`

Following the existing `docs/*.md` style (source, what was run, output shapes, decisions, gotchas):
- PubChem fetch: endpoint used, drug list source, resolved/unresolved counts and which drugs failed, output file size.
- STRING aliases fetch: source URL/version, output file size.
- STRING PPI: score threshold chosen and why, resulting node/edge counts after filtering.
- Drug-target mapping: how many distinct target symbols existed, how many resolved, notable unresolved cases (e.g. receptor family names).
- Final `HeteroData` object: full metadata printout, memory footprint, output path.

## Verification

1. Run `fetch_drug_smiles.py`, confirm `data/raw/pubchem/gdsc_drug_smiles.csv` is created with a plausible resolution rate, check the printed size/coverage report.
2. Run `fetch_string_aliases.py`, confirm the alias file downloads and its printed size.
3. Run `03_graph_construction.py` end-to-end; confirm the printed metadata (node/edge counts per type, memory size) looks sane (e.g. `cell_line` ≈ 532, `drug` ≈ resolved-SMILES count, `protein` in the tens of thousands, `interacts_with` edge count consistent with a score≥700 STRING filter).
4. `torch.load("data/processed/hetero_graph.pt", weights_only=False)` and re-print metadata to confirm the round-trip matches what was saved.
5. Confirm total process RSS stayed within the 4–8GB budget during the STRING streaming step (e.g. spot-check with `/usr/bin/time -v` or a quick `resource.getrusage` print).
