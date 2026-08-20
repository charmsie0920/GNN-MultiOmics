# Graph Construction — `03_graph_construction.py` Real-Data Run

This documents the real run of the Phase 4 graph-construction pipeline
(`fetch_drug_smiles.py` -> `fetch_string_aliases.py` -> `03_graph_construction.py`),
what each output contains, and one real API bug the run surfaced.

## What the pipeline does

Builds a heterogeneous `cell_line` / `drug` / `protein` graph
(`torch_geometric.data.HeteroData`) from the PCA-compressed omics outputs
(`docs/preprocessing_and_fusion_module.md`), a local PubChem SMILES cache, and
STRING's PPI network — purely local-file -> `HeteroData`, per the plan
(`docs/plan/graph_construction_plan.md`). The two fetch scripts do the one-time
network acquisition/caching so the graph-construction script itself never
depends on a live API being reachable.

## 1. `fetch_drug_smiles.py`

Reads `data/raw/gdsc/screened_compounds_rel_8.5.csv` (621 rows / 542 unique
`DRUG_NAME`s), queries PubChem PUG REST once per **unique drug name** (not
per `DRUG_ID` — 71 names are re-screened under a second `DRUG_ID`, e.g.
`Erlotinib` -> IDs `1` and `1168`), and fans each resolved (SMILES, CID) pair
out to every `DRUG_ID` sharing that name. Output:
`data/raw/pubchem/gdsc_drug_smiles.csv`.

**Real bug hit and fixed:** PubChem's PUG REST API no longer returns a
`CanonicalSMILES` key even when that property is explicitly requested in the
URL — it silently substitutes `ConnectivitySMILES` instead, which the script
wasn't reading, so the first run resolved 0/542 names despite every HTTP
request succeeding with `200 OK`. Confirmed via direct `curl` against the
live endpoint. Fixed by requesting PubChem's current `SMILES` property name
instead (`.../property/SMILES/JSON`, reads `props["SMILES"]`), which returns
the same canonical structure under PubChem's current schema.

**Run:**
```
.venv/bin/python src/data/fetch_drug_smiles.py
```
- Unique drug names in GDSC: 542
- PubChem queries issued: 542
- **Names resolved: 420 / 542 (77.5%)**
- Names unresolved: 122 — almost entirely internal/proprietary compound
  codenames not registered in PubChem under that exact string (e.g. GSK
  internal codes like `GSK2256098C`, Bayer codes like `BAY-MPS-combo-1
  (paclitaxel 5 uM)`, kinase-panel codes like `JAK1_3715`, unlabeled
  numeric IDs like `123829`). None of the synonym fallbacks resolved these
  either, since GDSC's `SYNONYMS` column is empty for most of them.
- Output: `data/raw/pubchem/gdsc_drug_smiles.csv`, 621 rows, **47.2 KB**.

## 2. `fetch_string_aliases.py`

Downloads `9606.protein.aliases.v12.0.txt.gz` from
`stringdb-downloads.org` (same v12.0 release as the existing
`9606.protein.links.v12.0.txt.gz`).

**Run:**
```
.venv/bin/python src/data/fetch_string_aliases.py
```
- Source: `https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz`
- Output: `data/raw/string/9606.protein.aliases.v12.0.txt.gz`, **18.9 MB**.

Note: `data/raw/string/9606.protein.links.v12.0.txt.gz` (the PPI edge file
`ingest_and_align.py` already references) was also missing from this
checkout — `data/raw/` is gitignored and the file simply hadn't been
re-fetched onto this machine — so it was re-downloaded from the same
`stringdb-downloads.org` source (83 MB) to unblock this run. It's not a new
pipeline deliverable, just restoring an already-referenced input.

## 3. STRING PPI edges

- **Score threshold: `combined_score >= 700`** (STRING's standard
  "high confidence" cutoff on its 0-999 scale), set as a named constant
  (`PPI_SCORE_THRESHOLD` in `03_graph_construction.py`) for easy tuning.
- Streamed in 1M-row chunks (`pandas.read_csv(..., chunksize=1_000_000)`) —
  13,715,404 total rows never materialized at once.
- **Verified STRING already lists both directions of every interaction**
  (spot-checked: row 2's `A B score` has a matching `B A score` row
  elsewhere in the file, and the total row count matches STRING's documented
  total exactly) — so no manual edge mirroring was needed, unlike the plan's
  contingency for the case where it wasn't already symmetric.
- Result: **473,860 / 13,715,404 rows kept** at the score≥700 filter, spanning
  **16,203 distinct proteins**.

## 4. Drug-target mapping

- GDSC's `TARGET` column (comma-separated free text) yields **380 distinct
  target symbol tokens** across all drugs.
- **266 / 380 (70%) resolved** to an `ENSP` via the STRING alias file
  (loaded once, filtered to only these 380 referenced symbols — not the full
  ~40M-row alias file).
- **114 unresolved**, almost entirely non-gene-symbol free text rather than
  ambiguous gene families: mechanism/class descriptions (`"Alkylating
  agent"`, `"DNA crosslinker"`, `"Proteasome"`), receptor-family shorthand
  that doesn't exact-match a single gene (`"FGFR"`, `"BCL-XL"`, `"CDK"`,
  `"IKK"`), and a few genuinely malformed tokens from multi-target cells
  that don't split cleanly on comma (`"TTK and microtubules"`,
  `"Tankyrase 1/2 (PARP5a"`). No fuzzy/family matching was attempted, per
  the plan's explicit "don't over-engineer" scope — these are logged, not
  guessed.
- Result: **683 `(drug, targets, protein)` edges.**

## 5. Final `HeteroData` object

**Run:**
```
.venv/bin/python src/data/03_graph_construction.py
```

```
node 'cell_line': num_nodes=532, feature_dim=384
node 'drug': num_nodes=498, feature_dim=2048
node 'protein': num_nodes=16203, feature_dim=128
edge ('protein', 'interacts_with', 'protein'): num_edges=473860
edge ('drug', 'targets', 'protein'): num_edges=683
total tensor memory: 19.82 MB
```

- `cell_line`: all 532 PCA rows matched across transcriptomics/genomics/
  proteomics (full overlap, none dropped) — inner join on
  `standard_model_id`. (Each PCA CSV has 532 data rows, not 533 — the
  `wc -l` count used during planning included the header line.)
- `drug`: 498 / 621 SMILES rows parsed successfully by RDKit (123 dropped:
  122 unresolved names + rows with malformed/empty SMILES); 2048-bit Morgan
  fingerprints (`radius=2`), `float32`.
- `protein`: zero-initialized `(16203, 128)` placeholder embeddings, to be
  learned by an `nn.Embedding` in the GNN training stage.
- Integrity check passed: every edge type's indices are within range for
  their endpoint node type.
- Saved to `data/processed/hetero_graph.pt`. Round-tripped with
  `torch.load(..., weights_only=False)` — reloaded shapes/edge counts match
  exactly.
- **Peak RSS during the run: ~876 MB** (measured via `/proc/<pid>/status`
  `VmHWM`, since `/usr/bin/time -v` isn't available in this environment) —
  comfortably inside the 4-8GB budget, driven mainly by the STRING alias
  file streaming pass and PPI chunked read, not the final tensors (19.82 MB).

## What's next

`data/processed/hetero_graph.pt` is ready to be paired with GDSC's IC50
labels (`data/processed/aligned/gdsc2_response_master.csv`, keyed by
`standard_model_id` + `drug_id`, both of which match `cell_line.node_ids`
and `drug.node_ids` here) for the GNN training stage.
