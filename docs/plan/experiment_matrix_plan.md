# Plan: Full Experiment Matrix — Omics × Drug Representation × Architecture

## Status (update as phases complete)

- [x] **Phase 0 — Shared infrastructure** (done: `experiment_utils.py` built,
      fusion generalized to modality subsets, `FingerprintEncoder` added,
      `random_state=42` pinned in both baselines, `scipy`/`xgboost` added to
      requirements. **Verified**: proteomics-only RF re-run through the new
      shared path reproduces the documented test RMSE 2.0000 / PCC 0.7460 /
      AUC 0.7520 / F1 0.6435 and floor 2.7097 exactly.)
- [x] **Phase 1 — Flat-feature matrix (RF + MLP), 42 runs** (21 + 21; grew from
      28 when the population-control arm was added). RF 482.8s, MLP 903.7s.
      Docs: `rf_ablation_results.md`, `mlp_ablation_results.md`.
- [x] **Phase 2 — Cross-attention matrix, 12 runs** (grew from 8). 1466.0s.
      Doc: `cross_attention_ablation_results.md`.
- [x] **Phase 3 — Graph linkage + GNN matrix, 4 runs.** 669.3s. Linkage added
      4,341 driver-mutation edges linking 526/532 cell lines.
      Doc: `gnn_ablation_results.md`.
- [x] **Phase 4 — XGBoost ensemble refinement.** 225.7s. **Negative result** —
      refinement made both best models slightly worse (−0.9%, −0.7%), opposite
      to the paper's +19.7%. Doc: `ensemble_refinement_results.md`.
- [x] **Phase 5 — Consolidation.** `docs/results.md` auto-generated from the
      result CSVs by `experiments/build_results_table.py` (58 runs).

**Mid-execution amendment (approved):** the fingerprint arm covers only 111,799
pairs vs one-hot's 134,764 (123 GDSC drugs never resolved to a SMILES), so
comparing them directly confounded representation with population. A third
`onehot_restricted` arm was added — one-hot features on the fingerprint-covered
rows — making the drug-representation comparison valid. Total runs 42 → 58.

## Context

`MoGraphDRP` (the paper we're benchmarked against) hits RMSE 0.66 / PCC 0.97
using a multi-branch architecture (4 omics types, dual drug encoder,
bilinear attention, XGBoost refinement) validated by extensive ablations
(GNN-architecture choice, single-fingerprint-only, single-omics-removed,
with/without-graph, with/without-ensemble). Our best result so far
(cross-attention fusion + MLP, tri-omics, one-hot drugs) is RMSE 1.27 —
good, but every experiment run to date shares two placeholders the paper
already moved past: one-hot drug identity (never Morgan fingerprints, despite
the graph pipeline already producing them) and no graph/relational structure
at all. The goal of this plan is to close that gap with a systematic ablation
matrix, mirroring the paper's ablation structure, so we know exactly which
piece (omics breadth, drug representation, fusion mechanism, graph structure)
is driving improvement — not just that "the fancier model is better."

Three research passes (Explore agents) confirmed the current state:
**already run** (all one-hot drug rep, all `GroupShuffleSplit`-by-cell-line
70/15/15): fused-tri-omics RF (RMSE 2.19), proteomics-only RF (RMSE 2.00,
unexplained beat of the fused model), cross-attention+MLP tri-omics (RMSE
1.27, best so far). **Exist as code but never logged**: MLP-on-fused,
genomics-only RF, transcriptomics-only RF, late-fusion RF ensemble.
**Real gaps**: zero Morgan-fingerprint experiments anywhere;
`MultiOmicsCrossAttentionFusion` is hardcoded to exactly 3 modalities (can't
run on a 2-omics subset without generalizing it); no GCN/GAT code exists in
the repo at all; `hetero_graph.pt` has no edges linking `cell_line` nodes to
`drug`/`protein` nodes (a disconnected component — a GNN can't propagate any
signal into a cell-line prediction as built today); no `docs/results.md`
running comparison table exists (proposed twice in existing docs, never
created).

**Decisions confirmed with the user before this plan was written:**
1. GNN scope: heterograph only (cell_line/drug/protein + STRING PPI, drug
   nodes = Morgan fingerprint vectors) — **not** a separate per-drug
   atom-level molecular graph like the paper's Fig 2.
2. Cell-line linkage: add `(cell_line, has_mutation, protein)` structural
   edges from GDSC mutation data (independent of IC50 — no label leakage),
   so PPI/drug-target signal can reach cell-line embeddings via message
   passing.
3. Matrix scope: the full literal matrix — all 7 omics subsets × both drug
   representations × every applicable model.
4. Docs: consolidated per model-family (one doc per model type covering all
   its matrix cells), plus one master `docs/results.md`.
5. GNN omics sweep: tri-omics only (one graph) — GNN's own ablation axis is
   architecture (GCN vs GAT) × with/without the new cell_line→protein edges,
   not omics subset (that's covered by RF/MLP/cross-attention already).
6. XGBoost ensemble refinement: included, applied post-hoc to the
   best-performing model(s) only, not as a full matrix axis.

## The experiment matrix (~42 runs total)

**Population held fixed across every cell** (for apples-to-apples
comparability): the same 532 cell lines with full tri-omics PCA coverage,
same 134,764 (cell line, drug) pairs, same GroupShuffleSplit-by-cell-line
70/15/15 with `random_state=42` pinned everywhere (fixes the two currently
unpinned scripts, `rf_baseline.py` and `mlp_baseline.py`) — every cell sits
on an identical split so the comparison table is valid. Mean-only floor RMSE
is computed once (same population/split → same floor for every cell), not
per-cell.

**Omics subsets (7):** GE, Mut_CNV, Proteomics, GE+Mut_CNV, GE+Proteomics,
Mut_CNV+Proteomics, GE+Mut_CNV+Proteomics (tri-omics).

**Drug representations (2):** one-hot (295-dim identity) vs. Morgan
fingerprint (2048-dim, from `data/raw/pubchem/gdsc_drug_smiles.csv` via
RDKit — reuses the exact parsing already in `03_graph_construction.py`'s
`build_drug_nodes()`).

| Model | Omics × drug-rep coverage | # runs |
|---|---|---|
| RF (flat concat) | all 7 omics subsets × both drug reps | 14 |
| MLP (flat concat) | all 7 omics subsets × both drug reps | 14 |
| Cross-attention fusion + MLP head | omics subsets with ≥2 modalities (4 of 7) × both drug reps | 8 |
| GNN (GCN, GAT) over heterograph | tri-omics only × {with, without} new cell_line→protein edges | 4 |
| **Subtotal** | | **40** |
| XGBoost ensemble refinement | with/without, on best 1-2 models from above (decided after Phases 1-3 land) | ~2-4 |
| **Total** | | **~42-44** |

Cross-attention needs ≥2 modalities to attend across (a single modality has
nothing to cross-attend with — falls back to plain MLP, already covered by
the MLP row). GNN drug representation is inherently the fingerprint (that's
what the graph's drug nodes already are); GNN doesn't get a separate one-hot
run.

## Metrics (uniform across all ~42 runs)

Regression: **RMSE, MAE, R², PCC (Pearson), SCC (Spearman)** — matches the
paper's own metric set (§3.2) exactly, superset of what any current script
computes. Classification (binarized via median `ln_ic50` split, matching
existing convention): **AUC, F1**. Plus reference numbers already used in
this repo's docs: mean-only floor RMSE, trainable param count, fit
time/peak RSS.

## Phase 0 — Shared infrastructure (blocks every other phase)

1. **Pin `random_state=42`** in `src/data_engineering/models/rf_baseline.py`
   and `experiments/Early Fusion & MLP/mlp_baseline.py` (currently `None` —
   the reproducibility gap flagged by the earlier RF results doc).
2. **New `src/data/experiment_utils.py`** — extract the boilerplate currently
   duplicated near-verbatim across `rf_baseline.py`, `mlp_baseline.py`,
   `experiments/Cross Attention Fusion/cross_attention_baseline.py`:
   - `load_omics_subset(modalities: list[str]) -> dict[str, pd.DataFrame]` —
     reads the requested subset of `{transcriptomics,genomics,proteomics}_pca.csv`.
   - `load_drug_features(mode: Literal["onehot","fingerprint"]) -> dict[str, np.ndarray]`
     — one-hot branch = existing `pd.factorize` logic; fingerprint branch =
     Morgan-fingerprint builder factored out of `03_graph_construction.py`'s
     `build_drug_nodes()` so it's not copy-pasted a second time.
   - `build_pair_matrix(...)`, `grouped_split(...)` (reuse existing logic,
     just centralized).
   - `evaluate(y_true, y_pred) -> dict` returning all 7 metrics above.
   - `peak_rss_gb()`, `print_report(...)` (standard formatted console output
     matching the existing docs' style, so writing the results docs is a
     copy-paste of the printed block like today).
3. **Generalize `MultiOmicsCrossAttentionFusion`**
   (`src/models/cross_attention_fusion.py`) to accept a
   `modalities: list[str]` constructor arg instead of the hardcoded
   `MODALITIES` triple — build `self.pairs` (all ordered pairs within the
   given list) and `output_proj`'s input width (`d_model * len(pairs)`)
   dynamically. Existing tri-omics callers keep working unchanged (default
   `modalities=list(MODALITIES)`).
4. **Small fingerprint encoder** for MLP/cross-attention experiments using
   the Morgan-fingerprint drug rep: `Linear(2048→128) → ReLU → Dropout →
   Linear(128→128)` before concatenation with the fused omics embedding
   (raw 2048 sparse bits directly concatenated with a 256-dim continuous
   embedding would be badly unbalanced). RF doesn't need this — trees handle
   the raw 2048-dim input directly.

## Phase 1 — Flat-feature matrix (RF + MLP), 28 runs

New `experiments/Full Matrix/rf_matrix.py` and
`experiments/Full Matrix/mlp_matrix.py`, each looping the 7 omics subsets ×
2 drug reps via Phase 0's shared utils. Build each of the 14 distinct
(omics-subset, drug-rep) pair-matrices **once**, reuse for both RF and MLP
(halves data-prep work). Log every run's metrics into a running
in-script dict/CSV for Phase 5.

**Doc:** `docs/rf_ablation_results.md` and `docs/mlp_ablation_results.md`
— one consolidated doc per model family, each with one results table (14
rows) covering all its omics/drug-rep cells, following the existing
results-doc style (what ran, why, numbers, interpretation).

## Phase 2 — Cross-attention matrix, 8 runs

New `experiments/Full Matrix/cross_attention_matrix.py`, using the
generalized fusion class from Phase 0, looping the 4 omics subsets with
≥2 modalities × 2 drug reps.

**Doc:** `docs/cross_attention_ablation_results.md` (8-row table).

## Phase 3 — Graph linkage + GNN matrix, 4 runs

1. **New `src/data/04_link_cell_lines.py`** — loads `hetero_graph.pt`, adds
   `(cell_line, has_mutation, protein)` edges: for each cell line, its
   mutated genes (from `data/processed/aligned/mutations_aligned.csv`
   presence, already used for the Mut_CNV wide matrix) mapped to `ENSP` via
   the same gene-symbol→ENSP alias dict already built in
   `03_graph_construction.py`'s `build_gene_symbol_to_ensp()` (factor that
   function out to be importable, not copy-pasted). Re-saves
   `data/processed/hetero_graph.pt` with the new edge type. Matches the
   repo's existing small-prerequisite-script pattern.
2. **New `src/models/hetero_gnn.py`** — a GCN variant and a GAT variant,
   both `torch_geometric.nn.HeteroConv` wrapping per-edge-type conv layers
   (`SAGEConv`/`GCNConv` for the GCN variant, `GATConv` for the GAT variant)
   over `interacts_with`, `targets`, and the new `has_mutation` edges. A
   prediction head takes `(cell_line_embedding, drug_embedding)` →
   MLP → `ln_ic50`. Protein node features start as the existing
   zero-initialized placeholder, now genuinely learned via `nn.Embedding`
   during training (the graph doc's already-flagged next step).
3. **New `experiments/GNN Ablation/gnn_baseline.py`** — joins GDSC IC50
   labels onto `hetero_graph.pt`'s `cell_line.node_ids` /
   `drug.node_ids`, trains {GCN, GAT} × {with, without the new
   `has_mutation` edges} = 4 runs, same split/metrics as every other phase.

**Doc:** `docs/gnn_ablation_results.md` (4-row table), explicitly including
whether the new cell_line→protein linkage helped (mirrors the paper's own
"without graph" ablation finding that graph structure mattered most).

## Phase 4 — XGBoost ensemble refinement (post-hoc, after Phases 1-3 results are in)

**New `src/models/ensemble_refinement.py`** — `XGBRegressor` per the paper's
§2.4 hyperparams (`n_estimators=100, max_depth=6, learning_rate=0.05,
subsample=0.8`), trained on `[fusion_vector ⊕ initial_prediction] →
residual-corrected ln_ic50`, applied to whichever 1-2 models score best
across Phases 1-3 (data-driven choice — can't be pre-specified before those
results land). Adds `xgboost` to `requirements.txt`.

**Doc:** `docs/ensemble_refinement_results.md` — with-vs-without comparison
on the chosen model(s), matching the paper's Table 7 style.

## Phase 5 — Consolidation

**`docs/results.md`** — the master comparison table across all ~42-44 runs
(finally building what's been proposed twice and never created). One row
per run: Experiment ID, Omics Subset, Drug Rep, Model, Architecture Details,
RMSE, MAE, R², PCC, SCC, AUC, F1, Params, Fit Time, Peak RSS, link to its
model-family doc. This single well-formatted markdown table is also the
literal artifact to copy/paste into Notion for the requested comprehension
table — Notion imports markdown tables directly, so no separate export
format is needed.

## Files to create

- `src/data/experiment_utils.py` (shared utilities)
- `src/data/04_link_cell_lines.py` (graph edge addition)
- `src/models/hetero_gnn.py` (GCN + GAT model classes)
- `src/models/ensemble_refinement.py` (XGBoost wrapper)
- `experiments/Full Matrix/rf_matrix.py`
- `experiments/Full Matrix/mlp_matrix.py`
- `experiments/Full Matrix/cross_attention_matrix.py`
- `experiments/GNN Ablation/gnn_baseline.py`
- `docs/rf_ablation_results.md`, `docs/mlp_ablation_results.md`,
  `docs/cross_attention_ablation_results.md`, `docs/gnn_ablation_results.md`,
  `docs/ensemble_refinement_results.md`, `docs/results.md`

## Files to modify

- `src/models/cross_attention_fusion.py` (generalize to arbitrary modality
  subsets)
- `src/data_engineering/models/rf_baseline.py`,
  `experiments/Early Fusion & MLP/mlp_baseline.py` (pin `random_state=42`)
- `src/data/03_graph_construction.py` (factor `build_gene_symbol_to_ensp()`
  and the Morgan-fingerprint builder into importable functions, no behavior
  change)
- `requirements.txt` (add `xgboost`)

## Verification

1. Phase 0: unit-check `experiment_utils.evaluate()` against one already-known
   result (re-run proteomics-only RF through the new shared path, confirm it
   reproduces RMSE 2.00/PCC 0.746 from the existing doc, within floating-point
   tolerance) — proves the refactor didn't change behavior.
2. Each phase's matrix script: confirm it produces exactly the expected
   number of rows (14/14/8/4) with no `NaN` metrics, and that every run's
   split shares the same cell-line group assignment (assert train/val/test
   cell-line sets are identical across runs, not just same *sizes*).
3. Phase 3: confirm the new `has_mutation` edges pass the same
   `edge_index.max() < num_nodes` integrity check already used in
   `03_graph_construction.py`; confirm cell_line nodes are no longer a
   disconnected component (a simple reachability check from any cell_line
   node to at least one protein node).
4. Phase 5: confirm `docs/results.md`'s table row count matches the total
   run count from Phases 1-4, and spot-check 2-3 rows against their
   individual model-family doc numbers for consistency.
