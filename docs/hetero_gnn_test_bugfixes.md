# `src/models/test/hetero_gnn.py` — Bug Fixes and Test Runs

**Files:** [`src/models/test/hetero_gnn.py`](../src/models/test/hetero_gnn.py),
[`src/models/test/dataset.py`](../src/models/test/dataset.py),
[`src/models/test/05_train.py`](../src/models/test/05_train.py)
**Reference model:** [`src/models/hetero_gnn.py`](../src/models/hetero_gnn.py) (the tracked `HeteroGNN`, backing E12/E24/E56/E58 in [`results.md`](./results.md))
**Date:** 2026-08-28

## Background

`src/models/test/` is a standalone, uncommitted-to-the-matrix prototype of a
heterogeneous GNN (`HeteroIC50GNN`) added in commit `a943717` ("Heterogenous
GNN model uploaded, run 05_train.py"). It duplicates the same
`cell_line`/`drug`/`protein` graph (`src/graph/hetero_graph.pt`) the tracked
`HeteroGNN` uses, but with its own simplified architecture, its own dataset
prep script, and its own training loop — none of which had been reviewed or
had results logged anywhere. Auditing it surfaced four problems, three in
code and one in the training methodology. All four are fixed below and the
model was re-run after each fix to quantify the effect.

## Bugs found and fixed

### 1. `cell_line` nodes never received a graph message

**File:** `hetero_gnn.py`

Both `HeteroConv` layers only wired up `protein↔protein`,
`drug→protein`, and `protein→drug`. No edge type touched `cell_line` at all,
so `h_dict["cell_line"]` was set once by `cell_proj` and passed through
untouched by both conv layers — the "GNN" only ever refined the drug and
protein embeddings. The cell-line side was, in effect, a plain linear layer
wearing a graph model's clothes.

**Why it matters:** the whole point of the graph (per the proposal, §4.1.2,
and `graph_construction.md`) is that a cell line's driver mutations connect it
to specific proteins, which connect to the proteins a drug targets — that's
the mechanism the GNN is supposed to learn over. Without it, structural
information about *why* a given drug might work on a given cell line's
mutated pathway never reaches the cell-line embedding.

**Fix:** added the `("cell_line", "has_mutation", "protein")` edge (already
present in the saved graph, confirmed via `torch.load` — 4,341 edges) and its
reverse `("protein", "rev_has_mutation", "cell_line")` to both `conv1` and
`conv2`, mirroring how `05_train.py` already reverses the `drug→protein` edge.
`05_train.py` now also builds that reverse tensor at load time, the same way
it already did for `rev_targets`.

### 2. Protein nodes were fed all-zero features through a `Linear`

**File:** `hetero_gnn.py`

`prot_proj = nn.Linear(128, hidden_dim)` was applied to `x_dict["protein"]`.
Checked directly against the saved graph: every one of the 16,214 protein
node feature vectors is all-zero (`torch.all(x==0)` → `True`) — a documented
placeholder (see the docstring in the tracked `hetero_gnn.py`), not real
protein features. Projecting an all-zero vector through a `Linear` layer
produces the *same* output vector for every protein (only the bias term
survives), so every protein started message passing with an identical,
non-identifying embedding.

**Why it matters:** with every protein indistinguishable at layer 0, the only
signal `SAGEConv` had to differentiate proteins was aggregate neighbor count
— no notion of *which* protein a drug targets or *which* protein a mutation
hits could propagate. This silently defeated the entire PPI-network rationale
for including proteins as nodes at all.

**Fix:** replaced `prot_proj` with `nn.Embedding(num_proteins, hidden_dim)`,
matching the tracked `HeteroGNN`'s existing approach to the same problem.
`HeteroIC50GNN.__init__` now takes `num_proteins` as a required argument, and
`05_train.py` passes `x_dict["protein"].shape[0]`.

### 3. The train/val/test split leaked both cell lines and drugs

**File:** `dataset.py`

`train_test_split(df_valid, test_size=0.2, random_state=42)` — a plain random
split over every `(cell_line, drug)` pair. Because most cell lines and drugs
appear in many pairs, the same cell line (and the same drug) routinely landed
in both the training set and the test set, just paired with different
partners. The model could partly "look up" a cell line or drug it had already
seen elsewhere in training, rather than generalize to a truly new one.

**Why it matters:** this is exactly the leakage pattern `split_protocol_comparison.md`
documents inflating MoGraphDRP's published numbers relative to this project's
own `GroupShuffleSplit`-by-cell-line protocol. Every other tracked experiment
in the matrix (`experiment_utils.grouped_split`, 70/15/15,
`random_state=42`) uses a grouped split specifically to avoid this, so this
script's numbers were not comparable to anything else in the project even
before the other bugs.

**Fix:** replaced the split with a `GroupShuffleSplit` grouped by
`sanger_model_id`, 70/15/15, `random_state=42` — the same parameters
`experiment_utils.grouped_split` uses — plus the same post-split assertion
that no cell line appears in more than one split. Regenerating
`data/raw/aligned_ic50_pairs.csv` with this fix produced 111,799 matched
pairs (77,668 / 17,291 / 16,840 train/val/test), matching the project's
official fingerprint-resolvable population size exactly — a good sanity
check that the join logic itself was fine.

### 4. Fixed 30 epochs, checkpointed on the wrong metric

**File:** `05_train.py`

The training loop ran a hardcoded 30 epochs with no early stopping, and
saved the "best" checkpoint whenever validation **Pearson r** improved, not
validation RMSE.

**Why it matters:** after fixing bugs 1–3 and running once with this loop
unchanged (see Run 1 below), the log showed val RMSE bottoming out at
**epoch 3–7** and then climbing steadily through epoch 30 as train loss kept
falling — textbook overfitting. But val Pearson r kept *improving* well past
that point, because the model's predictions stayed well correlated in
direction while becoming increasingly miscalibrated in magnitude. Selecting
the checkpoint on r instead of RMSE meant the script kept a badly overfit
epoch instead of the actually-best one, and a fixed epoch count meant it had
no way to stop once the metric that matters started getting worse.

**Fix:** switched to `max_epochs=200` with early stopping (`patience=15`) and
checkpointing on best validation RMSE — the same criterion and patience
`experiments/GNN Ablation/gnn_baseline.py` already uses for the tracked
`HeteroGNN` runs, so the two are now trained under comparable stopping rules.

## Test runs

Both runs used the fixed graph edges, protein embedding, and grouped split
(bugs 1–3 fixed in both); they differ only in the training-loop fix (bug 4).

| Run | Checkpoint criterion | Stop | Best epoch | Test RMSE | Test Pearson r |
|---|---|---|---|---|---|
| 1 — before bug 4 fix | best val Pearson r | fixed 30 epochs | 26 | 1.5654 | 0.8780 |
| 2 — after bug 4 fix | best val RMSE, patience 15 | early-stopped @ 18 | 3 | **1.4050** | 0.8716 |

Run 2's log: val RMSE bottomed at epoch 3 (1.3537), never improved again by
more than the 1e-4 threshold over the next 15 epochs, and training stopped at
epoch 18 as designed.

### Comparison to the tracked matrix ([`results.md`](./results.md))

| Model | Omics | Drug rep | RMSE | PCC |
|---|---|---|---|---|
| MLP (E04) | Proteomics | fingerprint | 1.2843 | 0.8866 |
| **GNN-GCN (E12, tracked)** | GE+Mut_CNV+Proteomics, +mutation edges | fingerprint | **1.3513** | **0.8757** |
| **`HeteroIC50GNN`, fixed (this doc, Run 2)** | GE+Mut_CNV+Proteomics, +mutation edges | fingerprint | **1.4050** | **0.8716** |
| GNN-GCN, no mutation edges (E24) | GE+Mut_CNV+Proteomics | fingerprint | 1.4010 | 0.8731 |

Fixing all four issues brings this prototype from a broken, non-comparable
result into the same tier as the tracked `GNN-GCN` — but it still trails E12
by ~0.05 RMSE. The remaining gap is architectural, not a bug: E12's
`HeteroGNN` (`src/models/hetero_gnn.py`) uses a `ReduceLROnPlateau` LR
scheduler, `Adam` with `weight_decay=1e-5` (vs. plain `AdamW` at `1e-4`
here), and a deeper prediction head with `BatchNorm` and more dropout —
`HeteroIC50GNN` has none of those.

## Follow-up: closing the gap to E12 (Experiment A)

The ~0.05 RMSE gap above was attributable to three remaining differences
against the tracked `HeteroGNN`, all since applied to the test model:

1. **`ReduceLROnPlateau` scheduler** (`factor=0.5, patience=5`), stepped on
   val RMSE — matching `gnn_baseline.py::train`.
2. **Deeper prediction head** — `(256, 128)` with `BatchNorm1d` and
   `dropout=0.3`, replacing the shallow no-BatchNorm `(128, 64)` head.
3. **Weight decay** — made CLI-overridable and both arms tested.

`05_train.py` now also sets `drop_last=True` on the train loader so BatchNorm
never sees a 1-sample final batch.

| Run | Weight decay | Stop | Best epoch | Test RMSE | Test PCC |
|---|---|---|---|---|---|
| Baseline (bugs 1–4 fixed only) | 1e-4 | early stop @ 18 | 3 | 1.4050 | 0.8716 |
| + scheduler + deep head | 1e-5 | early stop @ 47 | 32 | 1.3368 | 0.8769 |
| + scheduler + deep head | **1e-4** | early stop @ 28 | 13 | **1.3035** | **0.8837** |

**The scheduler was the decisive change, not the weight decay.** The original
`1e-4` turned out to be the better setting; `1e-5` was worse in this
configuration. What mattered was that the LR schedule let training run to
epoch 28–47 with a meaningful best epoch (13–32), instead of peaking at epoch
3 and immediately overfitting as every prior run had.

### Which number is the comparable one: 1.3035 vs 1.3301

The 1.3035 above comes from `05_train.py`, which builds its own split in
`dataset.py`. That split uses the same *protocol* as the matrix (70/15/15
`GroupShuffleSplit` by cell line, `random_state=42`) and lands on the same
111,799-pair population — but it is computed independently, over a different
row ordering, so it is not the identical partition.

Re-running the same tuned model through the matrix's own shared code path
(`gnn_baseline.load_graph_and_pairs` -> `experiment_utils.grouped_split` ->
`experiment_utils.evaluate`) gives **RMSE 1.3301 / PCC 0.8805**, via
[`experiments/GNN Ablation/hetero_ic50_gnn_matrix.py`](../experiments/GNN%20Ablation/hetero_ic50_gnn_matrix.py).

**1.3301 is the number to quote.** The 0.027 difference between the two is
split-luck, not model improvement, and only the shared-split run is comparable
to the rest of [`results.md`](./results.md). This is exactly why the row in
that table comes from the matrix-protocol script rather than from
`05_train.py`.

**Final (shared protocol): RMSE 1.3301 / PCC 0.8805.** That beats the tracked
GNN-GCN (E12, 1.3513) by 0.021, making this the best graph model in the
project — but it does *not* beat CrossAttention/fingerprint (E10, 1.3205) or
the flat MLP (E04, 1.2843), as the more favourable 1.3035 figure would have
suggested.

## Conclusion

The fixes take `src/models/test/hetero_gnn.py` from a broken, non-comparable
prototype to **the best-performing GNN in the project** on the cell-line
split (RMSE 1.3301, vs. tracked E12's 1.3513).

It still does not beat the flat MLP baseline (E04, 1.2843), nor
cross-attention on the same fingerprints (E10, 1.3205). That matters more than
the intra-GNN win: the proposal's stated success criterion (§3, line 177) is
beating a flat concatenation baseline on RMSE, and no graph or attention model
in the project currently clears that bar when using the SMILES-derived drug
representation the proposal requires.

See [`leave_drugs_out_results.md`](./leave_drugs_out_results.md) §GNN for the
companion result on unseen drugs, where the graph performs *worse* than every
flat baseline.

## Files changed

- `src/models/test/hetero_gnn.py` — bugs 1–2
- `src/models/test/dataset.py` — bug 3
- `src/models/test/05_train.py` — reverse-edge construction (bug 1), `num_proteins` plumbing (bug 2), training loop (bug 4)

No files outside `src/models/test/` were modified. `data/raw/aligned_ic50_pairs.csv`
was regenerated by the fixed `dataset.py` and is not tracked in git (`data/raw/` output).
