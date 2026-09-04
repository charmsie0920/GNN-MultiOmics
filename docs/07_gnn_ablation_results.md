# GNN Ablation — Architecture × Cell-Line Graph Linkage

**Scripts:** [`experiments/07_gnn_ablation/gnn_baseline.py`](../experiments/07_gnn_ablation/gnn_baseline.py), [`src/models/hetero_gnn.py`](../src/models/hetero_gnn.py), [`src/data/04_link_cell_lines.py`](../src/data/04_link_cell_lines.py)
**Raw results:** `experiments/07_gnn_ablation/gnn_results.csv`
**Companions:** [`06_rf_ablation_results.md`](./06_rf_ablation_results.md), [`06_mlp_ablation_results.md`](./06_mlp_ablation_results.md), [`06_cross_attention_ablation_results.md`](./06_cross_attention_ablation_results.md), [`results.md`](./results.md)
**Date:** 2026-08-24 · **Runtime:** 4 runs in 669.3s

## What this covers

The first models in the project to use **relational structure**. Two axes,
4 cells: {GCN, GAT} × {with, without cell_line→protein edges}. This is the
project's analogue of the MoGraphDRP paper's two graph ablations (their Table 3
GNN-architecture sweep, and their Table 6 "without graph" row).

Omics is fixed at tri-omics (the graph's `cell_line` nodes carry all three PCA
blocks concatenated, 384-dim) — the omics-subset axis is covered by the other
three model families. Drug representation is inherently the 2048-bit Morgan
fingerprint, since that is what the graph's drug nodes store; there is no
one-hot variant.

## Prerequisite: the graph had no cell-line edges

As built by
[`03_graph_construction.py`](../src/data/03_graph_construction.py) and
documented in [`graph_construction.md`](./graph_construction.md), the
`HeteroData` object contained `(protein, interacts_with, protein)` and
`(drug, targets, protein)` edges — but **nothing incident to `cell_line`
nodes**. Cell lines were a disconnected component, so message passing could
not move any PPI or drug-target signal into a cell-line embedding. A GNN over
that graph would have been a GNN in name only.

[`04_link_cell_lines.py`](../src/data/04_link_cell_lines.py) closes this by
adding `(cell_line, has_mutation, protein)` edges from **cancer-driver
mutations**:

```
[mutations] scanned 14,370,318 rows -> 4,341 (cell_line, gene) pairs
            across 526 cell lines, 537 distinct genes (driver_only=True)
[aliases]   537 / 537 gene symbols resolved to an ENSP
[protein]   added 11 protein nodes referenced only by mutations
[has_mutation] 4,341 edges linking 526 / 532 cell lines
```

**Why driver mutations only:** the mutation table averages ~6,300 mutations per
cell line, overwhelmingly intronic/non-coding. Linking all of them would add
millions of mostly-noise edges and swamp the 473,860 high-confidence PPI edges.
The `cancer_driver` flag (~0.3% of rows) is the biologically motivated sparse
subset and creates exactly the intended multi-hop path:

```
cell_line --has_mutation--> protein --interacts_with--> protein <--targets-- drug
```

These edges derive from mutation status only, never from IC50, so they
introduce **no label leakage**. Six cell lines carry no driver mutation at all
and remain unlinked — expected, not a defect.

## Model

`HeteroConv` wrapping one conv per edge type, every relation also mirrored so
messages flow both directions. 2 layers, hidden dim 128, LayerNorm + ReLU +
dropout 0.2 between layers, then a `[256,128]` MLP head over
`concat(cell_line_embedding, drug_embedding)`.

Protein nodes carry zero-filled placeholder features from graph construction;
here they become a learnable `nn.Embedding(16214, 128)` — the step
[`graph_construction.md`](./graph_construction.md) flagged as the GNN stage's
job.

The "GCN" variant uses `SAGEConv` rather than `GCNConv`: `GCNConv` normalizes
by degree over a single node set and requires a symmetric same-type graph, but
our `cell_line→protein` and `drug→protein` relations are bipartite. `SAGEConv`
handles `(src, dst)` tuple dims natively.

## Results

| Model | Mutation edges | Pairs | RMSE | MAE | R² | PCC | SCC | AUC | F1 | Params |
|---|---|---|---|---|---|---|---|---|---|---|
| **GNN-GCN** | **yes** | 111,799 | **1.3513** | **1.0128** | **0.7608** | **0.8757** | **0.8280** | **0.9081** | **0.8221** | 2,816,257 |
| GNN-GCN | no | 111,799 | 1.4010 | 1.0530 | 0.7429 | 0.8731 | 0.8236 | 0.9078 | 0.8195 | 2,684,673 |
| GNN-GAT | no | 111,799 | 2.2546 | 1.7285 | 0.3341 | 0.5865 | 0.5504 | 0.7857 | 0.6612 | 2,686,209 |
| GNN-GAT | yes | 111,799 | 2.3558 | 1.8151 | 0.2730 | 0.5359 | 0.5033 | 0.7482 | 0.6206 | 2,818,817 |

Mean-only floor: **2.7690** (111,799-pair population).

## 1. The cell-line linkage works — for GCN

GCN improves **1.4010 → 1.3513** (−0.050 RMSE, R² 0.7429 → 0.7608) when the
driver-mutation edges are present. That is a real gain from 4,341 edges — about
0.9% of the graph's total edge count — and it confirms the linkage is carrying
signal rather than being decorative.

It also mirrors the paper's own strongest ablation finding: their Table 6
"without graph" row collapsed to RMSE 2.5642 (from 0.6622), by far their
largest single-component effect. Direction agrees; magnitude here is far
smaller because our "without" arm still has a PPI graph, just one the cell
lines cannot reach.

## 2. GAT fails badly, and the linkage makes it worse

GAT lands at **2.2546 / PCC 0.5865** — worse than the RF baseline and barely
better than the 2.7690 mean-only floor. Adding the mutation edges *degrades* it
further (2.3558, PCC 0.5359), the opposite of GCN's response, at near-identical
parameter count (2.69M vs 2.68M).

**This independently reproduces the paper's Table 3 finding.** They tested four
drug-graph encoders and found GCN best (0.9497) with GAT notably worse
(1.0311), attributing it to GAT overfitting noisy biological data. We see the
same ordering on a completely different graph topology — a 16K-node PPI network
rather than small molecular graphs — which strengthens the result: it is not
specific to molecule-shaped graphs.

A plausible mechanism for the degradation-with-more-edges: GAT learns
per-edge attention coefficients, so the 4,341 new heterogeneous edges add
attention parameters that must be fit from the same data. On a graph where
`interacts_with` already supplies 473,860 edges of highly variable
informativeness, learned edge weighting appears to overfit where SAGE's uniform
neighborhood aggregation stays robust.

**Direct consequence for the proposal:** `MCS16 - Final Proposal Report.docx`
§4.2 specifies **GAT** for the final architecture. On this data, that choice is
worse than GCN by 0.90 RMSE. This should be raised with the team before the
architecture is frozen — either switch to GCN/SAGE, or treat GAT's
underperformance as itself an experimental finding worth reporting.

## 3. GNN does not beat cross-attention

Best GNN (1.3513) versus the flat-model families on the matched 111,799-pair
population:

| Model | Best config | RMSE |
|---|---|---|
| Cross-attention | GE+Proteomics, onehot_restricted | **1.2847** |
| MLP | Proteomics, fingerprint | 1.2843 |
| **GNN-GCN** | tri-omics, +mutation edges | 1.3513 |
| RF | GE+Proteomics, fingerprint | 1.3712 |

The GNN sits between RF and the neural flat models. Relational structure helps
relative to trees, but on this data it does not yet beat learned fusion over
flat features.

Two honest caveats before concluding "graphs don't help here":
- The GNN is constrained to the fingerprint population (111,799 pairs) and
  tri-omics, whereas cross-attention's best cell uses one-hot on the full
  134,764 pairs and drops Mut_CNV. Matched-population comparison above uses
  `onehot_restricted` to be fair.
- Only 4 GNN configurations were tried, versus 12 for cross-attention, with no
  hyperparameter search (layers, hidden dim, and dropout were fixed at first
  reasonable values). The GNN is the least-tuned family in the matrix.

## Caveats

- **Drug-target edges are extremely sparse** — 683 edges over 498 drug nodes,
  because only 266 of 380 GDSC target symbols resolved to an ENSP (see
  [`graph_construction.md`](./graph_construction.md) §4). Many drugs have no
  path into the protein subgraph at all, limiting how much the PPI network can
  inform their embeddings.
- Full-batch message passing over the whole graph each step; no neighbor
  sampling. Fine at this scale (~20MB of tensors) but would not extend to a
  substantially larger graph.
- Protein embeddings are learned from scratch with no biological
  initialization (no sequence or functional features), so the model must infer
  16,214 protein representations from IC50 supervision alone.
