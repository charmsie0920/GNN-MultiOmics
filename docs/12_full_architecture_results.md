# Full Proposed Architecture — Cross-Attention + PPI Graph + Molecular Graph

**Scripts:** [`experiments/12_full_architecture/full_architecture_matrix.py`](../experiments/12_full_architecture/full_architecture_matrix.py), [`src/models/full_architecture.py`](../src/models/full_architecture.py)
**Raw results:** `experiments/12_full_architecture/full_architecture_results.csv`
**Companions:** [`07_gnn_ablation_results.md`](./07_gnn_ablation_results.md), [`11_molecular_graph_results.md`](./11_molecular_graph_results.md), [`results.md`](./results.md)
**Date:** 2026-09-22 · **Runtime:** 2 runs in 400s

## What this covers

Every prior model implements a *part* of the architecture in
[`Architecture Simplified.jpeg`](./Architecture%20Simplified.jpeg). This is the
first run of all of it at once:

```
GE / Mut_CNV / Proteomics (PCA-128 each)
        |
        +-- cross-attention fusion ------> cell_line node features (256)
SMILES  +-- molecular graph GCN ---------> drug node features (128)
        |
        +-- HeteroConv message passing over the PPI graph
        |     (protein<->protein 473,860 | drug->protein 683 | cell_line->protein 4,341)
        |
        +-- cell_line embedding ++ drug embedding --> MLP --> ln(IC50)
```

**Controlled comparison: E13** (GNN-GCN, tri-omics, fingerprint, +mutation
edges, RMSE 1.3513), which uses the same `HeteroGNN` over the same graph. Two
things change:

- `cell_line` features: raw concatenated PCA blocks (384) → cross-attention
  fused embedding (256)
- `drug` features: 2048-bit Morgan fingerprint → molecular graph GCN (128)

Everything else — graph, edge types, hidden dim, layer count, optimizer,
split, seed, batch size — is pinned to
[`experiments/07_gnn_ablation/gnn_baseline.py`](../experiments/07_gnn_ablation/gnn_baseline.py).

Population: 111,799 pairs, identical to E13.

## Node ordering

The GNN indexes `cell_line` and `drug` embeddings by the graph's own node
indices, but the PCA CSVs and the SMILES cache have their own orderings. A
silent mismatch would pair one cell line's biology with another's label — and
the model would still train, still converge, and still produce a plausible
RMSE.

`align_omics_to_graph` and `align_drug_graphs_to_graph` reorder explicitly and
assert coverage. The module smoke test perturbs one cell line's omics and
confirms **exactly one** fused row changes.

## Results

| Model | Test RMSE | MAE | R² | PCC | SCC | AUC | F1 | Params | Best epoch | Fit (s) |
|---|---|---|---|---|---|---|---|---|---|---|
| **FullArchitecture-GCN** | **1.3397** | 1.0062 | 0.7649 | 0.8774 | 0.8309 | 0.9087 | 0.8337 | 3,319,425 | 65 | 267 |
| FullArchitecture-GAT | 2.3269 | 1.6907 | 0.2907 | 0.5500 | 0.5111 | 0.7592 | 0.6485 | 3,321,985 | 10 | 133 |

Validation: GCN 1.2705 (R² 0.8038), GAT 2.3765 (R² 0.3136).
Mean-only floor: 2.7690.

## 1. The encoder swaps help the graph model

| | Test RMSE |
|---|---|
| E13 — GNN-GCN, concat PCA + fingerprint | 1.3513 |
| **FullArchitecture-GCN** — fusion + molecular graph | **1.3397** |

An improvement of 0.0116. Cross-attention fusion and molecular graphs are
better inputs to the heterogeneous GNN than concatenated PCA blocks and Morgan
fingerprints. The gap is small and untested across seeds, so it should be
reported as directional rather than established.

## 2. The PPI graph does not pay for itself

The more important comparison is against the *same* fusion and the *same*
molecular graph **without** the PPI branch
([`11_molecular_graph_results.md`](./11_molecular_graph_results.md)):

| | Val RMSE | Test RMSE |
|---|---|---|
| Molecular graph, no PPI graph | 1.2679 | **1.3021** |
| Full architecture (adds PPI message passing) | 1.2705 | 1.3397 |

Adding message passing over the biological graph **costs 0.038 RMSE**. Note the
two models are 0.003 apart on validation and 0.038 apart on test — a gap that
opens only at test time, which is characteristic of unstable single-seed
results. Given the ±0.015 seed spread measured in
[`13_seed_variance_results.md`](./13_seed_variance_results.md) this is around
2.5× the noise: probably real, not confirmed.

The honest statement: **the PPI graph branch did not improve performance and
appears to cost roughly 0.04 RMSE under a cell-line-grouped split.**

This does not mean biological networks are useless for drug response. It means
that *this* graph — STRING PPI edges, driver-mutation cell-line links,
zero-initialized protein embeddings — does not add information the PCA-compressed
omics did not already carry.

## 3. GAT fails, for the third independent time

R² 0.29 against GCN's 0.76, with training effectively stalling at epoch 10.

This is the third GAT failure in the project ([E57](./07_gnn_ablation_results.md)
2.2546, [E59](./07_gnn_ablation_results.md) 2.3558, now 2.3269), and MoGraphDRP's
own Table 3 independently ranks GAT (1.0311) below GCN (0.9497), attributing
the gap to "overfitting in noisy biological data."

A mechanism was visible before training: during the runner's smoke test, GAT's
gradient into the drug GCN measured **0.34 against GCN's 8.97** after three
steps — a 27× difference. The attention appears to starve the drug branch.
This is a single observation on an untrained model and is offered as a
hypothesis, not a finding.

**Recommendation: formally deviate from the proposal's GAT choice.** Three
independent runs plus the benchmark's own ablation is sufficient evidence, and
a documented deviation with evidence is a stronger report position than an
unexplained gap between proposal and implementation.

## Caveats

- Single seed (42). Both comparisons above are within ~2.5× the measured
  seed spread; neither has been repeated.
- Full-graph message passing runs over all 473,860 PPI edges every step, which
  is why `BATCH_SIZE` is 1024 here against 256 for the flat models.
- `cell_line` nodes are only reachable through the 4,341 driver-mutation edges
  added by [`04_link_cell_lines.py`](../src/data/04_link_cell_lines.py). The
  runner refuses to run without them — without those edges `cell_line` is a
  disconnected component and message passing cannot reach it.
- Only tri-omics was run. Every prior arm found GE+Proteomics beats tri-omics,
  and that subset remains unmeasured here. Note the PPI graph carries mutation
  information through its edges regardless, so dropping `Mut_CNV` from the
  fusion does not remove genomics from this model — it changes the channel.
