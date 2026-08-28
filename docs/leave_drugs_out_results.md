# Leave-Drugs-Out — The Experiment That Justifies Morgan Fingerprints

**Script:** [`experiments/Leave Drugs Out/leave_drugs_out.py`](../experiments/Leave%20Drugs%20Out/leave_drugs_out.py)
**Raw results:** `experiments/Leave Drugs Out/leave_drugs_out_results.csv`
**Companions:** [`results.md`](./results.md), [`cross_attention_ablation_results.md`](./cross_attention_ablation_results.md), [`split_protocol_comparison.md`](./split_protocol_comparison.md)
**Date:** 2026-08-24 · **Runtime:** 6 runs

## The problem this solves

The main matrix showed Morgan fingerprints **hurting** the best model
(cross-attention GE+Proteomics: 1.2847 one-hot vs 1.3205 fingerprint — see
[`cross_attention_ablation_results.md`](./cross_attention_ablation_results.md)
§3). That is an awkward result to defend when the proposal's §4.1.2 specifies
2048-bit RDKit Morgan fingerprints for the final architecture.

The reason is that **every experiment so far held out cell lines, never drugs**.
Under that protocol a test drug has always been screened during training, so
one-hot identity works perfectly well — it looks up a drug the model already
knows. The comparison never tested the axis where the two representations
actually differ.

This experiment holds out **drugs** instead: 36 of 240 compounds never appear in
training (`drug_overlap = 0`, asserted in code), and the model must score them
cold.

| | One-hot identity | Morgan fingerprint |
|---|---|---|
| Held-out drug's representation | A column never activated in training — its weight stays at initialization | 2048 structural bits computed from SMILES |
| Can the model use it? | **No.** Not "predicts poorly" — structurally cannot represent the compound | Yes, if it shares substructure with training compounds |
| Falls back on | Cell-line features alone | Structural similarity transfer |

Both arms run on the fingerprint-resolvable population (111,799 pairs) so the
comparison is not confounded by differing row sets — the same control used
throughout the main matrix. Omics is fixed at GE+Proteomics, the best subset
from [`results.md`](./results.md) (E01).

## Results

Mean-only floor: **2.5500**. Train 168 drugs / val 36 / test 36, zero overlap.

| Model | Drug rep | RMSE | MAE | R² | PCC | SCC | AUC | F1 |
|---|---|---|---|---|---|---|---|---|
| RF | onehot | 2.3638 | 1.8325 | 0.1242 | 0.3533 | 0.3493 | 0.6608 | 0.5984 |
| RF | **fingerprint** | **1.8719** | **1.5328** | **0.4508** | **0.6733** | **0.5305** | **0.7382** | **0.6352** |
| MLP | onehot | 2.4298 | 1.9103 | 0.0746 | 0.3233 | 0.3207 | 0.6482 | 0.4292 |
| MLP | **fingerprint** | **1.8516** | **1.4777** | **0.4626** | **0.6883** | **0.5631** | **0.7600** | **0.6868** |
| CrossAttention | onehot | 2.4959 | 1.9809 | 0.0236 | 0.3308 | 0.3272 | 0.6511 | 0.3928 |
| CrossAttention | **fingerprint** | **1.9200** | **1.5459** | **0.4222** | **0.6618** | **0.5258** | **0.7396** | **0.6733** |

## 1. One-hot does not merely underperform on unseen drugs — it collapses

Measured as improvement over the mean-only floor (2.5500), which is what a model
with *no usable information* would achieve:

| Model | one-hot beats floor by | R² | fingerprint beats floor by | R² |
|---|---|---|---|---|
| RF | +0.186 | 0.124 | **+0.678** | **0.451** |
| MLP | +0.120 | 0.075 | **+0.698** | **0.463** |
| CrossAttention | **+0.054** | **0.024** | **+0.630** | **0.422** |

Cross-attention with one-hot reaches **R² = 0.024** — essentially
indistinguishable from predicting the training mean. Its F1 of 0.3928 is the
worst classification score anywhere in the project.

The small residual signal it does retain comes entirely from the **cell-line**
side: some cell lines are broadly drug-sensitive or broadly resistant, and the
model can still learn that. The drug side contributes nothing, exactly as the
mechanism predicts — there is no trained weight for a one-hot column that was
never activated.

Note the collapse is *worst* for the most expressive model. RF retains the most
one-hot signal (R² 0.124), cross-attention the least (R² 0.024), which is
consistent with the stronger models relying more heavily on drug identity when
it is available, and therefore having more to lose when it disappears.

## 2. Fingerprints transfer, consistently, across all three architectures

R² 0.42–0.46 and PCC 0.66–0.69 on compounds never screened during training. The
gain over one-hot is stable at **+0.49 to +0.58 RMSE** regardless of model
family, which is what you would expect from a representation effect rather than
an architecture-specific artifact.

This is the generalization argument made concrete: a held-out compound sharing
substructures with training drugs inherits what the model learned about them.
The model predicts *"this looks like an EGFR inhibitor, and EGFR inhibitors
behave this way on this cell line"* rather than needing to have seen the exact
molecule.

## 3. The trade-off, now measured rather than asserted

| Protocol | Task | one-hot | fingerprint | Winner |
|---|---|---|---|---|
| Leave-cell-lines-out | New cell line, **known** drug | **1.2442** | 1.3205 | one-hot by 0.076 |
| Leave-drugs-out | Known cell line, **new** drug | 2.4959 | **1.9200** | **fingerprint by 0.576** |

Fingerprints cost **~0.08 RMSE** in-distribution and buy **~0.58 RMSE** on
unseen compounds — roughly a **7× return** on the trade.

**This is the experimental justification for the proposal's §4.1.2 fingerprint
choice.** It also resolves the tension in
[`cross_attention_ablation_results.md`](./cross_attention_ablation_results.md)
§3, which could only note that fingerprints looked worse and argue theoretically
that unseen-drug generalization was the missing test. It was, and this is it.

## What to report

Present both protocols together. The one-line version:

> One-hot drug identity slightly outperforms Morgan fingerprints when the test
> drug was seen during training (RMSE 1.2442 vs 1.3205), but collapses to near
> the mean-only baseline on compounds never screened (R² 0.024 vs 0.422).
> Fingerprints are therefore retained in the final architecture: the small
> in-distribution cost buys the ability to score novel compounds, which one-hot
> encoding structurally cannot do.

Reporting only the cell-line-grouped comparison would make the fingerprint
choice look unjustified. Reporting only this one would overstate its
in-distribution value. Both are needed.

## GNN — does the graph help on unseen drugs?

**Script:** [`experiments/Leave Drugs Out/leave_drugs_out_gnn.py`](../experiments/Leave%20Drugs%20Out/leave_drugs_out_gnn.py)
**Raw results:** `experiments/Leave Drugs Out/leave_drugs_out_gnn_results.csv`
**Date:** 2026-08-28

The three architectures above are all flat models over a concatenated feature
vector. None of them use the `cell_line`/`drug`/`protein` graph. Since the
whole rationale for the graph is that PPI topology and drug–target edges carry
biological signal a flat vector cannot, the obvious question is whether the
graph helps on precisely the axis it should help most: scoring compounds never
seen in training.

It does not. Running the tracked `HeteroGNN` (variant `gcn`, with mutation
edges — the E12 configuration, RMSE 1.3513 on the cell-line split) through
`leave_drugs_out_split`:

| Model | Drug rep | Test RMSE | R² | PCC |
|---|---|---|---|---|
| MLP | fingerprint | **1.8516** | **0.463** | — |
| RF | fingerprint | 1.8719 | 0.451 | — |
| CrossAttention | fingerprint | 1.9200 | 0.422 | — |
| HeteroIC50GNN (tuned) | fingerprint | 1.9704 | 0.391 | 0.658 |
| **GNN-GCN (tracked, E12)** | fingerprint | **2.0636** | **0.333** | 0.628 |
| *mean-only floor* | — | *2.5500* | — | — |

Both GNNs land **below every flat baseline**. The tracked GNN-GCN is worst,
0.21 RMSE behind the flat MLP on identical Morgan fingerprints, clearing the
mean-only floor by only 0.49 (vs the MLP's 0.70). There is no one-hot arm:
the graph's drug nodes carry fingerprint features by construction, so this is
inherently fingerprint-only.

The second GNN row is the tuned `HeteroIC50GNN`
(`src/models/test/hetero_gnn.py`, script:
[`leave_drugs_out_gnn_test_model.py`](../experiments/Leave%20Drugs%20Out/leave_drugs_out_gnn_test_model.py)),
included because it *beats* the tracked HeteroGNN on the cell-line split
(1.3301 vs 1.3513 — see
[`hetero_gnn_test_bugfixes.md`](./hetero_gnn_test_bugfixes.md)). That
advantage does carry over here (+0.093 RMSE, +0.058 R² over E12), so it is
genuinely the better GNN on both protocols — but it is still 0.12 RMSE behind
the plain MLP and does not overturn the finding.

**That the gap survived a deliberate tuning effort is itself the evidence.**
If the GNN were merely undertrained, adding an LR scheduler, a deeper
BatchNorm head, and proper early stopping — which bought 0.048 RMSE on the
cell-line split — should have closed some of it. It did not. The limitation
is in what the graph carries, not in how the model is fit.

**The result is stronger than it looks, because the setup favours the GNN.**
This is full-batch transductive message passing over one fixed graph, so a
held-out drug's node — its fingerprint *and* its `drug→protein` target edges —
is present throughout training; only its IC50 labels are withheld from the
loss. The MLP and RF are strictly inductive by comparison, seeing a held-out
drug's features for the first time at test time. The GNN had a structural
advantage and still lost. A like-for-like inductive comparison would require
masking held-out drug nodes out of the graph during training, and would likely
be *worse* for the GNN, not better.

### Why the graph probably isn't contributing

Two properties of the constructed graph are the more plausible culprits than
the model architecture:

- **Protein nodes carry no features.** All 16,214 protein feature vectors are
  all-zero placeholders (`03_graph_construction.py` zero-fills them);
  `HeteroGNN` substitutes a learnable `nn.Embedding`, so protein identity is
  learned purely from topology with no biological prior attached.
- **Drug–target edges are extremely sparse.** 683 `drug→targets→protein`
  edges across 498 drug nodes — averaging ~1.4 known targets per compound,
  and many drugs likely have none. A held-out drug with no target edge is
  connected to the rest of the graph by nothing at all, leaving the GNN with
  only its fingerprint — the same information the MLP has, routed through more
  parameters and more opportunity to overfit.

A held-out drug with no target edge is connected to the rest of the graph by
nothing at all — so the GNN has exactly the information the MLP has (the
fingerprint), routed through 2.8M parameters with more opportunity to overfit.
No amount of tuning manufactures signal the graph does not carry, which is
what the tuned model's result above demonstrates.

Combined with both GNNs also trailing the flat MLP on the cell-line split
(1.3301 for the best-tuned variant, vs E04's 1.2843 — see
[`hetero_gnn_test_bugfixes.md`](./hetero_gnn_test_bugfixes.md)), the graph is
currently not earning its place on **either** axis. Improving graph
construction (real protein features, denser drug–target coverage) is a more
promising direction than further tuning the GNN architecture.

## Caveats

- **Absolute numbers here are not comparable to the main matrix.** Predicting
  for an unseen *drug* is a harder task than predicting for an unseen *cell
  line*, and the mean-only floor differs (2.5500 vs 2.7097) because the
  population and split differ. Compare within this document, not across.
- **Cell lines leak here by design.** This protocol holds out drugs, so a test
  cell line does appear in training. That is correct for the question being
  asked (can we score a new compound on cell lines we have profiled?), but it
  means these numbers do not measure the doubly-hard case of a new drug on a
  new cell line. A leave-both-out split remains unrun.
- Only 36 test drugs, so the test-set estimate is noisier than the main
  matrix's 80 held-out cell lines. Directionally the result is unambiguous —
  a 0.49–0.58 RMSE gap reproduced across three independent architectures — but
  the precise magnitudes should not be over-read.
- Omics fixed at GE+Proteomics; the omics-subset axis was not re-swept under
  this protocol.
