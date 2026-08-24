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
