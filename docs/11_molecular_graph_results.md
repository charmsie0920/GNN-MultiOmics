# Molecular Graph Drug Representation — Atom-Level GCN vs Morgan Fingerprint

**Scripts:** [`experiments/11_molecular_graph/molecular_graph_matrix.py`](../experiments/11_molecular_graph/molecular_graph_matrix.py), [`src/data/drug_graphs.py`](../src/data/drug_graphs.py), [`src/models/drug_gcn.py`](../src/models/drug_gcn.py)
**Raw results:** `experiments/11_molecular_graph/molecular_graph_matrix_results.csv`
**Companions:** [`06_cross_attention_ablation_results.md`](./06_cross_attention_ablation_results.md), [`13_seed_variance_results.md`](./13_seed_variance_results.md), [`results.md`](./results.md)
**Date:** 2026-09-22

## What this covers

The third drug representation in the project, after one-hot identity and the
2048-bit Morgan fingerprint: the drug as a **graph of atoms and bonds**,
encoded by a 3-layer GCN with mean pooling. This is the representation
MoGraphDRP uses (their §2.2.2), and the component their Table 6 attributes the
largest single effect to (`without graph` = RMSE 2.5642).

A fingerprint enumerates *which* substructures are present; a molecular graph
keeps *how the atoms connect*. The question is whether that topology is worth
anything for drug response prediction.

Atom features (39-dim) follow the descriptor set the benchmark names — element
type, degree, hydrogen count, aromaticity — plus hybridization and ring
membership. Across the 498 resolvable GDSC drugs this yields 15,762 atoms and
34,484 directed bond edges; the 240 drugs that actually appear in GDSC2
response pairs account for 7,575 atoms.

**Controlled comparison: E10** (CrossAttention, GE+Proteomics, fingerprint).
Identical fusion, omics, head, split and seed — only the drug encoder changes.
Population parity is asserted before training: both arms resolve exactly the
same 498 drugs and the same 111,799 pairs.

One-hot arms are deliberately absent. Drug identity cannot generalize to unseen
compounds (R² 0.024 in [`10_leave_drugs_out_results.md`](./10_leave_drugs_out_results.md)),
so it is not a representation this project builds on.

## Encode-once design

There are 240 distinct drugs but 111,799 pairs. Re-encoding a molecule per pair
would repeat the same computation thousands of times per epoch, so all
molecular graphs are collated into **one disjoint graph** and encoded once per
forward pass; the resulting `(n_drugs, 128)` table is then indexed by an
integer drug code.

The code is a lookup key, never a feature — it does not reach the prediction
head. The drug is represented purely by its structure, which is what makes this
arm distinct from the one-hot arm rather than a disguised version of it.

Because the collated graph is block-diagonal, this is mathematically identical
to encoding each molecule separately. That equivalence is asserted in
[`src/models/drug_gcn.py`](../src/models/drug_gcn.py)'s smoke test — if it ever
failed, molecules would leak structure into each other and every number below
would be silently wrong.

## Results

All rows: CrossAttention, GE+Proteomics, molecular graph, 111,799 pairs,
501,377 parameters, seed 42.

| Variant | Val RMSE | Test RMSE | Test R² | Test PCC | Best epoch | Fit (s) |
|---|---|---|---|---|---|---|
| `n_tokens=1` (first run) | 1.2679 | 1.3021 | 0.7779 | 0.8837 | 46 | 273 |
| `n_tokens=1` (repeat) | 1.2738 | 1.3099 | 0.7752 | 0.8825 | 50 | 300 |
| `n_tokens=2` | 1.2839 | 1.3190 | 0.7721 | 0.8798 | 37 | 246 |
| `n_tokens=4` | **1.2632** | 1.3173 | 0.7727 | 0.8812 | 85 | 475 |
| `n_tokens=8` | 1.2723 | 1.3131 | 0.7741 | 0.8828 | 70 | 407 |
| `--standardize-targets` | **1.2595** | 1.3331 | 0.7672 | 0.8802 | 23 | 159 |

Mean-only floor: **2.7690**. E10 (fingerprint) as recorded: **1.3205**.

> **The CSV holds only the most recent invocation.** Each run overwrites
> `molecular_graph_matrix_results.csv`, so the table above is transcribed from
> the run logs. Re-running any variant regenerates only that row.

## 1. Molecular graphs do not beat fingerprints

The first run (1.3021 vs E10's recorded 1.3205) looked like a 1.4% improvement.
It was not. Five-seed repeats in
[`13_seed_variance_results.md`](./13_seed_variance_results.md) give:

| Arm | Mean test RMSE | Std | Seeds won |
|---|---|---|---|
| Fingerprint (E10 config) | **1.3029** | 0.0126 | 4 / 5 |
| Molecular graph | 1.3194 | 0.0153 | 1 / 5 |

Paired per-seed difference: **+0.0165 ± 0.0231** in favour of the fingerprint.
The two arms are not separable at n=5, but the point estimate favours the
*simpler* representation, and the apparent win in the single-seed run was noise
pointing the wrong way.

What survives is a parameter claim: the molecular graph matches a Morgan
fingerprint using **501,377 parameters against 742,017** — the GCN encoder is
38,144 parameters versus the `FingerprintEncoder`'s 278,784. Equivalent
accuracy at 32% fewer parameters is a real, if modest, result.

## 2. Why topology buys nothing here — and where it might

Under the cell-line-grouped split, **71% of the variance in ln(IC50) is drug
identity** (see §4). Every test drug appears in training with thousands of
measurements, so the model can memorize each compound's response profile and
bank most of its accuracy there. The drug *representation* is nearly
irrelevant: there is nothing for structure to do.

This predicts that structure should matter on the axis where memorization is
impossible — unseen compounds. That is the leave-drugs-out protocol, where
one-hot collapses to R² 0.024 while fingerprints hold at 0.422. Whether
molecular graphs improve on fingerprints there is the open question this
experiment does **not** answer.

## 3. Making the attention non-degenerate changed nothing

`MultiOmicsCrossAttentionFusion` with `n_tokens=1` gives each cross-attention
call exactly one query and one key, so softmax returns 1.0 by construction —
the block projects, it does not attend. `n_tokens > 1` splits each modality's
128-dim PCA vector into contiguous chunks so the attention distribution is
genuinely learned.

Test RMSE across token counts spans **1.3099–1.3190, a range of 0.009** —
smaller than the difference between two runs of the *identical* configuration
(0.0078), and non-monotonic (ordering 1, 8, 4, 2). **No detectable effect.**

Higher token counts did produce the best validation scores while widening the
validation–test gap (`n_tokens=4`: best-ever val 1.2632, val→test gap 0.054 vs
0.036 at `n_tokens=1`), indicating the added flexibility fits the validation
cell lines rather than generalizing.

This is consistent with MoGraphDRP §2.1, which reports that direct
concatenation outperformed learned fusion weights in their own experiments.
**Caveat:** token count also reduces attention projection width (501,377 →
371,009 parameters), so granularity and capacity move together and a null
result could in principle be two effects cancelling.

`n_tokens=1` is retained as the default and is bit-identical to the
pre-existing implementation, so every recorded cross-attention result in
[`results.md`](./results.md) remains valid.

## 4. Per-drug target standardization: faster, not better

Standardizing ln(IC50) within each drug (training-set statistics only, inverted
before scoring so metrics stay in ln(IC50) units) converged in **23 epochs
against 39–55**, a 40% compute saving, and produced the best validation RMSE of
any run in this experiment (1.2595). Test RMSE (1.3331) sits inside the
unstandardized seed range (1.3044–1.3408): **no detectable effect.**

The diagnostic that motivated it is worth more than the result. Decomposing the
training target:

| Baseline | Test RMSE |
|---|---|
| Global training mean | 2.7690 |
| **Per-drug training mean** (no omics at all) | **1.4889** |
| Best model in this experiment | 1.3021 |

**Drug identity alone accounts for 71% of the reducible error.** A model with
no omics input whatsoever, predicting each compound's average response, reaches
1.4889. The multi-omics pipeline improves on that by roughly a quarter of the
remaining variance — a much more honest framing than "R² 0.78" against the
global-mean floor, and one that explains why one-hot drug identity scored so
well (it encodes that 71% directly) and why every architecture in this project
clusters within ~0.03 RMSE.

**The per-drug mean belongs in [`results.md`](./results.md) as a baseline row.**

## Caveats

- Every row above is a single seed except where
  [`13_seed_variance_results.md`](./13_seed_variance_results.md) is cited.
  Seed-to-seed spread for this configuration is ±0.015, and cross-session
  spread (different hardware, different PyTorch build) is larger still — E10's
  recorded 1.3205 re-measures as 1.2997 at the same seed. Treat any difference
  under ~0.03 as undetectable.
- GPU training is not bit-reproducible: two runs of the identical
  configuration at seed 42 gave 1.3021 and 1.3099. `torch.manual_seed` fixes
  initialization and shuffling, not CUDA reduction order.
- 123 of 621 GDSC drug IDs have no resolvable SMILES and are dropped, as in the
  fingerprint arm. Parity is asserted, not assumed.
- Only the GE+Proteomics subset was run. The tri-omics molecular-graph cell
  remains unmeasured.
