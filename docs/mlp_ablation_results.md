# MLP Ablation — Full Omics × Drug-Representation Matrix

**Script:** [`experiments/Full Matrix/mlp_matrix.py`](../experiments/Full%20Matrix/mlp_matrix.py)
**Raw results:** `experiments/Full Matrix/mlp_matrix_results.csv`
**Companions:** [`rf_ablation_results.md`](./rf_ablation_results.md), [`cross_attention_ablation_results.md`](./cross_attention_ablation_results.md), [`gnn_ablation_results.md`](./gnn_ablation_results.md), [`results.md`](./results.md)
**Date:** 2026-08-24 · **Runtime:** 21 runs in 903.7s

## What this covers

The same 21 cells as [`rf_ablation_results.md`](./rf_ablation_results.md) —
7 omics subsets × 3 drug-representation arms — with the learner swapped from a
tree ensemble to a feed-forward network. Everything else is held fixed
(identical `random_state=42` cell-line-grouped 70/15/15 split, identical
flat-concatenation input, identical metrics), so **RF vs MLP at a fixed cell
isolates "gradient-trained network vs. tree ensemble" cleanly.**

That comparison is the control
[`cross_attention_fusion_results.md`](./cross_attention_fusion_results.md) §2
called for and never got: it separates "neural nets fit this problem better"
from "learned fusion helps."

Model: `[Linear → BatchNorm → ReLU → Dropout] × 3` with hidden dims
`[512, 256, 128]`, dropout 0.3, Adam (lr 1e-3, weight decay 1e-5),
`ReduceLROnPlateau`, early stopping on val RMSE (patience 15, max 200 epochs),
batch 256, `torch.manual_seed(42)`.

See [`rf_ablation_results.md`](./rf_ablation_results.md) for why the three drug
arms exist (`onehot` = 134,764 pairs; `onehot_restricted` and `fingerprint` =
111,799 pairs, the SMILES-resolvable subset).

## Results

| Omics | Drug rep | Pairs | RMSE | MAE | R² | PCC | SCC | AUC | F1 |
|---|---|---|---|---|---|---|---|---|---|
| GE | onehot | 134,764 | 1.3995 | 1.0727 | 0.7321 | 0.8602 | 0.8122 | 0.8990 | 0.8161 |
| GE | onehot_restricted | 111,799 | 1.3906 | 1.0570 | 0.7467 | 0.8667 | 0.8226 | 0.9040 | 0.8254 |
| GE | fingerprint | 111,799 | 1.3639 | 1.0324 | 0.7563 | 0.8722 | 0.8273 | 0.9075 | 0.8271 |
| Mut_CNV | onehot | 134,764 | 2.1662 | 1.6036 | 0.3582 | 0.6473 | 0.5995 | 0.7886 | 0.7079 |
| Mut_CNV | onehot_restricted | 111,799 | 2.2588 | 1.6761 | 0.3316 | 0.6341 | 0.5902 | 0.7833 | 0.7064 |
| Mut_CNV | fingerprint | 111,799 | 1.9904 | 1.4624 | 0.4811 | 0.7253 | 0.6676 | 0.8229 | 0.7646 |
| Proteomics | onehot | 134,764 | 1.2921 | 0.9822 | 0.7716 | 0.8787 | 0.8325 | 0.9093 | 0.8393 |
| Proteomics | onehot_restricted | 111,799 | 1.3305 | 1.0096 | 0.7681 | 0.8765 | 0.8320 | 0.9085 | 0.8332 |
| **Proteomics** | **fingerprint** | 111,799 | **1.2843** | **0.9721** | **0.7839** | **0.8866** | **0.8461** | **0.9153** | **0.8328** |
| GE+Mut_CNV | onehot | 134,764 | 1.8493 | 1.3855 | 0.5323 | 0.7404 | 0.7052 | 0.8411 | 0.7771 |
| GE+Mut_CNV | onehot_restricted | 111,799 | 2.1084 | 1.4724 | 0.4176 | 0.7026 | 0.7062 | 0.8403 | 0.7715 |
| GE+Mut_CNV | fingerprint | 111,799 | 1.6010 | 1.2245 | 0.6642 | 0.8280 | 0.7769 | 0.8772 | 0.8004 |
| GE+Proteomics | onehot | 134,764 | 1.5128 | 1.1487 | 0.6870 | 0.8309 | 0.7887 | 0.8852 | 0.8087 |
| GE+Proteomics | onehot_restricted | 111,799 | 1.5559 | 1.1851 | 0.6829 | 0.8302 | 0.7826 | 0.8824 | 0.8044 |
| GE+Proteomics | fingerprint | 111,799 | 1.4515 | 1.1021 | 0.7240 | 0.8575 | 0.8128 | 0.8986 | 0.8161 |
| Mut_CNV+Proteomics | onehot | 134,764 | 1.9207 | 1.4440 | 0.4954 | 0.7273 | 0.7008 | 0.8406 | 0.7432 |
| Mut_CNV+Proteomics | onehot_restricted | 111,799 | 2.0156 | 1.5015 | 0.4678 | 0.7071 | 0.6809 | 0.8281 | 0.7330 |
| Mut_CNV+Proteomics | fingerprint | 111,799 | 1.6241 | 1.2340 | 0.6545 | 0.8248 | 0.7728 | 0.8764 | 0.7854 |
| GE+Mut_CNV+Proteomics | onehot | 134,764 | 1.8202 | 1.3737 | 0.5469 | 0.7549 | 0.7300 | 0.8591 | 0.7668 |
| GE+Mut_CNV+Proteomics | onehot_restricted | 111,799 | 1.9372 | 1.4623 | 0.5084 | 0.7361 | 0.7120 | 0.8473 | 0.7465 |
| GE+Mut_CNV+Proteomics | fingerprint | 111,799 | 1.5795 | 1.2079 | 0.6732 | 0.8357 | 0.7890 | 0.8840 | 0.7882 |

Mean-only floor: **2.7097** (134,764 pairs) / **2.7690** (111,799 pairs).

## 1. The MLP-vs-RF control: a network beats trees by a wide margin on one-hot

Test RMSE at matched cells, `onehot` arm:

| Omics | RF | MLP | MLP advantage |
|---|---|---|---|
| GE | 1.9823 | **1.3995** | 0.583 |
| Mut_CNV | 2.1060 | 2.1662 | −0.060 |
| Proteomics | 2.0000 | **1.2921** | 0.708 |
| GE+Mut_CNV | 2.1060 | **1.8493** | 0.257 |
| GE+Proteomics | 2.1093 | **1.5128** | 0.597 |
| Mut_CNV+Proteomics | 2.1527 | **1.9207** | 0.232 |
| tri-omics | 2.1721 | **1.8202** | 0.352 |

**This substantially reframes the earlier cross-attention result.**
[`cross_attention_fusion_results.md`](./cross_attention_fusion_results.md)
attributed a ~42% RMSE drop (2.1933 → 1.2658) to cross-attention fusion, while
noting it couldn't separate "learned fusion" from "gradient descent vs. trees."
The matrix now separates them:

```
RF tri-omics onehot            2.1721
MLP tri-omics onehot           1.8202   <- 0.352 from switching learner
CrossAttn tri-omics onehot     1.2658   <- 0.554 more from learned fusion
```

So of the total 0.906 improvement, roughly **39% comes from using a neural net
at all** and **61% from cross-attention fusion specifically**. Fusion is doing
real work — but the headline "42% better than RF" was crediting fusion with the
learner's share too.

## 2. Fingerprints help far less for a network than for a tree

| Omics | RF gain (restricted→fp) | MLP gain (restricted→fp) |
|---|---|---|
| GE | 0.549 | 0.027 |
| Proteomics | 0.547 | 0.046 |
| tri-omics | 0.787 | 0.358 |

RF gains 0.55–0.79; MLP gains 0.03–0.51, and on its best cells (GE,
Proteomics) the gain is nearly nil. A network can already learn a dense
embedding of drug identity from one-hot input via its first weight matrix — it
doesn't need structural features handed to it the way a `max_features="sqrt"`
tree does (see [`rf_ablation_results.md`](./rf_ablation_results.md) §1).

**Implication for the proposal:** the argument for Morgan fingerprints in the
final architecture cannot rest on in-distribution accuracy — for neural models
it barely moves. It rests on *generalization to unseen drugs*, which one-hot
structurally cannot do and which this matrix does not test (the split is
grouped by cell line, not by drug).

That experiment has since been run — see
[`leave_drugs_out_results.md`](./leave_drugs_out_results.md). Holding out drugs
instead of cell lines, MLP one-hot drops to **R² 0.075** (barely above the
mean-only floor) while fingerprints hold at **R² 0.463**, a +0.578 RMSE gap.
The fingerprint's value is real; it simply is not visible under a
cell-line-grouped split.

## 3. The population control changes sign versus RF

Unlike RF (where the restricted population was uniformly *easier*), the MLP
finds it **harder** in 6 of 7 subsets — up to 0.259 RMSE worse on GE+Mut_CNV.
Had we compared raw `onehot` against `fingerprint` without the control, the
MLP's fingerprint benefit would have been *overstated* on every subset. This is
the concrete payoff of the third arm.

## 4. "More omics hurts" reproduces — and is therefore not a tree artifact

```
Proteomics alone     1.2921   <- best
GE alone             1.3995
GE+Proteomics        1.5128
tri-omics            1.8202
Mut_CNV+Proteomics   1.9207
Mut_CNV alone        2.1662   <- worst
```

Same monotonic degradation as RF. Since MLP has no `max_features` parameter,
the `sqrt`-dilution explanation offered in
[`rf_ablation_results.md`](./rf_ablation_results.md) §2 **cannot be the whole
story**. Two candidate explanations remain, untested:

1. **Signal-to-noise.** Mut_CNV is weak on its own (2.1662) and drags down
   every combination it enters. Concatenation gives every modality equal
   width, so a noisy modality contributes proportional noise regardless of its
   informativeness — exactly the deficiency learned fusion is supposed to fix.
   Cross-attention's results support this reading: it recovers most of the
   loss (tri-omics 1.2658 vs MLP's 1.8202).
2. **Fixed capacity spread thinner.** Hidden dims are constant at
   `[512,256,128]` regardless of input width, so tri-omics asks the same
   parameter budget to model 3× the input.

Distinguishing these would need a capacity-matched sweep (scale hidden dims
with input dim). Not run — flagged as a follow-up.

## Caveats

- Single-modality cells still use only the 532 tri-omics-complete cell lines,
  so they are not "what you'd get training on all cell lines with proteomics."
- `mlp_baseline.py`'s original run was never logged; this matrix supersedes it
  and additionally pins its previously-unset `RANDOM_STATE`.
