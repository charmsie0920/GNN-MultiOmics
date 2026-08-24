# Cross-Attention Fusion Ablation — Omics Subsets × Drug Representations

**Script:** [`experiments/Full Matrix/cross_attention_matrix.py`](../experiments/Full%20Matrix/cross_attention_matrix.py)
**Raw results:** `experiments/Full Matrix/cross_attention_matrix_results.csv`
**Companions:** [`rf_ablation_results.md`](./rf_ablation_results.md), [`mlp_ablation_results.md`](./mlp_ablation_results.md), [`gnn_ablation_results.md`](./gnn_ablation_results.md), [`results.md`](./results.md)
**Date:** 2026-08-24 · **Runtime:** 12 runs in 1466.0s

## What this covers

The 4 omics subsets with ≥2 modalities × 3 drug-representation arms = 12 cells.
Single-modality cells are excluded by construction: cross-attention needs at
least two modalities to attend *between*, and that case is already covered by
[`mlp_ablation_results.md`](./mlp_ablation_results.md).

This extends the single tri-omics point in
[`cross_attention_fusion_results.md`](./cross_attention_fusion_results.md)
(RMSE 1.2658), which reproduces here exactly.

**Module change required:** `MultiOmicsCrossAttentionFusion`
([`src/models/cross_attention_fusion.py`](../src/models/cross_attention_fusion.py))
was hardcoded to all 3 modalities and raised `KeyError` on a subset. It now
takes a `modalities` argument and derives its attention pairs and projection
width from it — a 2-modality fusion builds 2 pairwise blocks instead of 6.
Existing 3-modality callers are unaffected (default unchanged).

For fingerprint-mode cells the 2048-bit block passes through a
`FingerprintEncoder` (`Linear(2048→128) → ReLU → Dropout → Linear(128→128)`)
before concatenation with the 256-dim fused embedding — concatenating raw
sparse bits would let the drug block dominate the head's input by width alone.

## Results

| Omics | Drug rep | Pairs | RMSE | MAE | R² | PCC | SCC | AUC | F1 |
|---|---|---|---|---|---|---|---|---|---|
| GE+Mut_CNV | onehot | 134,764 | 1.2961 | 0.9701 | 0.7702 | 0.8810 | 0.8376 | 0.9132 | 0.8368 |
| GE+Mut_CNV | onehot_restricted | 111,799 | 1.3667 | 1.0153 | 0.7553 | 0.8713 | 0.8248 | 0.9057 | 0.8303 |
| GE+Mut_CNV | fingerprint | 111,799 | 1.3949 | 1.0466 | 0.7451 | 0.8668 | 0.8217 | 0.9031 | 0.8178 |
| **GE+Proteomics** | **onehot** | 134,764 | **1.2442** | **0.9336** | **0.7883** | **0.8897** | **0.8480** | **0.9191** | **0.8423** |
| GE+Proteomics | onehot_restricted | 111,799 | 1.2847 | 0.9606 | 0.7838 | 0.8861 | 0.8430 | 0.9157 | 0.8386 |
| GE+Proteomics | fingerprint | 111,799 | 1.3205 | 0.9891 | 0.7716 | 0.8837 | 0.8424 | 0.9145 | 0.8305 |
| Mut_CNV+Proteomics | onehot | 134,764 | 1.2814 | 0.9482 | 0.7754 | 0.8827 | 0.8393 | 0.9134 | 0.8447 |
| Mut_CNV+Proteomics | onehot_restricted | 111,799 | 1.2987 | 0.9693 | 0.7790 | 0.8840 | 0.8413 | 0.9148 | 0.8394 |
| Mut_CNV+Proteomics | fingerprint | 111,799 | 1.3787 | 1.0367 | 0.7510 | 0.8725 | 0.8296 | 0.9070 | 0.8171 |
| GE+Mut_CNV+Proteomics | onehot | 134,764 | 1.2658 | 0.9383 | 0.7809 | 0.8849 | 0.8428 | 0.9160 | 0.8432 |
| GE+Mut_CNV+Proteomics | onehot_restricted | 111,799 | 1.2992 | 0.9716 | 0.7789 | 0.8854 | 0.8426 | 0.9141 | 0.8384 |
| GE+Mut_CNV+Proteomics | fingerprint | 111,799 | 1.3589 | 1.0262 | 0.7581 | 0.8752 | 0.8297 | 0.9064 | 0.8165 |

Mean-only floor: **2.7097** (134,764 pairs) / **2.7690** (111,799 pairs).

## 1. Best result in the entire matrix — and it drops a modality

**GE+Proteomics with one-hot drugs: RMSE 1.2442 / PCC 0.8897 / R² 0.7883.**

This beats the previously published tri-omics configuration (1.2658) by
excluding Mut_CNV entirely. Across all three arms, GE+Proteomics either wins or
ties the tri-omics cell:

| Omics | onehot | onehot_restricted | fingerprint |
|---|---|---|---|
| **GE+Proteomics** | **1.2442** | **1.2847** | **1.3205** |
| tri-omics | 1.2658 | 1.2992 | 1.3589 |
| Mut_CNV+Proteomics | 1.2814 | 1.2987 | 1.3787 |
| GE+Mut_CNV | 1.2961 | 1.3667 | 1.3949 |

Consistent with [`rf_ablation_results.md`](./rf_ablation_results.md) §3 and
[`mlp_ablation_results.md`](./mlp_ablation_results.md) §4 identifying Mut_CNV as
the weakest modality — but here the effect is small (0.022 RMSE), not the
0.3–0.4 degradation flat concatenation suffers. **Learned fusion largely
absorbs a weak modality instead of being dragged down by it**, which is exactly
the behavior the architecture was proposed to deliver.

## 2. Cross-attention is far more robust to modality count than flat concatenation

Spread between best and worst omics subset, `onehot` arm:

| Model | Best | Worst | Spread |
|---|---|---|---|
| RF | 1.9823 | 2.1721 | 0.190 |
| MLP | 1.2921 | 2.1662 | 0.874 |
| **Cross-attention** | **1.2442** | **1.2961** | **0.052** |

Every cross-attention cell lands within 0.052 RMSE. The choice of modality
subset barely matters once fusion is learned — whereas for MLP it swings
performance by 0.874. This is the strongest evidence in the matrix that the
fusion module is doing something real rather than acting as extra capacity.

## 3. Fingerprints *hurt* cross-attention

Uniquely among the model families, the fingerprint arm is **worst** in all four
subsets (comparing against `onehot_restricted`, the matched-population
control):

| Omics | onehot_restricted | fingerprint | Δ |
|---|---|---|---|
| GE+Proteomics | 1.2847 | 1.3205 | −0.036 |
| tri-omics | 1.2992 | 1.3589 | −0.060 |
| GE+Mut_CNV | 1.3667 | 1.3949 | −0.028 |
| Mut_CNV+Proteomics | 1.2987 | 1.3787 | −0.080 |

Reading across all three model families, the pattern is monotone in model
capacity: RF gains 0.55–0.79 from fingerprints, MLP gains 0.03–0.51,
cross-attention *loses* 0.03–0.08. The more expressive the model, the less
structural priors help — and at the top end the 2048→128 encoder appears to be
a lossy bottleneck compared to a one-hot lookup the model can embed exactly.

**This does not mean fingerprints are the wrong choice for the final
architecture.** One-hot cannot represent a drug outside the 295-drug training
vocabulary at all; fingerprints can. What these numbers show is that the
fingerprint's value is *not* in-distribution accuracy, so the case for it must
be made on unseen-drug generalization — which requires a leave-drugs-out split
this matrix does not run. See
[`mlp_ablation_results.md`](./mlp_ablation_results.md) §2.

## 4. Reproduction check

The tri-omics/onehot cell reproduces
[`cross_attention_fusion_results.md`](./cross_attention_fusion_results.md)
exactly: RMSE **1.2658**, PCC **0.8849**, AUC **0.9160**, F1 **0.8432** —
confirming the fusion-module generalization did not alter 3-modality behavior.

## Caveats

- **Attention weights remain degenerate.** As documented in
  [`preprocessing_and_fusion_module.md`](./preprocessing_and_fusion_module.md),
  each modality is a single pooled vector, so every attention call has one
  query and one key and softmax is 1.0 by construction. The gains come from the
  learned Q/K/V projections and residual FFN stack, not from a meaningful
  attention distribution. Attention weights here are not interpretable and
  should not be presented as such.
- Cross-attention is the most expensive family in the matrix (~2 min/run vs
  ~23s for RF), driven by 6 pairwise attention blocks at 3 modalities.
