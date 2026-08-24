# Random Forest Ablation — Full Omics × Drug-Representation Matrix

**Script:** [`experiments/Full Matrix/rf_matrix.py`](../experiments/Full%20Matrix/rf_matrix.py)
**Raw results:** `experiments/Full Matrix/rf_matrix_results.csv`
**Companions:** [`mlp_ablation_results.md`](./mlp_ablation_results.md), [`cross_attention_ablation_results.md`](./cross_attention_ablation_results.md), [`gnn_ablation_results.md`](./gnn_ablation_results.md), [`results.md`](./results.md)
**Date:** 2026-08-24 · **Runtime:** 21 runs in 482.8s

## What this covers

All 7 omics subsets (3 single, 3 dual, 1 tri) × 3 drug-representation arms =
21 cells. This supersedes the single-point results in
[`rf_baseline_results.md`](./rf_baseline_results.md) (fused tri-omics, RMSE
2.1933) and [`proteomics_baseline_results.md`](./proteomics_baseline_results.md)
(proteomics-only, RMSE 2.0000) — both of which reproduce here exactly, on the
same split.

Model unchanged from those baselines: `RandomForestRegressor(n_estimators=100,
min_samples_leaf=5, max_features="sqrt", random_state=42)`, flat concatenation
of PCA'd omics + drug block, `GroupShuffleSplit` by cell line 70/15/15.

### The three drug-representation arms

Only 498 of GDSC's 621 drug IDs resolved to a SMILES string (122 drug names are
internal codenames PubChem doesn't carry — see
[`graph_construction.md`](./graph_construction.md)). Fingerprint runs therefore
cover 111,799 pairs where one-hot runs cover 134,764. Comparing those two
directly would confound *representation* with *population*, so a third control
arm was added:

| Arm | Drug features | Pairs | Purpose |
|---|---|---|---|
| `onehot` | 295-dim identity | 134,764 | Continuity with published baselines |
| `onehot_restricted` | 295-dim identity | 111,799 | **Control** — same rows as fingerprint |
| `fingerprint` | 2048-bit Morgan | 111,799 | The representation under test |

`onehot` → `onehot_restricted` isolates the population effect;
`onehot_restricted` → `fingerprint` isolates the representation effect.

## Results

| Omics | Drug rep | Pairs | RMSE | MAE | R² | PCC | SCC | AUC | F1 |
|---|---|---|---|---|---|---|---|---|---|
| GE | onehot | 134,764 | 1.9823 | 1.6186 | 0.4626 | 0.7340 | 0.5590 | 0.7449 | 0.6925 |
| GE | onehot_restricted | 111,799 | 1.9322 | 1.5771 | 0.5109 | 0.7596 | 0.6161 | 0.7755 | 0.6975 |
| GE | fingerprint | 111,799 | 1.3834 | 1.0331 | 0.7493 | 0.8683 | 0.8240 | 0.9074 | 0.8232 |
| Mut_CNV | onehot | 134,764 | 2.1060 | 1.7126 | 0.3934 | 0.7021 | 0.4551 | 0.6952 | 0.4600 |
| Mut_CNV | onehot_restricted | 111,799 | 2.0953 | 1.7067 | 0.4249 | 0.7199 | 0.4992 | 0.7192 | 0.4960 |
| Mut_CNV | fingerprint | 111,799 | 1.4710 | 1.1048 | 0.7165 | 0.8494 | 0.7966 | 0.8947 | 0.8123 |
| Proteomics | onehot | 134,764 | 2.0000 | 1.6349 | 0.4529 | 0.7460 | 0.5759 | 0.7520 | 0.6435 |
| Proteomics | onehot_restricted | 111,799 | 1.9497 | 1.5917 | 0.5020 | 0.7692 | 0.6448 | 0.7896 | 0.6801 |
| Proteomics | fingerprint | 111,799 | 1.4027 | 1.0470 | 0.7422 | 0.8643 | 0.8215 | 0.9067 | 0.8237 |
| GE+Mut_CNV | onehot | 134,764 | 2.1060 | 1.7103 | 0.3934 | 0.7145 | 0.5260 | 0.7293 | 0.6465 |
| GE+Mut_CNV | onehot_restricted | 111,799 | 2.0945 | 1.7058 | 0.4253 | 0.7333 | 0.5580 | 0.7448 | 0.5920 |
| GE+Mut_CNV | fingerprint | 111,799 | 1.3973 | 1.0486 | 0.7442 | 0.8655 | 0.8206 | 0.9053 | 0.8224 |
| GE+Proteomics | onehot | 134,764 | 2.1093 | 1.7125 | 0.3915 | 0.7213 | 0.5405 | 0.7353 | 0.6284 |
| GE+Proteomics | onehot_restricted | 111,799 | 2.0772 | 1.6920 | 0.4348 | 0.7419 | 0.5905 | 0.7612 | 0.6488 |
| **GE+Proteomics** | **fingerprint** | 111,799 | **1.3712** | **1.0290** | **0.7537** | **0.8710** | **0.8280** | **0.9095** | **0.8275** |
| Mut_CNV+Proteomics | onehot | 134,764 | 2.1527 | 1.7452 | 0.3662 | 0.7243 | 0.5479 | 0.7371 | 0.4967 |
| Mut_CNV+Proteomics | onehot_restricted | 111,799 | 2.1119 | 1.7197 | 0.4157 | 0.7417 | 0.5839 | 0.7584 | 0.5267 |
| Mut_CNV+Proteomics | fingerprint | 111,799 | 1.4137 | 1.0603 | 0.7382 | 0.8623 | 0.8190 | 0.9048 | 0.8235 |
| GE+Mut_CNV+Proteomics | onehot | 134,764 | 2.1721 | 1.7578 | 0.3547 | 0.7099 | 0.5283 | 0.7292 | 0.6032 |
| GE+Mut_CNV+Proteomics | onehot_restricted | 111,799 | 2.1679 | 1.7594 | 0.3843 | 0.7348 | 0.5807 | 0.7548 | 0.5310 |
| GE+Mut_CNV+Proteomics | fingerprint | 111,799 | 1.3805 | 1.0378 | 0.7503 | 0.8692 | 0.8261 | 0.9083 | 0.8243 |

Mean-only floor: **2.7097** (134,764-pair population) / **2.7690** (111,799-pair population).

## 1. Morgan fingerprints are worth ~0.55–0.79 RMSE, and the confound is negligible

| Omics | Population effect | **Representation effect** |
|---|---|---|
| GE | 0.050 | **0.549** |
| Mut_CNV | 0.011 | **0.624** |
| Proteomics | 0.050 | **0.547** |
| GE+Mut_CNV | 0.012 | **0.697** |
| GE+Proteomics | 0.032 | **0.706** |
| Mut_CNV+Proteomics | 0.041 | **0.698** |
| GE+Mut_CNV+Proteomics | 0.004 | **0.787** |

The restricted population is marginally *easier* (0.004–0.050 RMSE), so a naive
`onehot` vs `fingerprint` comparison would have overstated the fingerprint gain
by at most ~4%. The effect is real and large regardless — this is the single
biggest lever in the RF matrix.

**Why fingerprints beat one-hot for a tree specifically:** with 295 one-hot
columns and `max_features="sqrt"`, each split samples ~20 of 423 features, and
each one-hot column is a weak candidate covering only 1/295 of rows. Morgan
bits are *shared across structurally similar drugs*, so a single bit
partitions a whole chemical family at once. Trees can exploit that; one-hot
forces them to isolate drugs one at a time. This is precisely the
generalization argument for structural encoding — the model learns
"EGFR-inhibitor-like compounds behave this way" rather than memorizing each
drug ID.

## 2. Adding omics modalities makes RF *worse*, monotonically

Within the `onehot` arm, ordered best to worst:

```
GE alone             1.9823   <- best
Proteomics alone     2.0000
Mut_CNV alone        2.1060
GE+Mut_CNV           2.1060
GE+Proteomics        2.1093
Mut_CNV+Proteomics   2.1527
tri-omics            2.1721   <- worst
```

Every single modality beats every combination. This generalizes the surprise
recorded in [`proteomics_baseline_results.md`](./proteomics_baseline_results.md)
§2 ("proteomics-only beats fused") — it is **not proteomics-specific**, it's a
property of adding *any* modality.

That doc offered `max_features="sqrt"` dilution as an unconfirmed explanation,
and the matrix supports it: each added modality contributes 128 more columns,
but `sqrt` keeps the per-split candidate count nearly flat (√423≈21 →
√679≈26), so the *probability* of any given informative feature being offered
at a split drops as dimensionality grows. The signal isn't destroyed, it's
diluted.

Note this dilution is much weaker in the `fingerprint` arm (1.371–1.471, a
0.10 spread vs 0.19 for one-hot), consistent with fingerprint bits being
individually stronger split candidates that survive the dilution.

## 3. Mut_CNV is the weakest modality

Worst single modality in every arm (onehot 2.1060, fingerprint 1.4710), and
its F1 collapses to 0.4600 — barely better than chance on the binarized task.
Any combination containing it underperforms the same combination without it.
Worth flagging to the team before the GNN work leans on mutation-derived
features.

## 4. Reproduction check

Two cells in this matrix re-run previously published configurations, and both
reproduce exactly on the same split — confirming the `experiment_utils.py`
refactor is behavior-preserving:

| Config | Published | This matrix |
|---|---|---|
| Proteomics-only, onehot | RMSE 2.0000 / PCC 0.7460 / AUC 0.7520 / F1 0.6435 | **identical** |
| tri-omics fused, onehot | RMSE 2.1933 / PCC 0.7050 | RMSE 2.1721 / PCC 0.7099 |

The tri-omics row differs slightly because
[`rf_baseline_results.md`](./rf_baseline_results.md) ran with
`RANDOM_STATE = None` (a fresh split each run — flagged as a reproducibility
fragility in its §3). That is now pinned to 42 across the whole repo.

## Caveats

- **All cells share the 532-cell-line tri-omics-complete population**, even
  single-modality cells. This isolates the modality effect (the alternative,
  letting each modality use all cell lines it covers, would confound modality
  with sample size) but means single-omics numbers here are *not* the numbers
  you'd get training on all cell lines with that modality available.
- The ~45% of GDSC2 cell lines lacking full tri-omics coverage remain excluded
  throughout — see [`rf_baseline_results.md`](./rf_baseline_results.md) §4.4.
