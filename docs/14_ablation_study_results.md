# Ablation Study — Base + One Component at a Time

**Script:** [`experiments/14_benchmark_alignment/benchmark_alignment.py`](../experiments/14_benchmark_alignment/benchmark_alignment.py)
**Model:** [`src/models/mographdrp_aligned.py`](../src/models/mographdrp_aligned.py)
**Raw results:** `experiments/14_benchmark_alignment/benchmark_alignment_results.csv`
**Companions:** [`results.md`](./results.md), [`13_seed_variance_results.md`](./13_seed_variance_results.md), [`09_split_protocol_comparison.md`](./09_split_protocol_comparison.md)
**Date:** TBD · **Runtime:** TBD · **Hardware:** TBD
**Status:** DRAFT TEMPLATE. Every `TBD` is to be filled from the results CSV; no number in this file is a result yet.

> How to use this file: fill the tables from the CSV, then write the short
> prose under each heading. Anything marked *(report)* is meant to be lifted
> into the final report. Delete these notes when done.

---

## 1. What this study answers

*(report)* One or two sentences. Suggested:

> Starting from the simplest pipeline that functions, which individual
> components improve prediction of ln(IC50) on unseen cell lines, by how much,
> and does combining them give the best model?

Requirement from the lecturer: **Base + a**, **Base + b**, **Base + c**, ...,
where each letter is one component, and Base + everything is the final model.

## 2. What is held fixed across every row

These must be identical for every rung, otherwise a delta is not attributable.

| Item | Setting |
|---|---|
| Pairs | 111,799 (cell line, drug) pairs; 532 cell lines; 240 drug IDs |
| Split | grouped by cell line, 70/15/15, `random_state=42` (77,668 / 17,291 / 16,840) |
| Omics preprocessing | per-modality scaling + PCA to 128 |
| Optimiser | Adam, lr 1e-4, weight decay 1e-5 |
| Batch size / dropout | 128 / 0.4 |
| Epochs / early stopping | max 200, patience 15 on validation RMSE |
| LR schedule | ReduceLROnPlateau, factor 0.5, patience 5 |
| Checkpoint | best validation RMSE, in ln(IC50) units |
| Seeds | TBD (list them) |
| Metrics | RMSE, MAE, R², PCC, SCC; AUC and F1 at the shared median threshold |

Reference floors on the same test fold:

| Floor | Test RMSE |
|---|---|
| Global training mean | 2.7690 |
| **Per-drug training mean** (no omics, no model) | **1.4889** |

**Noise band:** TBD from this study's own seeds (std of `base`). Until then use
~0.03 RMSE from [13_seed_variance](./13_seed_variance_results.md). A delta
smaller than the band is reported as "no measurable effect".

## 3. The base

*(report)* Describe it as a pipeline, then justify each choice as "simplest".

| Stage | Base choice | Why this is the simplest functioning option |
|---|---|---|
| Omics inputs | GE + Mut_CNV | TBD (two inputs are the minimum for fusion to apply; closest to the benchmark's inputs) |
| Cell encoder | one 2-layer branch per omics, concatenated | TBD |
| Drug encoder | Morgan fingerprint (radius 2, 2048 bits) → MLP | TBD |
| Cell–drug interaction | concatenation | TBD |
| Predictor | MLP 256 → 128 → 1 | TBD |
| Targets | raw ln(IC50) | TBD |

Parameters: 544,001. Diagram: TBD (`docs/figures/`).

## 4. The components

One row per component. "Switch" is the single field it changes in `BASE`.

| ID | Component | Switch | Replaces / adds | Origin | Hypothesis (why it should help) |
|---|---|---|---|---|---|
| a | Cross-attention fusion | `fusion=cross_attention` | replaces concat of omics branches | this project | TBD |
| b | Molecular-graph GCN | `drug_mode=both` | adds a GCN over the drug's atoms beside the fingerprint | MoGraphDRP §2.2.3 | TBD |
| c | Bilinear head | `head=bilinear` | replaces concat of cell and drug vectors | MoGraphDRP §2.3 | TBD |
| d | Proteomics | `modalities=+Proteomics` | adds a third omics branch | this project | TBD |
| e | Per-drug target standardisation | `standardize_targets=True` | trains on per-drug z-scores, inverted before scoring | this project | TBD |
| f | PPI GNN + pair-specific attention | TBD | TBD | this project (novelty component) | TBD |

Known caveats to state honestly in the report:

- **a:** each modality is one vector, so attention weights are always 1.0. As
  built, this rung measures the extra projection layers, not attention. TBD:
  fixed with multi-token modalities, or reported with this caveat.
- **e:** not usable under leave-drugs-out (unseen drugs have no statistics).
- **f:** not built yet. Covers only pairs whose drug has a known target
  (TBD % of pairs); 130 drugs have none.

## 5. Main result — grouped split

*(report)* This is the central table. Mean ± std over TBD seeds, test fold.

| Config | Components | RMSE | Δ vs base | Beyond noise? | PCC | R² | Gain over per-drug mean | Params | Fit (s) |
|---|---|---|---|---|---|---|---|---|---|
| per-drug mean | — | 1.4889 | — | — | — | — | 0 | 0 | — |
| `base` | none | TBD ± TBD | 0 | — | TBD | TBD | TBD | 544,001 | TBD |
| `base+cross_attention` | a | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| `base+mol_graph` | b | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| `base+bilinear` | c | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| `base+proteomics` | d | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| `base+std_targets` | e | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| `base+pair` | f | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| `aligned` | b + c | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| `full` | all | TBD | TBD | TBD | TBD | TBD | TBD | 1,702,913 | TBD |

Negative Δ means better than base.

**Is each delta real?** Fill one line per component:

| Component | Δ RMSE | Seeds where it beat base | Test used (TBD: paired bootstrap over test cell lines, as in [16](./16_target_edge_ablation_results.md)) | 95% CI | Verdict |
|---|---|---|---|---|---|
| a | TBD | TBD / TBD | TBD | TBD | helps / no effect / hurts |
| b | TBD | | | | |
| c | TBD | | | | |
| d | TBD | | | | |
| e | TBD | | | | |
| f | TBD | | | | |

## 6. Per-component findings

*(report)* One short block per component. Keep to this shape so they read
consistently.

### a. Cross-attention fusion
- **Result:** TBD (Δ, verdict).
- **Interpretation:** TBD. What does the result say about the hypothesis?
- **Why (evidence, not speculation):** TBD. Point to a measurement.
- **Cost:** TBD parameters, TBD × training time.

### b. Molecular-graph GCN
- **Result / Interpretation / Why / Cost:** TBD
- Compare with [11_molecular_graph_results](./11_molecular_graph_results.md), where the single-run ranking reversed across seeds.

### c. Bilinear head
- **Result / Interpretation / Why / Cost:** TBD
- State whether the head matches the MoGraphDRP repo's implementation (TBD: checked on date).

### d. Proteomics
- **Result / Interpretation / Why / Cost:** TBD
- Earlier flat-model result to reconcile: Mut_CNV hurt (E32 vs E04 in [results.md](./results.md)).

### e. Per-drug target standardisation
- **Result / Interpretation / Why / Cost:** TBD
- Earlier single-seed result on HeteroIC50GNN: 1.3449 → 1.2917.

### f. PPI GNN + pair-specific attention
- **Result / Interpretation / Why / Cost:** TBD
- **GNN depth ablation** (does message passing matter?):

| Message-passing layers | RMSE (all pairs) | RMSE (pairs with a known target) |
|---|---|---|
| 0 (attention on raw embeddings, no GNN) | TBD | TBD |
| 1 | TBD | TBD |
| 2 | TBD | TBD |

- **Coverage:** TBD % of test pairs have a target; report both subsets.

## 7. Do the components combine?

*(report)* Components can overlap or interfere, so the sum of single gains
rarely equals the gain of the full model.

| Quantity | RMSE change vs base |
|---|---|
| Sum of the individual deltas (a + b + ... ) | TBD |
| `full` measured | TBD |
| Difference (negative = components overlap) | TBD |

- Does `full` beat the best single component beyond noise? TBD
- Is any component harmful inside `full`? Leave-one-out check (`full` minus one), if time allows:

| Config | RMSE | Δ vs full |
|---|---|---|
| `full` | TBD | 0 |
| `full − a` | TBD | TBD |
| `full − b` | TBD | TBD |
| ... | | |

## 8. Final model

*(report)*

- **Chosen configuration:** TBD (list the components kept).
- **Components dropped and why:** TBD (each with its measured delta).
- **Final test RMSE / PCC / R²:** TBD ± TBD over TBD seeds.
- **Gain over the per-drug mean floor:** TBD.

Selection must be made on **validation** RMSE; the test numbers above are
reported after the choice, not used to make it. TBD: confirm this was done.

## 9. Sanity check against the benchmark — random split

Not a ranking. It checks that the reproduction (`aligned`) is faithful.

| Config | Random-split RMSE | Reference |
|---|---|---|
| `aligned` | TBD | MoGraphDRP without XGBoost: 0.9497 |
| `base` | TBD | — |
| `full` | TBD | MoGraphDRP with XGBoost: 0.6622 |

If `aligned` is far from 0.95, list the likely causes: 2–3 omics vs their 4,
PCA vs COSMIC gene filtering, one fingerprint type vs three, head
implementation. TBD.

Protocol gap for the final model (random − grouped): TBD. Compare with the
0.347 measured in [09](./09_split_protocol_comparison.md).

## 10. Comparison with the other models on the same pairs

Same 111,799 pairs, same grouped split. TBD seeds each.

| Model | RMSE | PCC | Notes |
|---|---|---|---|
| Per-drug mean | 1.4889 | — | floor |
| RF | TBD | TBD | |
| MLP | TBD | TBD | |
| CrossAttention (E10) | 1.3029 ± 0.0126 | TBD | 5 seeds, [13](./13_seed_variance_results.md) |
| HeteroIC50GNN (E11) | TBD | TBD | single run 1.3301 |
| `aligned` (MoGraphDRP-style) | TBD | TBD | |
| **Final model** | TBD | TBD | |

## 11. Further checks (fill if run)

- **Leave-drugs-out** (unseen compounds, no standardisation): TBD vs the MLP's 1.8516 in [10](./10_leave_drugs_out_results.md).
- **Interpretability:** does attention on BRAF track measured sensitivity against untrained models? TBD, see [15](./15_interpretability_validation_results.md).

## 12. Limitations to state

Tick the ones that apply and add numbers.

- [ ] Seeds: TBD per row. The split itself is fixed at seed 42, so variance from the choice of test cell lines is not measured.
- [ ] No hyperparameter tuning per rung; all use the benchmark's Table 1 settings.
- [ ] Base uses PCA features, not gene-level features; differs from the benchmark's COSMIC filtering.
- [ ] Only 532 of GDSC2's 969 cell lines (tri-omics-complete subset).
- [ ] Component a is degenerate as built (see §4).
- [ ] Component f covers TBD % of pairs.
- [ ] The benchmark comparison is against a reimplementation, not the authors' code on their data.

## 13. Figures to produce

| Figure | Content | Status |
|---|---|---|
| Ladder bar chart | Δ RMSE vs base per component, error bars = std over seeds, noise band shaded | TBD |
| Architecture diagram | base in grey, each component highlighted | TBD |
| Additivity plot | sum of single gains vs `full` | TBD |
| Protocol comparison | grouped vs random for base / aligned / final | TBD |

## 14. Report wording (draft once results are in)

*(report)* Three to five sentences, each backed by a table above. Template:

> The base model reaches RMSE TBD on unseen cell lines, against TBD for a
> per-drug mean. Of the TBD components tested individually, TBD improved on the
> base beyond seed-to-seed variation (TBD), and TBD did not. Combining TBD
> gives the final model at RMSE TBD, a gain of TBD over the base and TBD over
> the per-drug mean.

## 15. Reproduction

```
python "experiments/14_benchmark_alignment/benchmark_alignment.py" --list
python "experiments/14_benchmark_alignment/benchmark_alignment.py" --protocols grouped --seeds 42 43 44
python "experiments/14_benchmark_alignment/benchmark_alignment.py" --configs aligned base full --protocols random
```

Commit hash of the code that produced the CSV: TBD.
