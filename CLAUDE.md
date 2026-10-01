# CLAUDE.md

Drug-response prediction (ln IC50) on GDSC2 from multi-omics + drug structure.
FIT3161 final-year project. This file is the current plan; update it when a
phase finishes or a decision changes.

## Direction (decided 2026-10-01)

- **Base:** the MoGraphDRP-style model in `src/models/mographdrp_aligned.py`.
  It is the reproduced benchmark and the bottom of the ablation ladder.
- **Our model:** a new model built on that base, differing by measured
  components (below). Working name `PairGraphDRP` — placeholder, rename freely.
- **Novelty component:** pair-specific cross-attention between a drug's target
  proteins and a cell line's driver-mutated proteins over the STRING PPI graph.
- **`HeteroIC50GNN` is no longer the final model.** It stays in the results
  table as a compared graph baseline. Its own diagnostics show why
  (`experiments/15_diagnostics/`): 95.8% of proteins unattached, 130 drugs with
  no target edge, and no pair-specific path between mutation and target.

### What makes the final model ours, not MoGraphDRP

| Component | MoGraphDRP | Ours |
|---|---|---|
| Omics | expression, mutation, methylation, pathway | GE, Mut+CNV, **proteomics** |
| Biological prior | none on the cell side | **STRING PPI + drug-target + driver-mutation edges** |
| Cell-drug interaction | bilinear attention on two pooled vectors | bilinear **plus pair-specific target x mutation attention** |
| Target scaling | raw ln(IC50) | **per-drug standardisation** (cell-line split only) |
| Evaluation | random pair split | **cell-line-grouped split**, seeds, per-drug-mean floor |
| Interpretability | none validated | attention/attribution tested against **untrained-model nulls** |

Each row must be a switch and a row in the ablation table. A component that
does not beat the noise band is reported as tested-and-rejected, not kept to
make the model look different. MoGraphDRP is cited as the base architecture.

## Rules for every experiment

- **Headline protocol:** `grouped_split` (by cell line, 70/15/15, seed 42) from
  `src/data/experiment_utils.py`. Random split is a sanity check only.
- **Never compare against MoGraphDRP's published 0.6622 / 0.9497 as a ranking.**
  Different pairs, different split. Compare against the aligned model run here.
- **3-5 seeds per reported row**, mean +/- std. Differences under ~0.03 RMSE
  are noise (`docs/13_seed_variance_results.md`).
- **Always report the per-drug-mean floor (1.4889)** beside RMSE.
- **Checkpoint on best validation RMSE**, never Pearson r.
- **One change per rung.** New components are flags on one model, not forks.
- **Do not change the defaults of `MoGraphDRPAligned`** once Phase 0 passes; it
  is the frozen baseline.
- **`docs/results.md` is generated.** Edit `experiments/build_results_table.py`
  and re-run; never hand-edit.
- Per-drug standardisation cannot be used under leave-drugs-out (unseen drugs
  have no statistics).
- Hardware is a GTX 1650 (4 GB). Keep the graph full-batch and small; batch-128
  runs are slow, so launch long sweeps in the background and write results
  incrementally.

## Tasks

### Phase 0 — Validate the base
- [ ] Run `python "experiments/14_benchmark_alignment/benchmark_alignment.py"`
      (never run so far; no results CSV exists).
- [ ] Compare `MultiHeadBilinearAttention`, split code and hyperparameters
      against the MoGraphDRP repo (not in this workspace; ask for its path).
      The bilinear head was written from the paper text and is the likeliest
      mismatch. Check their seed too.
- [ ] Read MIDI (bioRxiv 2025.03.31.646490) in full: if its cross-attention
      already uses target proteins on the drug side, the novelty claim must be
      narrowed to the evaluation.

**Done when:** `aligned` under the random split is near 0.95 (their
no-XGBoost figure), or the gap is explained by the documented data differences.

### Phase 1 — Additive ablation ladder (lecturer requirement)
The lecturer requires **Base + a, Base + b, Base + c, ...**, where Base is the
simplest pipeline that functions with none of the improvements, and Base plus
everything is the final model. This is built in
`experiments/14_benchmark_alignment/benchmark_alignment.py`:

- [x] `BASE`: per-omics branches -> concat, Morgan fingerprint encoder,
      concat -> MLP head, raw ln(IC50), GE + Mut_CNV.
- [x] `COMPONENTS`, one switch each: `cross_attention`, `mol_graph`,
      `bilinear`, `proteomics`, `std_targets`.
- [x] Any combination by name (`base+bilinear+proteomics`); aliases `aligned`
      (= `base+mol_graph+bilinear`) and `full` (= base + every component).
- [x] One row appended per finished run; finished runs are skipped on re-run
      (`--rerun` to repeat); `--out` for scratch runs; `--list` prints the ladder.
- [x] Smoke-tested for 1 epoch on CPU (`base`, `base+std_targets`, `full`).
- [ ] Register experiment 14 in `experiments/build_results_table.py`.

**To add a component** (e.g. the pair module): add one entry to `COMPONENTS`
that overrides one new `BASE` field, and thread that field into `run_one`.
Its `base+<name>` rung and its place in `full` follow automatically.

### Phase 2 — Run the ladder (3-5 seeds)
- [ ] `python "experiments/14_benchmark_alignment/benchmark_alignment.py" --protocols grouped --seeds 42 43 44`
      (8 configs x 3 seeds; run on the GPU machine, in the background).
- [ ] Random-split run of `aligned` for the Phase 0 sanity check.
- [ ] Decide whether `cross_attention` stays as-is: with one vector per
      modality its attention weights are always 1.0, so the rung currently
      measures the projection layers, not attention. Multi-token modalities
      would make it a real component.

**Done when:** there is a table of mean +/- std and delta-vs-base per rung.
A component that does not beat base by more than the noise band is reported as
such; Mut_CNV has hurt flat models before (E32 vs E04), so base may be weak.

### Phase 3 — Cheap gate for the pair-specific idea
- [ ] Precompute per-pair features from `src/graph/hetero_graph.pt`: target
      directly mutated (0/1), min PPI hops between mutated set and target set,
      count of mutated proteins within 1 hop of a target, plus has-target and
      has-mutation flags.
- [ ] Append them to the head input behind a flag; run 3 seeds.
- [ ] Report RMSE on all pairs and separately on pairs whose drug has a target.

**Done when:** we know whether the signal exists. If there is no gain even on
the has-target subset, stop and reassess before building Phase 4.

### Phase 4 — Pair-specific attention module
New file `src/models/pair_graph_drp.py`; reuse encoders from
`mographdrp_aligned.py` rather than copying them.
- [ ] Protein encoder: `nn.Embedding(16214, d)` then 1-2 `SAGEConv` layers over
      the `interacts_with` edges only, computed once per step (full-batch).
- [ ] Precompute padded index tensors with masks: mutated proteins per cell
      line, target proteins per drug.
- [ ] Cross-attention: drug target tokens as queries, cell mutated-protein
      tokens as keys/values, masked mean-pool to one pair vector.
- [ ] Empty sets (130 drugs with no target, 6 cell lines with no mutation):
      add a learned "none" token so no row is fully masked
      (`nn.MultiheadAttention` returns NaN on an all-masked row).
- [ ] Concatenate the pair vector with the bilinear interaction vector before
      the predictor. Flag: `pair_module in {none, features, attention}`.
- [ ] GNN ablation, because GNN is the project topic and must be shown to
      matter: the same attention over raw protein embeddings with **zero**
      message-passing layers, against 1 and 2 layers. Only the layered versions
      are a GNN; the Phase 3 hand-built features are not.
- [ ] Smoke test in `__main__` like the existing module: shapes, finite output,
      gradients reach the protein encoder, and output changes when the target
      set changes with the drug fingerprint held fixed.

Mutation and target edges are inputs, not labels, so having held-out cell
lines' edges in the graph is not leakage; say so in the report.

### Phase 5 — Full ablation
- [ ] Ladder on top of the best Phase 2 row: `+ pair features`,
      `+ pair attention`, and the full model; 3-5 seeds, both protocols.
- [ ] Leave-drugs-out run of the full model (`leave_drugs_out_split`), without
      standardisation.
- [ ] Subset metrics: has-target vs no-target pairs.
- [ ] Compared baselines on the same pairs: per-drug mean, RF, MLP,
      cross-attention (E10), `HeteroIC50GNN` (E11), aligned.

### Phase 6 — Interpretability on the new model
- [ ] Repeat the experiment 15 BRAF test using attention weights: does
      attention on BRAF track measured ln(IC50) across BRAF-mutant lines,
      against 50 untrained models, including held-out cell lines?
- [ ] Repeat the experiment 16 target-removal test if time allows.

### Phase 7 — Report
- [ ] Rewrite `docs/architecture.md` for the new model. Its section 1.3 is
      wrong for the old one: two layers do not carry PPI signal to the readout.
- [ ] Regenerate `docs/results.md`; update figures in `docs/figures/`.
- [ ] Related work: MoGraphDRP, AMOGEL, MIDI, DrugVNN, SubCDR, drGT,
      Stanfield et al. 2017 (network-proximity link prediction).

## Key paths

| Path | What |
|---|---|
| `src/models/mographdrp_aligned.py` | frozen baseline with `fusion` / `head` / `drug_mode` switches |
| `src/data/experiment_utils.py` | shared loading, splits, `evaluate`, floors |
| `src/graph/hetero_graph.pt` | graph: 532 cell lines, 498 drugs, 16,214 proteins |
| `experiments/14_benchmark_alignment/` | ladder script |
| `experiments/15_diagnostics/` | why the old graph model underperforms |
| `docs/09_split_protocol_comparison.md` | protocol vs architecture gap |
| `docs/13_seed_variance_results.md` | noise band |
| `docs/CURRENTPLAN.md` | older plan, superseded by this file |

## Open questions

- Where is the MoGraphDRP repo on disk, and can it be run on their data under a
  grouped split for a system-level comparison?
- Report deadline, which decides whether Phase 6 and the leave-drugs-out run fit.
