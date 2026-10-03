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
- **Closest prior work: MIDI** (bioRxiv 2025.03.31.646490). It already uses
  drug-target knowledge with attention over genes and tests mutated vs
  wild-type cell lines. Do not claim to be first at either. Ours differs by
  explicit target tokens, pair-specific weights and PPI topology. Limitation
  to state: our component needs known targets (130 drugs have none).
  See `docs/17_midi_review.md`.
- **`HeteroIC50GNN` is no longer the final model.** It stays in the results
  table as a compared graph baseline. Its own diagnostics show why
  (`experiments/15_diagnostics/`): 95.8% of proteins unattached, 130 drugs with
  no target edge, and no pair-specific path between mutation and target.

### What makes the final model ours, not MoGraphDRP

| Component | MoGraphDRP | Ours |
|---|---|---|
| Omics | expression, mutation, methylation, pathway | GE, Mut+CNV, **proteomics** |
| Biological prior | none on the cell side | **STRING PPI + drug-target + driver-mutation edges** |
| Cell-drug interaction | bilinear attention on two pooled vectors | bilinear **plus pair-specific attention between known target proteins and mutated proteins over the PPI graph** (MIDI ranks genes per drug from structure; ours is per pair and uses targets as inputs) |
| Target scaling | raw ln(IC50) | **per-drug standardisation** (cell-line split only) |
| Evaluation | random pair split | cell-line-grouped split (also used by MIDI, so not a novelty claim), **seeds, per-drug-mean floor** |
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
- **Always update `src/ABLATION.md`.** It is the ablation write-up we present
  from, and it is hand-written (unlike `docs/results.md`).
  Whenever a ladder rung finishes, a component is added or dropped, a protocol
  or seed set changes, or any number, table or conclusion that belongs in the
  ablation story moves, put it in that file in the same session as the run —
  mean +/- std, delta vs base, the per-drug-mean floor, and whether the rung
  beat the noise band. Rejected components stay in the file, marked
  tested-and-rejected. Never leave results only in a CSV, a log or a chat
  reply.
- Per-drug standardisation cannot be used under leave-drugs-out (unseen drugs
  have no statistics).
- Hardware is a GTX 1650 (4 GB). Keep the graph full-batch and small; batch-128
  runs are slow, so launch long sweeps in the background and write results
  incrementally.

## Tasks

### Phase 0 — Validate the base (passed 2026-10-01)
- [x] Run `benchmark_alignment.py` for `aligned` under the random split,
      seeds 42-44: test RMSE **0.978 +/- 0.002** (0.976 / 0.980 / 0.980).
- [x] Compared against the MoGraphDRP repo and ported the differences into
      `MoGraphDRPAligned`: BAN bilinear head (3 heads), `(512, 128)` predictor,
      gated-sum drug fusion, halving fingerprint encoder, widening GCN with sum
      pooling, no weight decay, no early stopping (200 epochs, constant LR).
      Their training code sets no seed, so 0.9497 is one unreproducible draw.
- [x] Read MIDI in full: see `docs/17_midi_review.md`. The novelty component
      stands, with narrower wording (Direction section above).

**Done when:** `aligned` under the random split is near 0.95 (their
no-XGBoost figure), or the gap is explained by the documented data differences.

**Result:** 0.978 vs their 0.9497, a gap of 0.029 (at the edge of the noise
band). The architecture now follows their code; the remaining known differences
are data (2 PCA omics vs 4 gene-filtered, Morgan only vs 3 fingerprints, our
atom features, different pairs). Not tested: whether those explain the gap.
`--seeds` varies initialisation only; the random split is fixed. Best epoch was
195-200, so the model is still improving when their 200-epoch budget ends.
**The base is now frozen.**

### Phase 1 — Additive ablation ladder (lecturer requirement)
The lecturer requires **Base + a, Base + b, Base + c, ...**, where Base is the
simplest pipeline that functions with none of the improvements, and Base plus
everything is the final model.

**Decided 2026-10-03: all final experiments live in `src/final_model/`.**
`experiments/` is history and is not extended; `experiments/14_benchmark_alignment/`
keeps only the Phase 0 reproduction. The ladder is
`src/final_model/run_ablation.py`, which imports the frozen `MoGraphDRPAligned`
(no copy) and trains with the Phase 0 settings: Adam, lr 1e-4, weight decay 0,
200 epochs, patience 200 for both early stopping and the LR scheduler (so
neither fires), checkpoint on best validation RMSE. The earlier copy
`src/final_model/model.py` (pre-Phase-0: 4-head simple bilinear, weight decay
1e-5, early stopping) was deleted; no ladder result was produced with it.

- [x] `BASE`: per-omics branches -> concat, Morgan fingerprint encoder,
      concat -> MLP head, raw ln(IC50), GE + Mut_CNV.
- [x] `COMPONENTS`, one switch each: `cross_attention`, `mol_graph`,
      `bilinear`, `proteomics`, `std_targets`.
- [x] Any combination by name (`base+bilinear+proteomics`); aliases `aligned`
      (= `base+mol_graph+bilinear`) and `full` (= base + every component).
- [x] One row appended per finished run; finished runs are skipped on re-run
      (`--rerun` to repeat); `--out` for scratch runs; `--list` prints the ladder.
- [x] Smoke-tested for 1 epoch on CPU (`base`, `base+std_targets`, `full`),
      re-run 2026-10-03 on the frozen model.
- [x] Registered in `experiments/build_results_table.py`: a separate
      "Additive ablation ladder" section (mean +/- std, delta vs base, gain over
      the per-drug floor) reading `src/final_model/results/ablation_results.csv`,
      plus the Phase 0 random-split row from experiment 14. Kept out of the
      single-run ranking because these rows are multi-seed and mix protocols.

**Phase 1 done 2026-10-03.**

**To add a component** (e.g. the pair module): add one entry to `COMPONENTS`
that overrides one new `BASE` field, and thread that field into `run_one`.
Its `base+<name>` rung and its place in `full` follow automatically.

### Phase 2 — Run the ladder (3-5 seeds)

**Status: ready to run, nothing run yet.** The next action is the runbook
below. Do not edit the model or `run_ablation.py` before running; Phase 1
froze them.

#### Runbook (Colab)

**What runs:** 8 configs x 3 seeds = **24 runs**, grouped split only (the
random split was Phase 0). `--list` prints the 8 configs. Each run is 200
epochs. Phase 0's `aligned` took 35-62 min per run on a GPU, so budget
**~10-20 h in total**. A free Colab session ends after about 12 h, so run
**one seed per session**.

**Data:** five gitignored files (~68 MB) must be copied in. If they aren't on the
shared Drive, get them from a teammate who has a local copy. Paths are relative to the repo root:

| File | Size |
|---|---|
| `data/processed/transcriptomics_pca.csv` | 0.7 MB |
| `data/processed/genomics_pca.csv` | 0.7 MB |
| `data/processed/proteomics_pca.csv` | 0.7 MB |
| `data/processed/aligned/gdsc2_response_master.csv` | 66 MB |
| `data/raw/pubchem/gdsc_drug_smiles.csv` | 0.05 MB |

The STRING graph (`src/graph/hetero_graph.pt`) is not needed until Phase 4.

**Colab cells** (Runtime -> Change runtime type -> GPU). Set `DRIVE` to the
folder holding the data and the results CSV:

```
from google.colab import drive; drive.mount('/content/drive')
DRIVE = '/content/drive/MyDrive/FIT3161'
!git clone -b dev https://github.com/charmsie0920/GNN-MultiOmics.git /content/repo
%cd /content/repo
!pip install -q torch_geometric rdkit      # keep Colab's own torch; do not pip install -r requirements.txt
!mkdir -p data/processed/aligned data/raw/pubchem
!cp {DRIVE}/data/processed/*_pca.csv data/processed/
!cp {DRIVE}/data/processed/aligned/gdsc2_response_master.csv data/processed/aligned/
!cp {DRIVE}/data/raw/pubchem/gdsc_drug_smiles.csv data/raw/pubchem/
!git rev-parse HEAD                         # record this hash in src/ABLATION.md section 15
!python src/final_model/run_ablation.py --list
```

Then, one seed per session (change `42` to `43`, then `44`):

```
!python src/final_model/run_ablation.py --protocols grouped --seeds 42 --out {DRIVE}/ablation_results.csv
```

- **`--out` must point at Drive.** One row is appended after each finished
  run, so a disconnect loses only the run in progress. Re-running the same
  command skips rows already in the CSV and resumes. Never delete or hand-edit
  that CSV; `--rerun` repeats runs on purpose.
- Every session must write to the **same** CSV and run from the **same
  commit**. If `dev` has moved, `git checkout <hash>` first.
- Expected sanity values: `per_drug_mean_rmse` = 1.4889 on every row, and
  `best_epoch` usually late (Phase 0's was 195-200). After 200 epochs `base`
  should be below 1.4889 (comparable earlier models reached ~1.28-1.33); if it
  is not, stop and investigate. A 1-epoch `base` scores ~2.75, which is normal.

**After all 24 rows exist** (on the local machine, from the repo root):

1. Copy the Drive CSV to `src/final_model/results/ablation_results.csv`.
2. `python experiments/build_results_table.py` regenerates `docs/results.md`
   (section "Additive ablation ladder").
3. Fill `src/ABLATION.md` by hand: section 5 (mean +/- std, delta vs base,
   beyond the ~0.03 noise band?, gain over 1.4889), section 2 (seeds,
   hardware), section 6 per component (`cross_attention` with its caveat),
   section 7 (additivity), and the commit hash in section 15. A component
   that does not beat the noise band is marked tested-and-rejected.
4. Tick the boxes below, write the result under "Done when", and commit the
   CSV, `docs/results.md`, `src/ABLATION.md` and this file together.

- [ ] Seed 42 (8 runs)
- [ ] Seed 43 (8 runs)
- [ ] Seed 44 (8 runs)
- [ ] Results table regenerated and `src/ABLATION.md` filled
- [x] Random-split run of `aligned` for the Phase 0 sanity check (done in
      Phase 0).
- [x] `cross_attention` stays in the ladder as-is (decided 2026-10-03). With
      one vector per modality its attention weights are always 1.0, so the rung
      measures the projection layers, not attention, and it *replaces* the
      per-omics branches rather than adding to them. Report it with that
      caveat; a gain cannot be credited to attention, and no gain does not
      show attention fails. Multi-token modalities would make it a real
      component (possible later rung).

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
| `src/final_model/run_ablation.py` | the ablation ladder; all final experiments run from here |
| `src/final_model/results/ablation_results.csv` | ladder results, one row per run |
| `experiments/14_benchmark_alignment/` | Phase 0 reproduction only (random split) |
| `experiments/15_diagnostics/` | why the old graph model underperforms |
| `docs/09_split_protocol_comparison.md` | protocol vs architecture gap |
| `src/ABLATION.md` | ablation write-up; keep it current (see Rules) |
| `docs/13_seed_variance_results.md` | noise band |
| `docs/17_midi_review.md` | MIDI vs our novelty component |
| `docs/CURRENTPLAN.md` | older plan, superseded by this file |

## Open questions

- Where is the MoGraphDRP repo on disk, and can it be run on their data under a
  grouped split for a system-level comparison?
- Report deadline, which decides whether Phase 6 and the leave-drugs-out run fit.
