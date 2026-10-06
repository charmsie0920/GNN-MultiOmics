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
  Since 2026-10-05 it is built and tested in two stages: the attention first,
  then the PPI message passing as its own rung (Phase 4). Stage A (attention
  alone) failed its pass rule on 2026-10-05: worse than `base` on all pairs,
  and the attention did not learn to find a mutated target. The component as
  built is not supported by the results; see Phase 4. Stage B (the message
  passing) is skipped for now (decided 2026-10-06), so the ladder has no
  trained rung that uses the PPI topology.
- **Closest prior work: MIDI** (bioRxiv 2025.03.31.646490). It already uses
  drug-target knowledge with attention over genes and tests mutated vs
  wild-type cell lines. Do not claim to be first at either. Ours differs by
  explicit target tokens, pair-specific weights and PPI topology. Limitation
  to state: our component needs known targets. Of the 240 drugs in the
  training pairs, 71 have none (28.25% of pairs); the graph's "130 drugs with
  no target" counts all 498 drug nodes.
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

**How the table stands against the results (2026-10-06).** It is the
intended design, not what the ablation supports:

- Omics: proteomics is borderline (-0.030, 3 of 3 seeds), the only gain.
- Biological prior: no supporting rung. Hop features had no effect; pair
  attention was rejected; message passing (Stage B) was not run.
- Cell-drug interaction: bilinear is rejected as configured; pair attention
  is rejected as configured.
- Target scaling: no effect (-0.012).
- Evaluation: done as described (grouped split, seeds, floor).
- Interpretability: on hold, see Phase 6.

The table is rewritten once the final configuration is chosen.

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

**Where things stand (2026-10-06).** Phases 0-3 are done. Phase 4: Stage A
was run and pair attention is tested-and-rejected as configured; **Stage B is
skipped for now** (decided 2026-10-06). Next is Phase 5, and it needs three
decisions first, all in "Open questions":

1. The final configuration: which rungs make up "Base + everything". Only
   `base+proteomics` has a consistent gain so far.
2. Whether to run the `growth_rate` rung (candidate, not run).
3. What the report's graph result is, now that the ladder has no trained
   rung with message passing. Planned answer: `HeteroIC50GNN` with 3-5 seeds
   as a compared baseline (Phase 5).

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
      (= `base+mol_graph+bilinear`) and `full` (= base + the five Phase 2
      components; pinned since 2026-10-04).
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
Its `base+<name>` rung follows automatically. Since 2026-10-05 two components
may set the same field when they are alternatives (`pair_features`,
`pair_attention`); a config combining them is refused. `full` does not: since
2026-10-04 it is pinned to the five Phase 2 components (`PHASE2_COMPONENTS`),
so adding a component cannot change what the existing `full` rows mean. The
final model gets its own name when Phase 5 defines it.

### Phase 2 — Run the ladder (3-5 seeds)

**Status: done 2026-10-03.** All 24 runs finished locally in 6.5 h (15:33 to
22:05) at commit `66205f8`. Results and verdicts are under "Done when" below
and in `src/ABLATION.md`.

#### Runbook (local GTX 1650; decided 2026-10-03, replaces the Colab plan)

**What runs:** 8 configs x 3 seeds = **24 runs**, grouped split only (the
random split was Phase 0). `--list` prints the 8 configs. Each run is 200
epochs.

**Why local, not Colab:** the model is small (peak GPU memory 129 MiB, GPU at
66-80% utilisation on `aligned`), so the laptop is not the bottleneck one
would expect. Timed on the real training loop (200 batches per config,
extrapolated to 200 epochs): `base` ~9 min, `aligned` ~21 min, `full` ~39 min,
**~2.2 h per seed, ~6.5 h for all 24** (budget 7-10 h for thermal throttling).
Colab was not timed; its free tier also needs the browser tab open, so it does
not remove the need to keep the laptop on. The data files and `.venv`
(CUDA torch, PyG, RDKit) are already on this machine.

**Launch** from the repo root, at one commit, on a clean tree:

```
setsid nohup systemd-inhibit --what=sleep:idle --why="Phase 2 ladder" \
  .venv/bin/python -u src/final_model/run_ablation.py \
  --protocols grouped --seeds 42 43 44 > phase2_ladder.log 2>&1 < /dev/null &
```

- Results go to the default `src/final_model/results/ablation_results.csv`.
  One row is appended after each finished run, so an interruption loses only
  the run in progress. Re-running the same command skips rows already in the
  CSV and resumes. Never delete or hand-edit that CSV; `--rerun` repeats runs
  on purpose.
- The laptop must stay powered on and awake (plugged in). Shutdown, hibernate
  or suspend ends the run in progress. `HandleLidSwitch=ignore` is set, but a
  lid close has not been tested.
- **Keep all 24 rows on this machine and this commit.** GPU results are not
  bit-identical across hardware, so do not mix in rows from Colab or another
  laptop.
- `--seeds` varies initialisation and batch order only; the grouped split is
  fixed at seed 42. The std over seeds is training noise, not split noise.
- Progress: `tail phase2_ladder.log` (gitignored) and the CSV row count.
- Expected sanity values: `per_drug_mean_rmse` = 1.4889 on every row. Under
  the grouped split `best_epoch` is mostly early (10-60; `base` 72-122,
  `base+proteomics` 102-180), unlike Phase 0's 195-200 under the random split;
  only `base+bilinear` peaked at the end. After 200 epochs `base`
  should be below 1.4889 (comparable earlier models reached ~1.28-1.33); if it
  is not, stop and investigate. A 1-epoch `base` scores ~2.75, which is normal.

**After all 24 rows exist** (from the repo root):

1. `python experiments/build_results_table.py` regenerates `docs/results.md`
   (section "Additive ablation ladder").
2. Fill `src/ABLATION.md` by hand: section 5 (mean +/- std, delta vs base,
   beyond the ~0.03 noise band?, gain over 1.4889), section 2 (seeds,
   hardware), section 6 per component (`cross_attention` with its caveat),
   section 7 (additivity), and the commit hash in section 15. A component
   that does not beat the noise band is marked tested-and-rejected.
3. Tick the boxes below, write the result under "Done when", and commit the
   CSV, `docs/results.md`, `src/ABLATION.md` and this file together.

- [x] Seed 42 (8 runs) — done 2026-10-03, 2.2 h on the GTX 1650
- [x] Seed 43 (8 runs) — done 2026-10-03
- [x] Seed 44 (8 runs) — done 2026-10-03
- [x] Results table regenerated and `src/ABLATION.md` filled (2026-10-03)
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

**Result (2026-10-03, grouped split, seeds 42-44, floor 1.4889):**

| Config | Test RMSE | Delta vs base | Seeds better than base | Verdict |
|---|---|---|---|---|
| `base` | 1.3232 +/- 0.0105 | 0 | | |
| `base+proteomics` | 1.2930 +/- 0.0169 | -0.0302 | 3 / 3 | borderline, at the band edge |
| `base+std_targets` | 1.3108 +/- 0.0195 | -0.0124 | 2 / 3 | tested-and-rejected (no effect) |
| `base+mol_graph` | 1.3319 +/- 0.0223 | +0.0087 | 1 / 3 | tested-and-rejected (no effect) |
| `base+cross_attention` | 1.4144 +/- 0.0040 | +0.0912 | 0 / 3 | tested-and-rejected (hurts) |
| `base+bilinear` | 1.8343 +/- 0.4946 | +0.5111 | 0 / 3 | tested-and-rejected as configured (unstable) |
| `aligned` | 1.3416 +/- 0.0088 | +0.0184 | 0 / 3 | no measurable difference from base |
| `full` | 1.3394 +/- 0.0081 | +0.0162 | 0 / 3 | not better than base |

- **`full` does not beat `base`**, and no single component clearly beats the
  noise band. Proteomics is the only consistent gain.
- **`base+bilinear` is not converged**: best epoch 187-200, per-seed RMSE
  2.40 / 1.62 / 1.48, PCC 0.85 with a large absolute error. The same head
  trains normally with the molecular graph (`aligned`). Cause not diagnosed.
- **Validation and test rank the rungs differently** (`aligned` and
  `std_targets` look good on validation only), so rankings from this one
  split are weak.
- **Open before Phase 5:** which row Phase 5 builds on. Best on test is
  `base+proteomics`; best on validation is `base+std_targets`. Not decided.
- `src/ABLATION.md` previously listed `base` at 544,001 parameters; the frozen
  model has 3,184,257 (`full`: 4,270,706).

### Phase 3 — Cheap gate for the pair-specific idea

**Status: run 2026-10-04, gate not passed. Reassessed 2026-10-05: Phase 4
goes ahead in two stages** (see Phase 4).

- [x] Precompute per-pair features from `src/graph/hetero_graph.pt`
      (`src/data/pair_features.py`): target directly mutated (0/1),
      1 / (1 + min PPI hops) between mutated set and target set, log count of
      mutated proteins within 1 hop of a target, plus has-target and
      has-mutation flags.
- [x] Append them to the head input behind a flag: component `pair_features`
      (`pair_module="features"`), model `PairGraphDRP` in
      `src/models/pair_graph_drp.py` (subclass of the frozen base; only the
      first predictor layer is widened). 5 seeds (42-46), with `base` re-run
      beside it, into `src/final_model/results/pair_gate_results.csv`.
- [x] Report RMSE on all pairs and separately on pairs whose drug has a target
      (`src/final_model/pair_gate_report.py`), plus the direct-hit subset.

**Done when:** we know whether the signal exists. If there is no gain even on
the has-target subset, stop and reassess before building Phase 4.

**Why the gate is read on a subset.** A target is directly mutated in 1.26% of
pairs (train 989, val 213, test 211), so all-pairs RMSE cannot move by more
than about 0.005. Rule fixed before the run: direct-hit RMSE better on every
seed on validation, and a test bootstrap interval that excludes 0.

**Result (2026-10-04, grouped split, commit `cc303ee`, 10 runs, 98 min).**
Seeds 42-44 first; 45-46 added the same day to settle the direct-hit result,
so 5 seeds per arm. Test fold; intervals are a paired bootstrap over test cell
lines and a t-interval over the per-seed differences:

| Subset | Pairs | Floor | `base` | `base+pair_features` | Delta | Seeds better | 95% CI, cell lines | 95% CI, seeds |
|---|---|---|---|---|---|---|---|---|
| all | 16,840 | 1.4889 | 1.3252 +/- 0.0119 | 1.3216 +/- 0.0161 | -0.0035 | 3 / 5 | -0.016 to +0.009 | -0.012 to +0.005 |
| drug has a target | 12,119 | 1.4494 | 1.2632 +/- 0.0163 | 1.2568 +/- 0.0170 | -0.0064 | 4 / 5 | -0.022 to +0.009 | -0.020 to +0.007 |
| drug has no target (control) | 4,721 | 1.5857 | 1.4722 +/- 0.0098 | 1.4751 +/- 0.0167 | +0.0028 | 3 / 5 | -0.011 to +0.018 | -0.014 to +0.020 |
| target directly mutated | 211 | 1.9609 | 1.7157 +/- 0.1212 | 1.5771 +/- 0.1528 | -0.1386 | 4 / 5 | -0.304 to +0.067 | -0.374 to +0.097 |
| nearest mutation 1 hop away | 8,688 | 1.4699 | 1.2924 +/- 0.0149 | 1.2899 +/- 0.0164 | -0.0025 | 3 / 5 | -0.020 to +0.014 | -0.010 to +0.005 |
| nearest mutation 2+ hops away | 3,051 | 1.3578 | 1.1455 +/- 0.0147 | 1.1421 +/- 0.0120 | -0.0034 | 3 / 5 | -0.017 to +0.010 | -0.019 to +0.013 |

- **No gain on all pairs or on the has-target subset** (-0.0064, inside both
  intervals). `base+pair_features` is tested-and-rejected on the headline
  metric.
- **`base` is biased on direct-hit pairs:** it predicts ln(IC50) about 0.4 too
  high on every seed (test +0.404, val +0.454). The features cut that to
  +0.028 / +0.097, on all 5 seeds and both folds.
- **The direct-hit RMSE gain is likely but not established.** Better on 4/5
  seeds on both folds (seed 43 worse on both). Test -0.139, both intervals
  include 0. Validation -0.236, cell-line interval -0.357 to -0.094, seed
  interval -0.518 to +0.047. The rule is not met. With 3 seeds the test delta
  was -0.046; the two added seeds both improved.
- **Nothing at 1 hop or more**, no different from the no-target control. The
  hand-built PPI proximity features add nothing; any signal is in "a target
  is itself mutated".
- The `base` re-run matched the Phase 2 rows exactly on seeds 42-44, so
  same-machine reruns are reproducible. The gate's 5-seed `base` is
  1.3252 +/- 0.0119 (Phase 2, 3 seeds: 1.3232 +/- 0.0105).
- `full` is now pinned to the five Phase 2 components; runs save their
  validation and test predictions to `<csv name>_predictions/`.

### Phase 4 — Pair-specific attention module (staged; decided 2026-10-05)

**Status: Stage A built and run 2026-10-05, pass rule failed; pair attention
is tested-and-rejected as configured. Stage B skipped for now (decided
2026-10-06): not built, not run.** Phase 4 is closed unless Stage B is
reopened.

**Decision (2026-10-05): build Phase 4 in two stages, not all at once.**
Phase 3 found signal where a drug's target is itself mutated, but none from
generic PPI proximity. So the attention between target and mutated proteins
is tested first, without message passing. The PPI message-passing layers are
added after it, as their own rung. Reasons:

- Building both together would make any gain or loss unattributable, which
  breaks "one change per rung".
- The zero-layer version is not extra work: it is the first row of the GNN
  depth ablation the original plan already required.
- It gives an early exit. If attention cannot beat the Phase 3 flag, the
  message-passing layers are unlikely to rescue it.
- Attention, padding and empty-set handling get debugged before the heavier
  full-batch pass over 473,860 PPI edges is added.

Why attention could beat the flag: the direct-hit effect differs by drug
(-1.8 to +0.4 z on training rows), which one yes/no flag cannot express. The
Phase 3 hop features were gene-agnostic counts, so they rule out generic
proximity, not a specific mutated neighbour mattering for a specific drug.

**Expectations, written before building.** Overall RMSE will not move: direct
hits are 1.26% of pairs. The bar at each stage is a subset comparison with the
stage below it. The likely outcome of Stage B is no difference between 0, 1
and 2 layers, which is then reported as message passing tested-and-rejected.

#### Shared setup (both stages)

Extend `src/models/pair_graph_drp.py` (`PairGraphDRP` already subclasses the
frozen base and widens only the first predictor layer). Reuse
`src/data/pair_features.py` for the graph loading and ID mapping, and keep
`src/models/mographdrp_aligned.py` untouched.

**Built 2026-10-05.** Not run for results yet.

- [x] Protein embeddings: `nn.Embedding(16216, 64)` in `PairAttention`, the
      16,214 proteins plus two learned "none" tokens.
- [x] Padded index tensors with masks: `build_pair_sets` in
      `src/data/pair_features.py`. Targets 240 x 7, mutated proteins
      531 x 108 (median 6), mapped by node ID. Checked row for row against
      the Phase 3 direct-hit flag on all 111,799 pairs.
- [x] Cross-attention: drug target tokens as queries, cell mutated-protein
      tokens as keys/values, masked mean-pool over targets to one pair vector.
      One `nn.MultiheadAttention` layer, 4 heads, attention dropout 0.
- [x] Empty sets: a drug with no target (71 of 240) queries with a "no
      target" token. The "no mutation" token is a key for **every** cell line,
      not only the 6 with no mutation: softmax weights sum to 1, so a target
      with no relevant mutation needs somewhere to attend.
- [x] Concatenate the pair vector with the base's interaction vector before
      the predictor.
- [x] Attention weights on request:
      `PairGraphDRP.attention_weights(cell_codes, drug_codes)`. Attention runs
      also save their best-epoch weights to `<csv name>_checkpoints/` (17 MB
      per run, committed with the results), because attention cannot be recovered from saved
      predictions and Phase 6 needs the trained models.
- [x] Registered in `run_ablation.py` as `pair_attention`. Components that set
      the same `BASE` field are alternatives: each is its own rung, and
      `resolve_config` refuses a config that combines two of them. Stage B
      variants fit this without a new mechanism.
- [x] Smoke test in `__main__` (`python -m src.models.pair_graph_drp`): shapes,
      finite output in train and eval mode, no NaN with an empty target set,
      an empty mutation set or both, weights sum to 1 and are 0 on padding,
      gradients reach exactly the proteins in the batch and both "none"
      tokens, and output changes when only the target set or only the mutated
      set changes.
- [x] 1-epoch CPU check: `base` still gives val 2.688133 / test 2.747617.
      `base+pair_features` gives 2.554799 / 2.649105 before and after (second
      reference, recorded 2026-10-05). The Phase 3 feature matrix is
      byte-identical after the refactor.

**What the graph gives the module (measured 2026-10-05):**

- `protein.x` is all zeros: there are no protein features. Embeddings are
  learned from the response loss alone.
- Only 634 of the 16,214 proteins are a target or a mutation of any pair
  (148 targets, 536 mutated, 50 both). In Stage A the other rows of the
  embedding table never receive a gradient; the table stays full size so
  Stage B changes only the layers. Report both counts: +1,087,232 parameters,
  at most 40,704 of the embedding's 1,037,824 trained.
- 15 of the 251 proteins mutated in test cell lines are mutated in no training
  cell line (102 of the 501 training ones are seen in one cell line only).
  Their tokens are untrained at test time. A Stage A limitation to state.
- Direct hits are concentrated: 85 drugs have any, drug `1931` has 372 of the
  1,413, and the median such drug has 6 in the training fold.
- Nothing tells the module that a target and a mutated protein are the same
  protein; it has to learn that from 989 training direct hits. Check it from
  the saved weights after the run, before interpreting a failure.
- On the rationale above: the predictor already sees the direct-hit flag
  beside the drug vector, so `base+pair_features` can express a drug-specific
  effect. What attention adds is which protein is mutated.

**Speed:** look tokens up with `F.embedding(index, states)`, not
`states[index]`. The indexing form's backward pass cost 7.6 ms a step on this
GPU against 0.8 ms, and made the rung 3.7x slower than `base`. Stage B must
keep the `F.embedding` form when `protein_states()` returns message-passed
states.

#### Stage A — attention over raw protein embeddings, no message passing

- [x] Config `base+pair_attention` (`pair_module="attention"`, 0 layers).
      4,271,489 parameters. Timed over 3 epochs on the GTX 1650: 5.4 s per
      epoch (`base` 2.8 s), so about 18 min per run and 1.5 h for 5 seeds;
      peak GPU memory 227 MiB.
- [x] Run 5 seeds (42-46) into the gate CSV, so `base` and
      `base+pair_features` are skipped and reused. Commit first: the runbook
      needs a clean tree at one commit. Done 2026-10-05 at `9c7ba4f`,
      17.2 min per run.

```
setsid nohup systemd-inhibit --what=sleep:idle --why="Phase 4 stage A" \
  .venv/bin/python -u src/final_model/run_ablation.py \
  --configs base base+pair_features base+pair_attention --protocols grouped \
  --seeds 42 43 44 45 46 --out src/final_model/results/pair_gate_results.csv \
  > phase4_stage_a.log 2>&1 < /dev/null &
```

- [x] Report against both baselines. `--tag` names the output files and
      `--configs` fixes which configs share the bootstrap's random stream, so
      neither call touches the committed Phase 3 files. Those are reproduced
      by `--configs base base+pair_features` (checked byte-identical
      2026-10-05). Without `--configs`, the third config shifts every
      bootstrap draw and the published Phase 3 intervals change.

```
python src/final_model/pair_gate_report.py --tag attention_vs_features \
  --configs base+pair_features base+pair_attention --baseline base+pair_features
python src/final_model/pair_gate_report.py --tag attention_vs_base \
  --configs base base+pair_attention
```

- [x] Before reading the result: from the checkpoints, check that on
      direct-hit pairs the target's attention goes to the matching mutated
      protein rather than to the "no mutation" token
      (`python src/final_model/pair_attention_check.py`).
- [x] **Stage A passes if**, on direct-hit pairs, `base+pair_attention` beats
      `base+pair_features` on at least 4 of 5 seeds on validation, with the
      validation cell-line interval excluding 0. Test is reported after, not
      used to decide. It must also be no worse than `base` on all pairs beyond
      the ~0.03 band. Also report the has-target and hop subsets: a gain there
      would be the first evidence beyond direct hits.
- [x] Write it up in `src/ABLATION.md` section 6f in the same session.

**Result (2026-10-05, grouped split, commit `9c7ba4f`, 5 runs, 87 min,
15:48 to 17:14): Stage A failed both conditions. `base+pair_attention` is
tested-and-rejected as configured.** Floor 1.4889. Test fold unless stated:

| Subset | Pairs | `base` | `base+pair_features` | `base+pair_attention` | Δ vs `base` | Seeds better | 95% CI, cell lines |
|---|---|---|---|---|---|---|---|
| all | 16,840 | 1.3252 +/- 0.0119 | 1.3216 +/- 0.0161 | 1.3861 +/- 0.0142 | +0.0609 | 0 / 5 | +0.026 to +0.095 |
| drug has a target | 12,119 | 1.2632 +/- 0.0163 | 1.2568 +/- 0.0170 | 1.3327 +/- 0.0214 | +0.0695 | 0 / 5 | +0.035 to +0.104 |
| drug has no target (control) | 4,721 | 1.4722 +/- 0.0098 | 1.4751 +/- 0.0167 | 1.5143 +/- 0.0120 | +0.0421 | 0 / 5 | -0.001 to +0.088 |
| target directly mutated | 211 | 1.7157 +/- 0.1212 | 1.5771 +/- 0.1528 | 1.6164 +/- 0.1369 | -0.0993 | 3 / 5 | -0.306 to +0.104 |
| nearest mutation 1 hop away | 8,688 | 1.2924 +/- 0.0149 | 1.2899 +/- 0.0164 | 1.3630 +/- 0.0343 | +0.0706 | 0 / 5 | +0.037 to +0.102 |
| nearest mutation 2+ hops away | 3,051 | 1.1455 +/- 0.0147 | 1.1421 +/- 0.0120 | 1.1997 +/- 0.0078 | +0.0541 | 0 / 5 | +0.015 to +0.092 |

- **The deciding comparison:** direct-hit pairs on validation, attention
  against the features: +0.131, better on 1 seed of 5, cell-line interval
  +0.050 to +0.212. Worse, not better. (Test: +0.039, 2 of 5, interval
  includes 0.)
- **Worse than `base` on all pairs** (+0.061, 5 of 5 seeds; validation
  +0.052), and on every subset except direct hits, including the no-target
  control. So the loss is not about target-mutation matching.
- **It overfits early:** best epoch 24 / 13 / 24 / 15 / 15 (`base`:
  122 / 111 / 72 / 22 / 67).
- **The attention did not learn the match.** A mutated target gives its own
  protein 0.145-0.167 of its attention; chance is 0.126-0.148 and untrained
  modules give the same as chance. It ranks first 21-24% of the time.
- **So this rejects the module as built, not the idea.** The idea was not
  tested, because the attention never found the mutated target. Do not write
  "pair-specific attention does not help".
- On direct hits the attention sits between `base` and the features: signed
  error +0.404 / +0.028 / +0.160 on test.
- **Why (diagnosed 2026-10-05, `src/ABLATION.md` 6f):** the attention stayed
  near uniform (entropy 0.963 of maximum) and the embeddings barely left
  their random start (row norm 7.96 used vs 7.97 never-updated; unit-variance
  init, lr 1e-4, checkpoint at epoch 13-24). The pair vector is then a random
  fingerprint of the cell line: 54% of its variance is cell-only, 31%
  drug-only, 16% pair-specific. The predictor fits that fingerprint, 61 of 80
  test cell lines get worse, and validation peaks before the attention has
  learned. The match signal is only 1.6 pairs per batch.
- Cost: +1,087,232 parameters, 1.81x training time (1,041 s vs 576 s).

**Stage A failed, so:** pair attention is reported as tested-and-rejected.
Stage B was optional under this plan; **decided 2026-10-06: skip it for now**
(next section).

#### Stage B — PPI message passing (the GNN depth ablation) — skipped for now

**Decision (2026-10-06): not built and not run.** Reasons:

- It would put message passing under an attention that did not learn (near
  uniform, 0.02 above chance on the match), so a null result would say little
  about the PPI graph. It would carry Stage A's caveat: it rejects the module,
  not the idea.
- About 8.5 h of GPU time for that (estimated).
- The ablation is valid without it: every built component has its rung.

**What the ablation lacks as a result, to state in the report:**

- No trained rung uses the STRING network with message passing. The 1- and
  2-layer rows of the depth table stay "not run".
- "Does the PPI graph matter?" is answered only indirectly: the Phase 3 hop
  features (no effect, 5 seeds) and `HeteroIC50GNN` as a compared baseline.
- The "Biological prior" and "Cell-drug interaction" rows of the Direction
  table are not supported by a passing rung.

**Reopen only if** the supervisor requires a measured message-passing rung.
The specification below is kept for that case; nothing in it has been done.

- [ ] 1 and 2 `SAGEConv` layers over the `interacts_with` edges only, computed
      once per step (full-batch), feeding the same Stage A attention. Configs,
      e.g., `base+pair_attention_gnn1` and `base+pair_attention_gnn2`. The
      entry point is `PairAttention.protein_states()`: it returns the raw
      embeddings now and is the only method the layers need to change.
- [ ] Time it first on the real training loop (as in the Phase 2 runbook)
      before launching 10 runs; message passing on every step may be much
      slower than Stage A. Estimate from a benchmark of the module alone
      (2026-10-05, `SAGEConv(64, 64)` + ReLU over the 473,860 edges, batch
      128): 4.7 / 14.7 / 24.9 ms per training step for 0 / 1 / 2 layers, peak
      GPU memory under 200 MiB. That is about 12 s per epoch and 40 min per
      run for 1 layer, 19 s and 60 min for 2: **about 8.5 h for the 10 runs**
      (Stage A: 1.5 h for 5). Not a timing of the real loop.
- [ ] 5 seeds each, same CSV, same report. Compare with Stage A (0 layers) on
      direct-hit and has-target pairs, using `--baseline base+pair_attention`.
- [ ] Fill the depth-ablation table in `src/ABLATION.md` section 6f. Only the
      layered versions are a GNN; the Phase 3 features and Stage A are not.

Mutation and target edges are inputs, not labels, so having held-out cell
lines' edges in the graph is not leakage; say so in the report.

### Phase 5 — Full ablation
**Not started. Blocked on the final configuration** (see "Where things
stand"). Changed by Phase 4's outcome (2026-10-06):

- [ ] Decide the final configuration first, on validation RMSE, and give it
      its own name in `run_ablation.py` (`full` stays pinned to Phase 2).
- [ ] Ladder on top of the best Phase 2 row: `+ pair features` and the final
      model; 3-5 seeds, both protocols. `+ pair attention` is dropped from
      this ladder: it is rejected as configured and already has its rung.
- [ ] Leave-drugs-out run of the final model (`leave_drugs_out_split`),
      without standardisation.
- [ ] Subset metrics: has-target vs no-target pairs.
- [ ] Compared baselines on the same pairs: per-drug mean, RF, MLP,
      cross-attention (E10), `HeteroIC50GNN` (E11), aligned.
      **`HeteroIC50GNN` now needs 3-5 seeds** (it has a single run, 1.3301):
      with Stage B skipped it is the report's only trained GNN.

### Phase 6 — Interpretability on the new model

**On hold (2026-10-06).** The BRAF test was planned on the pair attention's
weights. Stage A's attention is near uniform, so its weights carry nothing to
interpret. Run this phase only if the final model has a component whose
attention or attribution is worth testing; otherwise report it as not done.

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
| `src/data/pair_features.py` | per-pair target-vs-mutation features, the gate's row subsets, and the padded target / mutation sets (`build_pair_sets`) |
| `src/models/pair_graph_drp.py` | `PairGraphDRP`: frozen base + a per-pair vector into the predictor, from features or from `PairAttention` |
| `src/final_model/pair_gate_report.py` | subset RMSE and paired bootstrap from saved predictions (`--configs`, `--baseline`, `--tag`) |
| `src/final_model/pair_attention_check.py` | whether the trained attention finds a mutated target, against untrained modules |
| `src/final_model/results/pair_gate_results_*attention*.csv` | Stage A subset metrics, paired deltas and the attention check |
| `src/final_model/results/pair_gate_results*.csv` | Phase 3 gate runs, subset metrics, paired deltas |
| `src/final_model/results/*_checkpoints/` | best-epoch weights of attention runs, 17 MB each; committed with the results |
| `experiments/14_benchmark_alignment/` | Phase 0 reproduction only (random split) |
| `experiments/15_diagnostics/` | why the old graph model underperforms |
| `docs/09_split_protocol_comparison.md` | protocol vs architecture gap |
| `src/ABLATION.md` | ablation write-up; keep it current (see Rules) |
| `docs/13_seed_variance_results.md` | noise band |
| `docs/17_midi_review.md` | MIDI vs our novelty component |
| `docs/CURRENTPLAN.md` | older plan, superseded by this file |

## Open questions

- **Is a report with no passing graph component, and no trained
  message-passing rung, acceptable?** This file says the GNN "must be shown
  to matter". The Phase 3 gate found no effect, Stage A of Phase 4 failed
  (2026-10-05) and Stage B is skipped (2026-10-06). A question for the
  supervisor. If the answer is no, Stage B is reopened (about 8.5 h of GPU
  for 10 runs, estimated).
- **Which rungs make up the final model?** The lecturer's ladder ends in
  "Base + everything". `full` (the five Phase 2 components) is not better
  than `base`. Only `base+proteomics` has a consistent gain (-0.030, 3 of 3
  seeds, at the band edge); two more seeds would bring it to 5. The choice
  must be made on validation RMSE. Not decided; it blocks Phase 5.
- **Pair attention: closed as rejected (2026-10-06).** The cause was
  diagnosed (Phase 4, Stage A result). Any changed version of the module
  (regularisation, an explicit same-protein signal, a smaller table) would be
  a new rung with its own name and a pass rule written first, not a retry of
  `base+pair_attention`. None is planned: direct hits cap its all-pairs gain
  at about 0.005.
- **Candidate component: growth rate (not run, not decided).** 26% of
  `base`'s squared test error is one constant per cell line (sd 0.68, same
  across seeds), far more than any direct-hit component can reach (0.005).
  Log doubling time (`data/raw/gdsc/growth_rate_20220907.csv`, 524 of 531
  cell lines) correlates +0.41 / +0.50 with that offset, and one slope fitted
  on the other held-out fold cuts test RMSE by 0.029 and validation by 0.017,
  5 of 5 seeds each. Checked and negative: tissue, cancer type, omics PCs,
  and expression of the drug's target genes. Details in `src/ABLATION.md`
  section 11. Open: is an assay covariate that is not an omics layer
  acceptable to the supervisor? If run, it is one switch (`growth_rate`),
  5 seeds, about 50 min, with its pass rule written first.
- `base+bilinear` instability is still undiagnosed. Seed 42 sat at validation
  RMSE 3.0-4.1 for most of 200 epochs with PCC 0.85. Hypothesis only: dropout
  feeding BatchNorm gives a train/eval scale mismatch.
- Which row Phase 5 builds on (`base+proteomics` on test, `base+std_targets`
  on validation) is still not decided; Phase 3 attached to `base`.
- Where is the MoGraphDRP repo on disk, and can it be run on their data under a
  grouped split for a system-level comparison?
- Report deadline, which decides whether Phase 6 and the leave-drugs-out run fit.
