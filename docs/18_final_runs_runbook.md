# Runbook — the remaining runs, split across two laptops

Written 2026-10-11, status updated 2026-10-13. What is left to run before the
report, who runs it, on which machine, in what order, and into which file. The
plan and the reasons are in [`CLAUDE.md`](../CLAUDE.md); the results go into
[`src/ABLATION.md`](../src/ABLATION.md).

## 0. Status on 2026-10-13

**Done**

| Block | Machine | Result | Written up in |
|---|---|---|---|
| A. Graph decision, 15 runs | Laptop A | all three configs eligible; `+mol_graph` 1.2731, `+mol_graph+pair_features` 1.2737, `+mol_graph+pair_attention` 1.3225 (test) | `src/ABLATION.md` section 5 |
| B2. Anchor, 5 runs | Laptop B | `base` 1.3290 +/- 0.0169 against 1.3252 +/- 0.0119 on laptop A | section 2 |
| B3. `HeteroIC50GNN`, 6 runs | Laptop B | 1.3203 +/- 0.0007 raw, 1.3092 +/- 0.0147 standardised | section 10 |
| B4. Random split part 1, 12 runs | Laptop B | `base` 1.0162, `aligned` 0.9941, `base+proteomics` 1.0134, `base+std_targets` 0.9423 | section 9 |
| B5. Random split part 2, 3 runs | Laptop B | `final` 0.9622 +/- 0.0058 (-0.054 against `base`, 3 of 3 seeds); 0.020 worse than `base+std_targets` alone | section 9 |
| `ldo` protocol in `run_ablation.py` | code, `e4e3517` | 168 / 36 / 36 drugs, floor 2.5500, `std_targets` refused | section 11 |
| B6. Leave-drugs-out, 10 runs | Laptop B | `base` 2.1722 +/- 0.0369, `final` without standardisation 2.0038 +/- 0.0843 (-0.169, 4 of 5 seeds; interval over drugs includes 0) | section 11 |

**Final model, chosen 2026-10-12:** `final` =
`base+proteomics+std_targets+mol_graph+pair_features`. Test 1.2737 +/- 0.0144
against `base` 1.3252 +/- 0.0119, grouped split, 5 seeds, laptop A.

**Left: nothing.** Every block of this runbook is run and written up.

Done since the last update (2026-10-12 to 2026-10-13): the `final` alias was
pushed and pulled on laptop B, B5 was run, the missing prediction file
`base+std_targets__random__seed43.npz` was committed (`edb1bf8`), the `ldo`
protocol was added, B5 was written up, and B6 was run on laptop B
(`5efc470`) and written up.

**Optional:** `HeteroIC50GNN` seeds 45 and 46 on laptop B (about a minute
each); RF and MLP with seeds under the ladder's settings (their rows are
filled from the earlier experiments); the growth-rate rung, if the supervisor
agrees.

**Open with the supervisor:** the proposal asks for a heterogeneous GNN with
GAT layers and cross-attention fusion. The eligible configs have a GCN on the
drug molecule and concatenation. See "What the proposal asks for" in
`CLAUDE.md`.

## 1. Who runs what

| Block | Machine | Runs | Time at laptop A's speed | Status |
|---|---|---|---|---|
| A. Graph decision | **Laptop A** (GTX 1650, Aditya) | 15 | about 5.7 h | **done** 2026-10-12 (5.8 h) |
| B1. Timing and data check | **Laptop B** (teammate) | 1 short | 2 min | not reported; superseded by B2 |
| B2. Anchor: `base`, grouped split | Laptop B | 5 | about 0.8 h | **done** 2026-10-11 |
| B3. `HeteroIC50GNN` baseline | Laptop B | 6 | about 10.4 h | **done** 2026-10-11 |
| B4. Random split, part 1 | Laptop B | 12 | about 2.5 h | **done** 2026-10-11 |
| B5. Random split, part 2 (final model) | Laptop B | 3 | about 1 h | **done** 2026-10-12 (2.2 h on laptop B) |
| B6. Leave-drugs-out | Laptop B | 10 | about 2.4 h | **done** 2026-10-13 (3.7 h on laptop B) |

**Measured 2026-10-11: laptop B is slower on the ladder, 0.63x the speed of
laptop A** (871 s against 546 s per `base` run), so its ladder blocks take about
1.6x the times in this table. `HeteroIC50GNN` is the exception: about a minute
per run there. B6 is self-contained, so it can run on either laptop; laptop A is
the faster choice.

**Why this split is safe.** Results are not bit-identical across hardware. So
two numbers that are subtracted from each other, or sit in one table, must
come from the same machine. Each block above is complete on one machine: it
contains its own `base` rows. Block A is compared with the `base` and
candidate rows already produced on laptop A, so it stays on laptop A.

## 2. Rules for both laptops

- **One commit, clean tree.** Check `git status` before a launch. Write the
  commit hash down with the results.
- **Laptop B never writes into an existing results CSV.** Always pass `--out`
  with the file name given below. The ladder script skips rows that are
  already in the file, so running into `ablation_results.csv` or
  `pair_gate_results.csv` would silently mix the two machines.
- **Do not hand-edit or delete a results CSV.** Re-running the same command
  resumes after an interruption; only the run in progress is lost.
- **Keep the laptop plugged in and awake.** Locking the screen is fine. A
  suspend or shutdown ends the run in progress.
- **Every ladder row must show** `n_pairs` = 111799 and `per_drug_mean_rmse`
  = 1.4889 under the grouped split. If not, the data differs: stop.
- **When a block is done:** commit the CSV and its `_predictions/` folder,
  push, and say so. The write-up (`src/ABLATION.md`, `docs/results.md`,
  `CLAUDE.md`) is done on laptop A from those files.

## 3. Laptop B setup (once)

1. Pull the branch `dev` at the commit agreed with laptop A.
2. Python 3.14, then `pip install -r requirements.txt`. Laptop A has
   `torch 2.13.0+cu130` and `torch_geometric 2.8.0.post1`. Install the CUDA
   build of torch and check that this prints `True`:

   ```
   python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
   ```

3. Copy these five files from laptop A into the same paths. They are not in
   the repository. **Copy them; do not regenerate them.**

   | File | SHA-256, first 16 characters |
   |---|---|
   | `data/processed/transcriptomics_pca.csv` | `d5a750d7a27e2788` |
   | `data/processed/genomics_pca.csv` | `3b7d073910e26e77` |
   | `data/processed/proteomics_pca.csv` | `3dd58430fa88c922` |
   | `data/processed/aligned/gdsc2_response_master.csv` | `d72b916fd5a6d35d` |
   | `data/raw/pubchem/gdsc_drug_smiles.csv` | `670f44617858e71b` |

   `src/graph/hetero_graph.pt` comes with the repository (`d06f1c61609e7758`).
   Check with `sha256sum <file>` on Linux or macOS, or
   `certutil -hashfile <file> SHA256` on Windows.

All commands below are run from the repository root.

## 4. Laptop A — block A, graph decision (done 2026-10-12)

**Result.** Run 2026-10-11 21:51 to 2026-10-12 03:40 from commit `2e425a9`.
`C` is the candidate `base+proteomics+std_targets`; 5 seeds; floor 1.4889:

| Config | Validation | Test | Eligible (1.2754 or better on validation) |
|---|---|---|---|
| `C` | 1.2454 +/- 0.0037 | 1.2797 +/- 0.0192 | no graph component |
| `C+mol_graph` | 1.2559 +/- 0.0039 | 1.2731 +/- 0.0182 | yes |
| `C+mol_graph+pair_features` | 1.2441 +/- 0.0102 | 1.2737 +/- 0.0144 | yes, lowest validation |
| `C+mol_graph+pair_attention` | 1.2745 +/- 0.0159 | 1.3225 +/- 0.0285 | by 0.0009 only; worse than `C` on every seed |

**Decided 2026-10-12: `C+mol_graph+pair_features` is the final model**, the
rule's tie-break (lowest validation RMSE). It is named `final` in
`run_ablation.py`. The plan as written before the runs follows.

The project objectives need a graph component in the final model. The
candidate `base+proteomics+std_targets` (test 1.2797 +/- 0.0192) has none.
These runs add the built graph components to it:

| Config | Adds |
|---|---|
| `base+proteomics+std_targets+mol_graph` | a GCN over each drug's atoms and bonds (a GNN) |
| `base+proteomics+std_targets+mol_graph+pair_features` | that, plus the PPI network as 5 topology features |
| `base+proteomics+std_targets+mol_graph+pair_attention` | that, plus attention between target and mutated proteins (no message passing) |

```
setsid nohup systemd-inhibit --what=sleep:idle --why="Graph decision" \
  .venv/bin/python -u src/final_model/run_ablation.py \
  --configs base base+proteomics+std_targets \
            base+proteomics+std_targets+mol_graph \
            base+proteomics+std_targets+mol_graph+pair_features \
            base+proteomics+std_targets+mol_graph+pair_attention \
  --protocols grouped --seeds 42 43 44 45 46 \
  > graph_decision.log 2>&1 < /dev/null &
```

- Results go to `src/final_model/results/ablation_results.csv`. The 10
  `base` and candidate rows are already there and are skipped: 15 runs.
- About 68 min per seed (20 + 20 + 28 min). The runs go seed by seed, so
  after an interruption every finished seed is complete for all three.
- If time is short, drop the third config: it is the one expected to lose.
- The attention runs also save weights (about 19 MB each) in
  `ablation_results_checkpoints/`.

**Decision rule, written before the runs.** A graph-bearing config is
eligible if its validation RMSE is no more than 0.03 worse than the
candidate's (1.2454, so 1.2754 or better). Among eligible configs the team
picks the one that meets the project objective; if more than one does, the
one with the lower validation RMSE. Test RMSE is reported afterwards and is
not used to choose. If none is eligible, stop and reassess.

**Expected before the runs, from single-rung effects only:** the first two
near the candidate (about 1.28 to 1.29 on test), the third about 0.06 worse.
Measured: 1.2731, 1.2737 and 1.3225.

**After block A:** the chosen config gets the name `final` in
`run_ablation.py`, is committed, and laptop B pulls that commit before B5
and B6.

## 5. Laptop B — blocks B1 to B5 (B2 to B5 done)

### B1. Timing and data check (2 minutes) — not reported, superseded by B2

```
python src/final_model/run_ablation.py --configs base --protocols grouped \
  --max-epochs 3 --out ../timing_check.csv
```

- The output is outside the repository on purpose. Delete
  `../timing_check.csv` and `../timing_check_predictions/` after.
- Laptop A takes **2.8 s per epoch** for `base`. Divide the reported
  `fit=` seconds by 3 and compare: that ratio is the speed-up.
- The log must show 111,799 pairs and a split of 77,668 / 17,291 / 16,840.

### B2. Anchor: `base` under the grouped split — done 2026-10-11

**Result:** validation 1.2928 +/- 0.0151, test 1.3290 +/- 0.0169. The means
are within 0.004 of laptop A; single seeds differ by up to 0.028. Data
confirmed (111,799 pairs, floor 1.4889). 871 s per run against 546 s on
laptop A.

```
python -u src/final_model/run_ablation.py --configs base \
  --protocols grouped --seeds 42 43 44 45 46 \
  --out src/final_model/results/laptop_b_anchor_results.csv
```

- Purpose: it measures how far the two machines differ on the same model,
  and it is the reference for `HeteroIC50GNN` on this machine.
- Laptop A's `base`: validation 1.2961 +/- 0.0059, test 1.3252 +/- 0.0119
  (seed 42: 1.2944 / 1.3141). **Report the difference.** A gap of about 0.01
  is expected. A gap of 0.03 or more means the data or environment differs:
  stop and check section 3.

### B3. `HeteroIC50GNN` baseline — done 2026-10-11 (seeds 42 to 44)

**Result:** raw 1.3203 +/- 0.0007 (validation 1.2628), standardised
1.3092 +/- 0.0147 (validation 1.2462). About a minute per run on laptop B.
Seeds 45 and 46 are optional.

The older graph model over the PPI, drug-target and mutation graph. It is a
compared baseline, not a ladder rung. With message passing skipped in the
ladder, it is the report's only trained PPI GNN, so it needs seeds.

Two arms, raw targets and per-drug standardised targets, seeds 42 to 44.
**Always pass `--out`**: the script writes one row and overwrites its file,
and its default is a committed result.

```
for s in 42 43 44; do
  python -u experiments/07_gnn_ablation/hetero_ic50_gnn_matrix.py --seed $s \
    --out src/final_model/results/hetero_gnn/raw_seed$s.csv
  python -u experiments/07_gnn_ablation/hetero_ic50_gnn_matrix.py --seed $s \
    --standardize-targets \
    --out src/final_model/results/hetero_gnn/std_seed$s.csv
done
```

- On Windows, run the six commands one by one with the seed written in.
- About 104 min per run on laptop A. Seeds 45 and 46 can be added later.
- Laptop A, seed 42: raw validation 1.2699 / test 1.3449; standardised
  1.2506 / 1.2917. Expect values near these, not identical.
- This model trains differently from the ladder (AdamW, lr 1e-3, early
  stopping with patience 15). Say so wherever the two are compared.

### B4. Random split, part 1 — done 2026-10-11

**Result (test, 3 seeds):** `base` 1.0162 +/- 0.0011, `aligned` 0.9941 +/-
0.0047, `base+proteomics` 1.0134 +/- 0.0063, `base+std_targets` 0.9423 +/-
0.0080. Committed as `e34f390`. The one prediction file that was missing
(`base+std_targets__random__seed43.npz`) was added in `edb1bf8`.

The random pair split is a sanity check against the benchmark's protocol,
not a headline result. It is never ranked against the grouped rows.

```
python -u src/final_model/run_ablation.py \
  --configs base aligned base+proteomics base+std_targets \
  --protocols random --seeds 42 43 44 \
  --out src/final_model/results/random_split_results.csv
```

- About 50 min per seed on laptop A.
- Reference from Phase 0, same split, older script: `aligned` 0.978 +/-
  0.002. Expect best epochs near 200 here, unlike the grouped split.

### B5. Random split, part 2 — done 2026-10-12

**Result (3 seeds, committed as `c346fc9`):** `final` validation 0.9610 +/-
0.0022, test 0.9622 +/- 0.0058, against `base` 1.0162 +/- 0.0011 on the same
machine (-0.054, 3 of 3 seeds). `base+std_targets` alone is 0.020 better than
`final` on every seed. Best epoch 48 / 67 / 46; 2,649 s per run.

Pull the commit that defines `final`, then run the same command with
`final` added. Finished rows are skipped:

```
python -u src/final_model/run_ablation.py \
  --configs base aligned base+proteomics base+std_targets final \
  --protocols random --seeds 42 43 44 \
  --out src/final_model/results/random_split_results.csv
```

## 6. Block B6, leave-drugs-out — done 2026-10-13 (laptop B)

Tests the final model on drugs it has never seen. It is self-contained (it
has its own `base` rows), so it can run on either laptop; laptop A is about
1.6x faster. Both things it needed exist:

1. The final model is chosen (2026-10-12): `final`.
2. `run_ablation.py` has the `ldo` protocol since `e4e3517` (2026-10-12). It
   uses `leave_drugs_out_split` from `src/data/experiment_utils.py` and
   refuses any config with standardised targets, `final` included, because
   unseen drugs have no statistics. Checked on laptop A on 2026-10-13 with a
   1-epoch run into scratch: 168 / 36 / 36 drugs, 78,973 / 16,878 / 15,948
   pairs, mean-only floor 2.5500.

**Result (laptop B, 5 seeds, committed as `5efc470`, 3.7 h).** Test fold, 36
held-out drugs, mean-only floor 2.5500:

| Config | Validation | Test |
|---|---|---|
| `base` | 2.2021 +/- 0.0283 | 2.1722 +/- 0.0369 |
| `base+proteomics+mol_graph+pair_features` | 2.0715 +/- 0.0268 | 2.0038 +/- 0.0843 |

Delta on test -0.169, 4 of 5 seeds; the interval over the 36 drugs includes
0. Not better than the earlier MLP's single run (1.8516). Written up in
`src/ABLATION.md` section 11.

The command. The final model runs without target standardisation here:

```
python -u src/final_model/run_ablation.py \
  --configs base base+proteomics+mol_graph+pair_features \
  --protocols ldo --seeds 42 43 44 45 46 \
  --out src/final_model/results/ldo_results.csv
```

The per-drug mean floor does not apply under this split. The reference is
the earlier MLP result, 1.8516 (`docs/10_leave_drugs_out_results.md`).

## 7. Not in this runbook

- **PPI message passing under the pair attention (Stage B):** skipped,
  decided 2026-10-06. Not built.
- **The growth-rate rung:** not decided; waiting on the supervisor.
- **RF and MLP baselines:** minutes each; added when the comparison table in
  `src/ABLATION.md` section 10 is filled.

## 8. Checklist

Done:

- [x] Laptop A: results of 2026-10-07 committed and pushed
- [x] Laptop B: set up, data confirmed by the B2 anchor (111,799 pairs, floor
      1.4889)
- [x] A graph decision, 15 runs (2026-10-11 21:51 to 2026-10-12 03:40, commit
      `2e425a9`; all three configs eligible, see `src/ABLATION.md` section 5)
- [x] B2 anchor, 5 runs (2026-10-11): laptop B 1.3290 +/- 0.0169 against laptop A
      1.3252 +/- 0.0119, means within 0.004; laptop B is 0.63x the speed of A
- [x] B3 `HeteroIC50GNN`, 6 runs (2026-10-11): raw 1.3203 +/- 0.0007,
      standardised 1.3092 +/- 0.0147
- [x] B4 random split part 1, 12 runs (2026-10-11, on `origin/dev` at `e34f390`)
- [x] Write-up of A, B2, B3 and B4 in `src/ABLATION.md` (2026-10-12)
- [x] Laptop A: graph-decision results and write-up committed and pushed;
      `origin/dev` pulled (2026-10-12)
- [x] Final model chosen by the rule in section 4 and named `final`
      (2026-10-12): `base+proteomics+std_targets+mol_graph+pair_features`
- [x] Laptop A: the `final` alias and document updates committed and pushed
      (2026-10-12)
- [x] B5 random split part 2, 3 runs (laptop B, 2026-10-12): `final`
      0.9622 +/- 0.0058
- [x] Laptop B: `base+std_targets__random__seed43.npz` committed (`edb1bf8`)
- [x] `ldo` protocol added to `run_ablation.py` (`e4e3517`, 2026-10-12)
- [x] Write-up of B5 in `src/ABLATION.md` section 9 (2026-10-13)
- [x] B6 leave-drugs-out, 10 runs (laptop B, 2026-10-13, `5efc470`)
- [x] Write-up of B6 in `src/ABLATION.md` section 11; `docs/results.md`
      regenerated, unchanged (2026-10-13)

Left:

Nothing.

Not needed any more:

- B1 timing and data check: never reported; B2 measured both the data and the
  speed.
