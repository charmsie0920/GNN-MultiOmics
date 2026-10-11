# Runbook — the remaining runs, split across two laptops

Written 2026-10-11. What is left to run before the report, who runs it, on
which machine, in what order, and into which file. The plan and the reasons
are in [`CLAUDE.md`](../CLAUDE.md); the results go into
[`src/ABLATION.md`](../src/ABLATION.md).

## 1. Who runs what

| Block | Machine | Runs | Time at laptop A's speed | Can start |
|---|---|---|---|---|
| A. Graph decision | **Laptop A** (GTX 1650, Aditya) | 15 | about 5.7 h | now |
| B1. Timing and data check | **Laptop B** (teammate) | 1 short | 2 min | now |
| B2. Anchor: `base`, grouped split | Laptop B | 5 | about 0.8 h | now |
| B3. `HeteroIC50GNN` baseline | Laptop B | 6 | about 10.4 h | now |
| B4. Random split, part 1 | Laptop B | 12 | about 2.5 h | now |
| B5. Random split, part 2 (final model) | Laptop B | 3 | about 1 h | after A is decided |
| B6. Leave-drugs-out | Laptop B | 10 | about 2.4 h | after A is decided, and after the code in section 6 exists |

Laptop B's times will be shorter. How much is unknown until B1 is run; divide
by the speed-up it gives.

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

1. Pull the branch `aditya` at the commit agreed with laptop A.
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

## 4. Laptop A — block A, graph decision

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

**Expected, from single-rung effects only:** the first two near the
candidate (about 1.28 to 1.29 on test), the third about 0.06 worse. These
are estimates, not measurements.

**After block A:** the chosen config gets the name `final` in
`run_ablation.py`, is committed, and laptop B pulls that commit before B5
and B6.

## 5. Laptop B — blocks B1 to B5

### B1. Timing and data check (2 minutes)

```
python src/final_model/run_ablation.py --configs base --protocols grouped \
  --max-epochs 3 --out ../timing_check.csv
```

- The output is outside the repository on purpose. Delete
  `../timing_check.csv` and `../timing_check_predictions/` after.
- Laptop A takes **2.8 s per epoch** for `base`. Divide the reported
  `fit=` seconds by 3 and compare: that ratio is the speed-up.
- The log must show 111,799 pairs and a split of 77,668 / 17,291 / 16,840.

### B2. Anchor: `base` under the grouped split

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

### B3. `HeteroIC50GNN` baseline

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

### B4. Random split, part 1

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

### B5. Random split, part 2 — after block A is decided

Pull the commit that defines `final`, then run the same command with
`final` added. Finished rows are skipped:

```
python -u src/final_model/run_ablation.py \
  --configs base aligned base+proteomics base+std_targets final \
  --protocols random --seeds 42 43 44 \
  --out src/final_model/results/random_split_results.csv
```

## 6. Laptop B — block B6, leave-drugs-out (not runnable yet)

Tests the final model on drugs it has never seen. Two things must exist
first, and neither does on 2026-10-11:

1. The final model must be chosen (block A).
2. `run_ablation.py` only knows the `grouped` and `random` protocols. A
   `ldo` protocol has to be added (using `leave_drugs_out_split` from
   `src/data/experiment_utils.py`), with standardised targets refused under
   it, because unseen drugs have no statistics.

Planned command, once both are done. The final model runs without target
standardisation here:

```
python -u src/final_model/run_ablation.py \
  --configs base <final model without std_targets> \
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

- [ ] Laptop A: results of 2026-10-07 committed and pushed
- [ ] Laptop B: setup done, five data files verified (section 3)
- [ ] B1 timing and data check
- [ ] A graph decision, 15 runs
- [ ] B2 anchor, 5 runs; difference from laptop A reported
- [ ] B3 `HeteroIC50GNN`, 6 runs
- [ ] B4 random split part 1, 12 runs
- [ ] Final model chosen by the rule in section 4 and named `final`
- [ ] B5 random split part 2, 3 runs
- [ ] `ldo` protocol added to `run_ablation.py`
- [ ] B6 leave-drugs-out, 10 runs
- [ ] Write-up updated on laptop A
