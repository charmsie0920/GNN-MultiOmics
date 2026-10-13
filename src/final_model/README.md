# src/final_model — the final model and its ablation study

Everything for the model that goes in the final report lives here. Earlier
exploratory work stays in `experiments/` and is not edited from this folder;
`experiments/14_benchmark_alignment/` holds the Phase 0 reproduction only.

| File | What |
|---|---|
| `run_ablation.py` | Base + one component at a time; trains, evaluates, appends results |
| `results/ablation_results.csv` | One row per finished run (created on first run) |
| `results/<csv name>_predictions/` | Validation and test predictions of each run, saved since Phase 3 |
| `results/<csv name>_checkpoints/` | Best-epoch weights of runs with the attention module (17 MB each, committed), and of any run launched with `--save-weights` |
| `pair_gate_report.py` | Subset RMSE and paired bootstrap from saved predictions (`--configs`, `--baseline`, `--tag`) |
| `pair_attention_check.py` | Whether the trained attention finds a mutated target, against untrained modules; reads the checkpoints |
| `pair_feature_check.py` | Whether the trained final model uses the direct-hit signal, and in the measured direction, against untrained models; reads the checkpoints |
| `results/pair_gate_results.csv` | Pair-module runs: `base` re-run and `base+pair_features` (Phase 3), `base+pair_attention` (Phase 4 stage A) |
| `results/pair_gate_results_*.csv` | Subset metrics and paired deltas of those runs; the `_attention_*` files are Stage A |
| `results/ablation_results_*graph_decision*.csv` | Subset metrics and paired deltas of the graph-decision runs (2026-10-12) |
| `results/laptop_b_anchor_results.csv` | Laptop B: `base` under the grouped split, 5 seeds, to measure the difference between the two machines |
| `results/hetero_gnn/` | Laptop B: `HeteroIC50GNN` with seeds, raw and standardised targets, one file per run |
| `results/random_split_results.csv` | Laptop B: the random-split sanity check, `base`, `aligned`, the two single components and `final`, 3 seeds |
| `results/ldo_results.csv` | Laptop B: leave-drugs-out, `base` and `final` without standardisation, 5 seeds |
| `results/final_weights_results.csv` | `final` trained again on seeds 42-46 with `--save-weights` (2026-10-13), for `pair_feature_check.py`. New runs, not the ladder's rows |

The model itself is the frozen baseline `src/models/mographdrp_aligned.py`,
imported, never copied. Training follows Phase 0: Adam, lr 1e-4, no weight
decay, 200 epochs at a constant rate, checkpoint on best validation RMSE.

New model components go in `src/models/` beside the base, and are registered
in `COMPONENTS` in `run_ablation.py`: one entry overriding one `BASE` field,
threaded into `run_one`. Its `base+<name>` rung follows automatically. `full`
is pinned to the five Phase 2 components, so its existing rows keep their
meaning. `final` names the final model, chosen 2026-10-12:
`base+proteomics+std_targets+mol_graph+pair_features`. The pair-specific component is `src/models/pair_graph_drp.py`, in two
alternative forms that cannot be combined: `pair_features` (hand-built
features) and `pair_attention` (attention of a drug's target proteins over a
cell line's mutated proteins).

Shared code is imported, not copied: data loading, splits and metrics from
`src/data/experiment_utils.py`, target scaling from
`src/data/target_scaling.py`, and the encoders from `src/models/`.

Run from the repository root, on the local GTX 1650 (launch command and
runbook in `CLAUDE.md`). Re-running skips rows already in the CSV:

```
python src/final_model/run_ablation.py --list
python src/final_model/run_ablation.py --protocols grouped --seeds 42 43 44
python src/final_model/run_ablation.py --configs base base+bilinear
```

Three protocols: `grouped` (by cell line, the headline), `random` (pair
split, a sanity check against the benchmark) and `ldo` (leave drugs out).
`ldo` is not in the default sweep and refuses any config with `std_targets`,
`final` included: a held-out drug has no training statistics. Results from a
second machine, and every `random` or `ldo` block, go to their own file with
`--out`, never into `ablation_results.csv`.

After a run: `python experiments/build_results_table.py` regenerates
`docs/results.md`, and the numbers go into `src/ABLATION.md` by hand.
Plan: `CLAUDE.md`.
