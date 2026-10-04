# src/final_model — the final model and its ablation study

Everything for the model that goes in the final report lives here. Earlier
exploratory work stays in `experiments/` and is not edited from this folder;
`experiments/14_benchmark_alignment/` holds the Phase 0 reproduction only.

| File | What |
|---|---|
| `run_ablation.py` | Base + one component at a time; trains, evaluates, appends results |
| `results/ablation_results.csv` | One row per finished run (created on first run) |
| `results/<csv name>_predictions/` | Validation and test predictions of each run, saved since Phase 3 |
| `pair_gate_report.py` | Phase 3 gate: subset RMSE and paired bootstrap from saved predictions |
| `results/pair_gate_results.csv` | Phase 3 gate runs (`base` re-run and `base+pair_features`) |

The model itself is the frozen baseline `src/models/mographdrp_aligned.py`,
imported, never copied. Training follows Phase 0: Adam, lr 1e-4, no weight
decay, 200 epochs at a constant rate, checkpoint on best validation RMSE.

New model components (e.g. the PPI GNN with pair-specific attention, planned
as `src/models/pair_graph_drp.py`) go in `src/models/` beside the base, and are
registered in `COMPONENTS` in `run_ablation.py`: one entry overriding one new `BASE` field,
threaded into `run_one`. Its `base+<name>` rung follows automatically. `full`
is pinned to the five Phase 2 components, so its existing rows keep their
meaning.

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

After a run: `python experiments/build_results_table.py` regenerates
`docs/results.md`, and the numbers go into `src/ABLATION.md` by hand.
Plan: `CLAUDE.md`.
