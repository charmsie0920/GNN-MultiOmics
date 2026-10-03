# src/final_model — the final model and its ablation study

Everything for the model that goes in the final report lives here. Earlier
exploratory work stays in `experiments/` and `src/models/` and is not edited
from this folder.

| File | What |
|---|---|
| `model.py` | The model, with one switch per component (`fusion`, `head`, `drug_mode`) |
| `run_ablation.py` | Base + one component at a time; trains, evaluates, appends results |
| `results/ablation_results.csv` | One row per finished run (created on first run) |

New components for the final model (e.g. the PPI GNN with pair-specific
attention) go in this folder as new modules, and are registered in
`COMPONENTS` in `run_ablation.py`.

Shared code is imported, not copied: data loading, splits and metrics from
`src/data/experiment_utils.py`, target scaling from
`src/data/target_scaling.py`, and the fingerprint, cross-attention and
molecular-graph encoders from `src/models/`.

Run from the repository root:

```
python src/final_model/run_ablation.py --list
python src/final_model/run_ablation.py --protocols grouped --seeds 42 43 44
python src/final_model/run_ablation.py --configs base base+bilinear
```

Write-up template: `docs/14_ablation_study_results.md`. Plan: `CLAUDE.md`.
