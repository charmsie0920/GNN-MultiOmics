"""Backend wrapper for experiments/01_cross_attention_fusion/cross_attention_baseline.py.

Runs the script completely unmodified — loads it as a module, monkeypatches
its `print` to also stream lines to the UI, then drives its own top-level
functions (`load_omics_modalities`, `load_targets`, `build_pair_tensors`,
`grouped_split`, `make_loader`, `CrossAttentionRegressor`, `train`,
`predict`) directly, the same way `main()` does internally. All
training-orchestration and MC-dropout inference logic lives here, in the
backend, not in the script, so swapping in a different team model later is
just swapping which backend module gets registered — no science script
needs touching.
"""

from __future__ import annotations

import builtins
import gc
import importlib.util
import re
import types
from collections.abc import Callable
from pathlib import Path

import numpy as np
import torch

from backend.model_backends._common import enable_mc_dropout, rank_predictions
from backend.model_backends.base import ModelBackend
from backend.model_backends.registry import register

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT_PATH = _REPO_ROOT / "experiments" / "01_cross_attention_fusion" / "cross_attention_baseline.py"

# MC-dropout forward passes per (cell-line, drug) pair used to derive the
# confidence score — see _predict_drug_panel below.
_MC_SAMPLES = 30


def _load_module() -> types.ModuleType:
    spec = importlib.util.spec_from_file_location("cross_attention_baseline_run", _SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module spec from {_SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _prepare_and_train(module: types.ModuleType) -> dict:
    """Replicate `module.main()`'s data-prep/train steps, keeping the
    trained model and full-population tensors around for inference
    afterwards (which `main()` itself doesn't need, so it doesn't return
    them)."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(module.TORCH_SEED)

    omics, cell_ids = module.load_omics_modalities()
    keys = list(omics.keys())
    y_df = module.load_targets(cell_ids)
    gathered, drug_onehot, y, groups, n_drugs = module.build_pair_tensors(omics, cell_ids, y_df)
    del omics

    drug_levels = module.pd.factorize(y_df[module.COL_DRUG], sort=True)[1]

    # Descriptive metadata (names, not indices) for the drugs/cell lines this
    # run trains on, read straight from the same merged CSV used for targets
    # -- keyed by the same raw id strings as drug_levels/cell_ids, so there's
    # no separate id-mapping step that could drift out of sync.
    meta = module.pd.read_csv(
        module.TARGET_CSV,
        usecols=[module.COL_DRUG, module.COL_CELL_LINE, "drug_name", "putative_target", "pathway_name", "cell_line_name", "tissue", "cancer_type"],
        dtype=str,
        keep_default_na=False,
    )
    drug_meta = meta[[module.COL_DRUG, "drug_name", "putative_target", "pathway_name"]].drop_duplicates(module.COL_DRUG)
    drug_info = {
        row[module.COL_DRUG]: {
            "drug_name": row["drug_name"],
            "putative_target": row["putative_target"],
            "pathway_name": row["pathway_name"],
        }
        for row in drug_meta.to_dict("records")
    }
    cell_meta = meta[[module.COL_CELL_LINE, "cell_line_name", "tissue", "cancer_type"]].drop_duplicates(module.COL_CELL_LINE)
    cell_info = {
        row[module.COL_CELL_LINE]: {
            "cell_line_name": row["cell_line_name"],
            "tissue": row["tissue"],
            "cancer_type": row["cancer_type"],
        }
        for row in cell_meta.to_dict("records")
    }

    train_idx, val_idx, _test_idx = module.grouped_split(groups)
    omics_train = {k: v[train_idx] for k, v in gathered.items()}
    omics_val = {k: v[val_idx] for k, v in gathered.items()}

    train_loader = module.make_loader(
        omics_train, drug_onehot[train_idx], y[train_idx], keys, module.BATCH_SIZE, shuffle=True
    )

    model = module.CrossAttentionRegressor(
        n_drugs=n_drugs, hidden_dims=module.HEAD_HIDDEN_DIMS, dropout=module.HEAD_DROPOUT
    ).to(device)
    model, _best_epoch = module.train(
        model, train_loader, omics_val, drug_onehot[val_idx], y[val_idx], keys, device
    )

    # Deterministic (no MC-dropout) prediction on the held-out validation
    # set, purely to get a real, model-specific error scale (val_rmse) for
    # calibrating the confidence score against -- see `rank_predictions`.
    val_pred = module.predict(model, omics_val, drug_onehot[val_idx], keys, device, module.BATCH_SIZE)
    val_rmse = float(np.sqrt(np.mean((val_pred - y[val_idx]) ** 2)))

    return {
        "model": model,
        "keys": keys,
        "gathered": gathered,
        "cell_ids": cell_ids,
        "drug_levels": drug_levels,
        "drug_info": drug_info,
        "cell_info": cell_info,
        "n_drugs": n_drugs,
        "val_rmse": val_rmse,
        "device": device,
    }


@torch.no_grad()
def _predict_drug_panel(module: types.ModuleType, artifacts: dict, target_cell_line: str, mc_samples: int) -> list[dict]:
    """Predict ln_ic50 (mean + MC-dropout std) for every trained drug, for one cell line."""
    cell_ids = artifacts["cell_ids"]
    if target_cell_line not in cell_ids:
        raise ValueError(f"Unknown cell line for this run: {target_cell_line!r}")

    row = cell_ids.get_loc(target_cell_line)
    n_drugs = artifacts["n_drugs"]
    keys = artifacts["keys"]
    device = artifacts["device"]
    model = artifacts["model"]

    omics_row = {
        k: np.repeat(arr[row : row + 1], n_drugs, axis=0) for k, arr in artifacts["gathered"].items()
    }
    drug_onehot = np.eye(n_drugs, dtype=module.DTYPE)

    # module.predict() -- the unmodified script's own function -- calls
    # model.eval() as its first line, *inside* the call, right before the
    # forward pass. Re-enabling dropout before calling predict() (as a
    # previous version of this function did) doesn't survive that internal
    # eval() call, so dropout was silently off for every "MC" sample,
    # collapsing all of them to the same value. The only seam available
    # without touching the script is to replace what `model.eval()` itself
    # does for the duration of this sampling loop.
    original_eval = model.eval
    model.eval = lambda: enable_mc_dropout(model)
    try:
        samples = np.empty((mc_samples, n_drugs), dtype=np.float64)
        for i in range(mc_samples):
            samples[i] = module.predict(model, omics_row, drug_onehot, keys, device, module.BATCH_SIZE)
    finally:
        model.eval = original_eval
        model.eval()

    mean = samples.mean(axis=0)
    std = samples.std(axis=0)
    return [
        {"drug_id": str(drug_id), "ln_ic50_mean": float(mean[i]), "ln_ic50_std": float(std[i])}
        for i, drug_id in enumerate(artifacts["drug_levels"])
    ]


#  train()'s own epoch line: "[epoch  12] train_loss=... val_rmse=... val_pcc=..."
_EPOCH_LINE_RE = re.compile(r"^\[epoch\s*(\d+)\]")
_EPOCH_METRICS_RE = re.compile(
    r"^\[epoch\s*(\d+)\]\s+train_loss=([\d.eE+-]+)\s+val_rmse=([\d.eE+-]+)\s+val_pcc=([\d.eE+-]+)"
)


class CrossAttentionBackend(ModelBackend):
    name = "cross_attention"
    # Calibrated from the user's measured runs: ~7 min on CUDA, ~30 min on CPU.
    # Only used as the *initial* estimate before the first epoch line
    # arrives -- real progress comes from on_progress() below once training
    # starts logging epochs.
    expected_duration_seconds = {"cuda": 420.0, "cpu": 1800.0}

    def run(
        self,
        log: Callable[[str], None],
        target_cell_line: str,
        on_results: Callable[[list[dict]], None],
        on_progress: Callable[[int, int], None],
        on_training_history: Callable[[list[dict]], None],
    ) -> None:
        module = _load_module()
        on_progress(0, module.MAX_EPOCHS)
        history: list[dict] = []

        def tee_print(*args: object, sep: str = " ", end: str = "\n", **kwargs: object) -> None:
            builtins.print(*args, sep=sep, end=end, **kwargs)
            text = sep.join(str(arg) for arg in args)
            log(text)
            # train() only prints an epoch line every 5 epochs or on a val
            # improvement, so this is a best-effort, occasionally-stale
            # progress signal -- still far more honest than a fixed timer.
            match = _EPOCH_LINE_RE.match(text)
            if match:
                on_progress(int(match.group(1)), module.MAX_EPOCHS)
            metrics_match = _EPOCH_METRICS_RE.match(text)
            if metrics_match:
                history.append(
                    {
                        "epoch": int(metrics_match.group(1)),
                        "train_loss": float(metrics_match.group(2)),
                        "val_rmse": float(metrics_match.group(3)),
                        "val_pcc": float(metrics_match.group(4)),
                    }
                )

        # Every function in `module` resolves the bare name `print` against
        # the module's own globals, so this alone routes all of the
        # script's existing progress prints (load/split/epoch logs) to the
        # UI too, with no changes to the script itself.
        module.print = tee_print

        try:
            artifacts = _prepare_and_train(module)
            on_training_history(history)
            log(f"[confidence] validation RMSE (ln IC50) = {artifacts['val_rmse']:.4f} -- used as the confidence scale")
            cell_info = artifacts["cell_info"].get(target_cell_line, {})
            log(
                f"[cell-line] {target_cell_line} -> {cell_info.get('cell_line_name', 'unknown')} "
                f"({cell_info.get('tissue', 'unknown')}, {cell_info.get('cancer_type', 'unknown')})"
            )
            log(f"Running inference for target cell line: {target_cell_line}")
            raw_predictions = _predict_drug_panel(module, artifacts, target_cell_line, mc_samples=_MC_SAMPLES)
            on_results(rank_predictions(raw_predictions, artifacts["val_rmse"], artifacts["drug_info"]))
        finally:
            # The backend process stays alive across runs (each run just
            # loads a fresh copy of the script), so the trained model,
            # optimizer, and full omics tensors from *this* run would
            # otherwise linger referenced by the Python/CUDA caching
            # allocator until the next run's allocations happen to reclaim
            # them -- on repeated runs in the same process this shows up as
            # each run getting slower than the last. Drop everything and
            # force a real reclaim before returning.
            artifacts = None
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            gc.collect()


register(CrossAttentionBackend())
