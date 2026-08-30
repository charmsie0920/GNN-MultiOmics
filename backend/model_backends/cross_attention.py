"""Backend wrapper for experiments/Cross Attention Fusion/cross_attention_baseline.py.

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
from torch import nn

from backend.model_backends.base import ModelBackend
from backend.model_backends.registry import register

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT_PATH = _REPO_ROOT / "experiments" / "Cross Attention Fusion" / "cross_attention_baseline.py"

# MC-dropout forward passes per (cell-line, drug) pair used to derive the
# confidence score — see _predict_drug_panel below.
_MC_SAMPLES = 30

# Fixed IC50 (uM) cutoffs for the sensitivity ranking, standard
# pharmacological convention for how potent/sensitive a response is:
# sub-micromolar is considered a strong (HIGH SENSITIVITY) response, single-
# to-low-double-digit micromolar is MEDIUM, anything above is LOW. Unlike a
# percentile/tertile split, these don't shift based on what else happens to
# be in the candidate drug list for a given run.
_HIGH_SENSITIVITY_MAX_UM = 1.0
_MEDIUM_SENSITIVITY_MAX_UM = 10.0


def _sensitivity_ranking(ic50_um: float) -> str:
    if ic50_um < _HIGH_SENSITIVITY_MAX_UM:
        return "HIGH SENSITIVITY"
    if ic50_um < _MEDIUM_SENSITIVITY_MAX_UM:
        return "MEDIUM"
    return "LOW"


def _rank_predictions(raw: list[dict], val_rmse: float) -> list[dict]:
    """Turn raw per-drug {drug_id, ln_ic50_mean, ln_ic50_std} into display-ready results.

    - predicted_ic50_um: back-transform from ln(IC50 uM), the model's native scale.
    - confidence_percent: compares each drug's MC-dropout std against `val_rmse`
      -- the model's own measured error on held-out validation data -- via
      `100 * exp(-std / val_rmse)`. This is an *absolute* scale (comparable
      across separate runs), unlike a min-max normalization over the current
      run's candidate set, which always forces some drug to exactly 0% and
      another to exactly 100% regardless of whether the underlying spread is
      actually large or small.
    - ranking: fixed IC50 (uM) thresholds, see `_sensitivity_ranking`.
    """
    if not raw:
        return []

    ic50_um = np.array([np.exp(r["ln_ic50_mean"]) for r in raw])
    stds = np.array([r["ln_ic50_std"] for r in raw])

    # Guard against a degenerate (near-zero) val_rmse blowing up the ratio.
    scale = max(val_rmse, 1e-6)
    confidence = 100.0 * np.exp(-stds / scale)

    return [
        {
            "drug_name": raw[i]["drug_id"],
            "predicted_ic50_um": float(ic50_um[i]),
            "ranking": _sensitivity_ranking(float(ic50_um[i])),
            "confidence_percent": float(confidence[i]),
        }
        for i in range(len(raw))
    ]


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
    # calibrating the confidence score against -- see _rank_predictions.
    val_pred = module.predict(model, omics_val, drug_onehot[val_idx], keys, device, module.BATCH_SIZE)
    val_rmse = float(np.sqrt(np.mean((val_pred - y[val_idx]) ** 2)))

    return {
        "model": model,
        "keys": keys,
        "gathered": gathered,
        "cell_ids": cell_ids,
        "drug_levels": drug_levels,
        "n_drugs": n_drugs,
        "val_rmse": val_rmse,
        "device": device,
    }


def _enable_mc_dropout(model: nn.Module) -> None:
    """Put `model` in true eval mode, then re-enable train-mode only on
    Dropout layers, so BatchNorm keeps using its running stats (safe for
    any batch size) while dropout still injects stochasticity for MC
    sampling.

    Must call `nn.Module.eval(model)` (the unbound class method), not
    `model.eval()` -- this function is installed *as* `model.eval` itself
    (see `_predict_drug_panel`), so calling the bound method here would
    recurse into this same function forever.
    """
    nn.Module.eval(model)
    for submodule in model.modules():
        if isinstance(submodule, nn.Dropout):
            submodule.train()


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
    model.eval = lambda: _enable_mc_dropout(model)
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
    ) -> None:
        module = _load_module()
        on_progress(0, module.MAX_EPOCHS)

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

        # Every function in `module` resolves the bare name `print` against
        # the module's own globals, so this alone routes all of the
        # script's existing progress prints (load/split/epoch logs) to the
        # UI too, with no changes to the script itself.
        module.print = tee_print

        try:
            artifacts = _prepare_and_train(module)
            log(f"[confidence] validation RMSE (ln IC50) = {artifacts['val_rmse']:.4f} -- used as the confidence scale")
            log(f"Running inference for target cell line: {target_cell_line}")
            raw_predictions = _predict_drug_panel(module, artifacts, target_cell_line, mc_samples=_MC_SAMPLES)
            on_results(_rank_predictions(raw_predictions, artifacts["val_rmse"]))
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
