"""Model-agnostic helpers shared by every `ModelBackend` implementation.

These turn a backend's raw per-drug output into the exact shape the UI
expects. Keeping them here rather than duplicating them per backend means
every model reports sensitivity and confidence on the same scale, so two
backends' results pages are directly comparable -- and the ranking strings
stay in sync with the badge colors keyed off them in
`pages/final_results_page.py`.
"""

from __future__ import annotations

import numpy as np
from torch import nn

# Fixed IC50 (uM) cutoffs for the sensitivity ranking, standard
# pharmacological convention for how potent/sensitive a response is:
# sub-micromolar is considered a strong (HIGH SENSITIVITY) response, single-
# to-low-double-digit micromolar is MEDIUM, anything above is LOW. Unlike a
# percentile/tertile split, these don't shift based on what else happens to
# be in the candidate drug list for a given run.
HIGH_SENSITIVITY_MAX_UM = 1.0
MEDIUM_SENSITIVITY_MAX_UM = 10.0


def sensitivity_ranking(ic50_um: float) -> str:
    if ic50_um < HIGH_SENSITIVITY_MAX_UM:
        return "HIGH SENSITIVITY"
    if ic50_um < MEDIUM_SENSITIVITY_MAX_UM:
        return "MEDIUM"
    return "LOW"


def rank_predictions(raw: list[dict], val_rmse: float, drug_info: dict[str, dict]) -> list[dict]:
    """Turn raw per-drug {drug_id, ln_ic50_mean, ln_ic50_std} into display-ready results.

    - predicted_ic50_um: back-transform from ln(IC50 uM), the model's native scale.
    - confidence_percent: compares each drug's MC-dropout std against `val_rmse`
      -- the model's own measured error on held-out validation data -- via
      `100 * exp(-std / val_rmse)`. This is an *absolute* scale (comparable
      across separate runs), unlike a min-max normalization over the current
      run's candidate set, which always forces some drug to exactly 0% and
      another to exactly 100% regardless of whether the underlying spread is
      actually large or small.
    - ranking: fixed IC50 (uM) thresholds, see `sensitivity_ranking`.
    - drug_name/putative_target/pathway_name: looked up from `drug_info`,
      keyed by the same raw `drug_id` string the backend produced.
    """
    if not raw:
        return []

    ic50_um = np.array([np.exp(r["ln_ic50_mean"]) for r in raw])
    stds = np.array([r["ln_ic50_std"] for r in raw])

    # Guard against a degenerate (near-zero) val_rmse blowing up the ratio.
    scale = max(val_rmse, 1e-6)
    confidence = 100.0 * np.exp(-stds / scale)

    results = []
    for i in range(len(raw)):
        info = drug_info.get(raw[i]["drug_id"], {})
        results.append(
            {
                "drug_name": info.get("drug_name") or raw[i]["drug_id"],
                "putative_target": info.get("putative_target", ""),
                "pathway_name": info.get("pathway_name", ""),
                "predicted_ic50_um": float(ic50_um[i]),
                "ranking": sensitivity_ranking(float(ic50_um[i])),
                "confidence_percent": float(confidence[i]),
            }
        )
    return results


def enable_mc_dropout(model: nn.Module) -> None:
    """Put `model` in true eval mode, then re-enable train-mode only on
    Dropout layers, so BatchNorm keeps using its running stats (safe for
    any batch size) while dropout still injects stochasticity for MC
    sampling.

    Must call `nn.Module.eval(model)` (the unbound class method), not
    `model.eval()` -- in the cross-attention backend this function is
    installed *as* `model.eval` itself, so calling the bound method here
    would recurse into this same function forever.
    """
    nn.Module.eval(model)
    for submodule in model.modules():
        if isinstance(submodule, nn.Dropout):
            submodule.train()
