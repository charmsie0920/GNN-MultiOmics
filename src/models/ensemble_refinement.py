"""XGBoost residual-refinement stage, per MoGraphDRP §2.4.

The paper's final architectural component: rather than trusting the neural
model's output directly, its compressed interaction vector `f` and its initial
prediction `y_hat` are concatenated and handed to a gradient-boosted tree
ensemble, which learns to correct the residual errors the network leaves
behind. They report it as a substantial win (Table 7: RMSE 0.833 -> 0.674 on
their independent test set, ~19.7%).

This module applies the same idea to our models. The "interaction vector" here
is the penultimate hidden activation of the prediction head -- the closest
analogue to their 128-dim `f_ij`, being the last representation before the
scalar readout.

The refiner is fit on the *training* fold only and applied to val/test, so it
cannot see held-out labels; the cell-line-grouped split is inherited unchanged
from the base model's run.
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
from xgboost import XGBRegressor

# Hyperparameters as specified in MoGraphDRP §2.4 / Table 1.
N_ESTIMATORS = 100
MAX_DEPTH = 6
LEARNING_RATE = 0.05
SUBSAMPLE = 0.8
RANDOM_STATE = 42


def build_refiner_input(interaction: np.ndarray, prediction: np.ndarray) -> np.ndarray:
    """Concatenate the interaction vector with the base model's scalar prediction.

    Mirrors the paper's Z_ij = f_ij (+) y_hat_ij (their eq. 14).
    """
    return np.concatenate([interaction, prediction.reshape(-1, 1)], axis=1).astype(np.float32)


def fit_refiner(
    interaction_train: np.ndarray,
    prediction_train: np.ndarray,
    y_train: np.ndarray,
) -> XGBRegressor:
    """Fit the XGBoost refiner on the training fold's interaction vectors."""
    model = XGBRegressor(
        n_estimators=N_ESTIMATORS,
        max_depth=MAX_DEPTH,
        learning_rate=LEARNING_RATE,
        subsample=SUBSAMPLE,
        random_state=RANDOM_STATE,
        objective="reg:squarederror",
        n_jobs=-1,
    )
    model.fit(build_refiner_input(interaction_train, prediction_train), y_train)
    return model


def refine(model: XGBRegressor, interaction: np.ndarray, prediction: np.ndarray) -> np.ndarray:
    """Apply a fitted refiner, returning the corrected predictions."""
    return model.predict(build_refiner_input(interaction, prediction))


def feature_importance(model: XGBRegressor, top_n: int = 10) -> Dict[str, float]:
    """Top-N refiner feature importances, keyed f0..fN-1 plus 'predicted_ic50'.

    The paper's Fig 10 reports the base model's own prediction as one of the
    highest-importance features; this makes the same check available here.
    """
    scores = model.feature_importances_
    n_interaction = len(scores) - 1
    names = [f"f{i}" for i in range(n_interaction)] + ["predicted_ic50"]
    ranked = sorted(zip(names, scores), key=lambda kv: kv[1], reverse=True)
    return {name: float(score) for name, score in ranked[:top_n]}
