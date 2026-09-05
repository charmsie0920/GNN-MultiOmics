"""Swappable model-training backends. Importing this package registers all of them."""

# Only the hetero GNN is registered -- it is the better model (test RMSE
# 1.3302 vs the cross-attention arm) and it ships with trained weights, so a
# run does not have to retrain from scratch. `cross_attention.py` is still on
# disk and can be brought back by adding its import here.
from backend.model_backends import hetero_gnn  # noqa: F401
