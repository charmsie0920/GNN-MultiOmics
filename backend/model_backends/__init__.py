"""Swappable model-training backends. Importing this package registers all of them."""

# Only the final model is registered -- it is the model the report presents
# (test RMSE 1.2737 +/- 0.0144 against 1.3092 for the hetero GNN, both under
# the grouped split) and its five trained seeds are committed, so a run never
# retrains. `hetero_gnn.py` and `cross_attention.py` are still on disk and can
# be brought back by adding their imports here.
from backend.model_backends import final_model  # noqa: F401
