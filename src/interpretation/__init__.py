"""Biological interpretation of model predictions.

`genes.py` is model-agnostic: it turns a per-protein score vector into ranked,
annotated gene rows. `hetero_gnn_attribution.py` produces that score vector for
the heterogeneous GNN specifically. Keeping the two apart is what lets a future
model supply its own scorer without touching the API, the enrichment step, or
the UI.
"""
