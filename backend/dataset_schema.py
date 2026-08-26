"""Shared target-CSV column contract used to validate dataset uploads.

Mirrors the `sanger_model_id` / `drug_id` / `ln_ic50` columns that
`cross_attention_baseline.py` and its sibling experiment scripts under
`experiments/` read from `data/processed/aligned/gdsc2_response_master.csv`.
Kept independent of any specific experiment script so the backend never
needs to import from `experiments/` just to validate an upload.
"""

from __future__ import annotations

REQUIRED_COLUMNS = ["sanger_model_id", "drug_id", "ln_ic50"]
