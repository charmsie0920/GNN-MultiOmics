"""Shared target-CSV column contract used to validate dataset uploads.

Mirrors the `sanger_model_id` / `drug_id` / `ln_ic50` columns that
`cross_attention_baseline.py` and its sibling experiment scripts under
`experiments/` read from `data/processed/aligned/gdsc2_response_master.csv`.
Kept independent of any specific experiment script so the backend never
needs to import from `experiments/` just to validate an upload.
"""

from __future__ import annotations

from pathlib import Path

REQUIRED_COLUMNS = ["sanger_model_id", "drug_id", "ln_ic50"]

# Where an uploaded dataset is kept for the lifetime of the application. It is
# a copy of its own, never the experiments' master CSV, because it is deleted
# when the application closes (see `delete_uploaded_dataset`). The folder sits
# under data/processed/, which is gitignored.
UPLOADED_CSV_PATH = (
    Path(__file__).resolve().parents[1] / "data" / "processed" / "uploads" / "uploaded_response.csv"
)


def delete_uploaded_dataset() -> bool:
    """Delete the uploaded dataset, if there is one.

    Returns whether a file was removed. A file that is still being read (a run
    in progress on Windows) cannot be removed and is left for the next attempt.
    """
    try:
        UPLOADED_CSV_PATH.unlink()
    except (FileNotFoundError, PermissionError):
        return False
    return True
