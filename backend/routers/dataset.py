from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
from fastapi import APIRouter, HTTPException, UploadFile

from backend.dataset_schema import REQUIRED_COLUMNS, UPLOADED_CSV_PATH, delete_uploaded_dataset
from backend.schemas.dataset import DatasetUploadResponse

router = APIRouter(prefix="/api/v1/dataset", tags=["dataset"])

_REPO_ROOT = Path(__file__).resolve().parents[2]

# The final model reads each cell line's row of the per-modality omics PCA
# files, so only cell lines in them can be scored. The three files share one
# cell-line index (asserted by `load_omics_subset`), and every one of those
# cell lines is also a node of the graph the pair features come from.
_OMICS_REFERENCE_CSV = _REPO_ROOT / "data" / "processed" / "transcriptomics_pca.csv"

# Read once per process rather than on every upload. `None` means "not loaded
# yet"; an empty set is a legitimate cached result.
_reference_ids: set[str] | None = None


@router.post("/upload", response_model=DatasetUploadResponse)
async def upload_dataset(file: UploadFile) -> DatasetUploadResponse:
    if not file.filename or not file.filename.lower().endswith(".csv"):
        raise HTTPException(status_code=400, detail="Uploaded file must be a .csv file.")

    raw_bytes = await file.read()
    try:
        frame = pd.read_csv(io.BytesIO(raw_bytes))
    except Exception as exc:  # pandas raises various parser errors
        raise HTTPException(status_code=400, detail=f"Could not parse CSV: {exc}") from exc

    missing_columns = [column for column in REQUIRED_COLUMNS if column not in frame.columns]
    if missing_columns:
        raise HTTPException(
            status_code=400,
            detail=f"Uploaded CSV is missing required columns: {', '.join(missing_columns)}",
        )

    UPLOADED_CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    UPLOADED_CSV_PATH.write_bytes(raw_bytes)

    cell_line_ids = _valid_target_cell_lines(frame)

    return DatasetUploadResponse(
        filename=file.filename,
        saved_path=str(UPLOADED_CSV_PATH),
        row_count=len(frame),
        columns=list(frame.columns),
        cell_line_ids=cell_line_ids,
        message="Dataset uploaded and validated successfully.",
    )


@router.delete("/upload")
def delete_upload() -> dict[str, bool]:
    """Remove the uploaded dataset. The desktop app calls this as it closes."""
    return {"deleted": delete_uploaded_dataset()}


def _reference_cell_ids() -> set[str]:
    """Cell-line ids the active backend can score: those with an omics profile."""
    global _reference_ids

    if _reference_ids is None:
        if not _OMICS_REFERENCE_CSV.exists():
            return set()
        _reference_ids = set(pd.read_csv(_OMICS_REFERENCE_CSV, index_col=0, usecols=[0]).index.astype(str))
    return _reference_ids


def _valid_target_cell_lines(frame: pd.DataFrame) -> list[str]:
    """Cell lines the model can actually predict for: present in both the
    upload and the backend's reference set (see `_reference_cell_ids`)."""
    reference_ids = _reference_cell_ids()
    if not reference_ids:
        return []
    uploaded_ids = frame["sanger_model_id"].astype(str).unique()
    return sorted(cell_id for cell_id in uploaded_ids if cell_id in reference_ids)
