from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
from fastapi import APIRouter, HTTPException, UploadFile

from backend.dataset_schema import REQUIRED_COLUMNS
from backend.schemas.dataset import DatasetUploadResponse

router = APIRouter(prefix="/api/v1/dataset", tags=["dataset"])

# Canonical location every experiment script under experiments/ reads its
# target CSV from (e.g. cross_attention_baseline.py's TARGET_CSV).
_REPO_ROOT = Path(__file__).resolve().parents[2]
TARGET_CSV_PATH = _REPO_ROOT / "data" / "processed" / "aligned" / "gdsc2_response_master.csv"


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

    TARGET_CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    TARGET_CSV_PATH.write_bytes(raw_bytes)

    return DatasetUploadResponse(
        filename=file.filename,
        saved_path=str(TARGET_CSV_PATH),
        row_count=len(frame),
        columns=list(frame.columns),
        message="Dataset uploaded and validated successfully.",
    )
