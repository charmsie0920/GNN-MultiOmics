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

# The heterogeneous graph the active backend predicts over: only cell lines
# that exist as `cell_line` nodes in it can be scored, so it is the real
# source of truth for which uploaded cell lines are selectable.
_GRAPH_PATH = _REPO_ROOT / "src" / "graph" / "hetero_graph.pt"

# Fallback for when the graph file is missing: one of the fixed per-modality
# omics reference files (any of the three shares the same cell-line index).
_OMICS_REFERENCE_CSV = _REPO_ROOT / "data" / "processed" / "transcriptomics_pca.csv"

# The graph is ~21 MB, so the node id list is read once per process rather
# than on every upload. `None` means "not loaded yet"; an empty set is a
# legitimate cached result.
_graph_cell_ids: set[str] | None = None


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

    cell_line_ids = _valid_target_cell_lines(frame)

    return DatasetUploadResponse(
        filename=file.filename,
        saved_path=str(TARGET_CSV_PATH),
        row_count=len(frame),
        columns=list(frame.columns),
        cell_line_ids=cell_line_ids,
        message="Dataset uploaded and validated successfully.",
    )


def _reference_cell_ids() -> set[str]:
    """Cell-line ids the active backend can score, preferring the graph's own
    `cell_line` node ids and degrading to the omics reference CSV if the graph
    file is unavailable."""
    global _graph_cell_ids

    if _graph_cell_ids is None and _GRAPH_PATH.exists():
        try:
            import torch

            graph = torch.load(_GRAPH_PATH, weights_only=False)
            _graph_cell_ids = {str(cell_id) for cell_id in graph["cell_line"].node_ids}
        except Exception:  # noqa: BLE001 - fall back rather than fail the upload
            _graph_cell_ids = set()

    if _graph_cell_ids:
        return _graph_cell_ids

    if not _OMICS_REFERENCE_CSV.exists():
        return set()
    return set(pd.read_csv(_OMICS_REFERENCE_CSV, index_col=0, usecols=[0]).index.astype(str))


def _valid_target_cell_lines(frame: pd.DataFrame) -> list[str]:
    """Cell lines the model can actually predict for: present in both the
    upload and the backend's reference set (see `_reference_cell_ids`)."""
    reference_ids = _reference_cell_ids()
    if not reference_ids:
        return []
    uploaded_ids = frame["sanger_model_id"].astype(str).unique()
    return sorted(cell_id for cell_id in uploaded_ids if cell_id in reference_ids)
