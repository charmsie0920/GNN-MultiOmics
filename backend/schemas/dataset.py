from __future__ import annotations

from pydantic import BaseModel


class DatasetUploadResponse(BaseModel):
    filename: str
    saved_path: str
    row_count: int
    columns: list[str]
    cell_line_ids: list[str]
    message: str
