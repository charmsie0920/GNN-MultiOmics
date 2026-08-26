from __future__ import annotations

from pydantic import BaseModel


class StartRunResponse(BaseModel):
    run_id: str
    backend: str
    device: str
    expected_duration_seconds: float
    status: str


class RunStatusResponse(BaseModel):
    run_id: str
    status: str
    progress_percent: float
    elapsed_seconds: float
    expected_duration_seconds: float
    new_log_lines: list[str]
    next_since: int
    error_message: str | None = None
