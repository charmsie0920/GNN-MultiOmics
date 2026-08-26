from __future__ import annotations

import time

from fastapi import APIRouter, HTTPException

from backend.model_backends.registry import DEFAULT_BACKEND
from backend.schemas.model_run import RunStatusResponse, StartRunResponse
from backend.services.run_manager import run_manager

router = APIRouter(prefix="/api/v1/model", tags=["model-run"])


@router.post("/run/start", response_model=StartRunResponse)
def start_run(backend_name: str = DEFAULT_BACKEND) -> StartRunResponse:
    state = run_manager.start_run(backend_name)
    return StartRunResponse(
        run_id=state.run_id,
        backend=state.backend_name,
        device=state.device,
        expected_duration_seconds=state.expected_duration_seconds,
        status=state.status,
    )


@router.get("/run/status/{run_id}", response_model=RunStatusResponse)
def run_status(run_id: str, since: int = 0) -> RunStatusResponse:
    state = run_manager.get(run_id)
    if state is None:
        raise HTTPException(status_code=404, detail=f"Unknown run_id: {run_id}")

    elapsed = time.monotonic() - state.start_time
    if state.status == "running":
        state.extend_if_overdue(elapsed)
    progress_percent = state.progress_percent(elapsed)

    new_log_lines, next_since = state.log_lines_since(since)

    return RunStatusResponse(
        run_id=state.run_id,
        status=state.status,
        progress_percent=progress_percent,
        elapsed_seconds=elapsed,
        expected_duration_seconds=state.expected_duration_seconds,
        new_log_lines=new_log_lines,
        next_since=next_since,
        error_message=state.error_message,
    )
