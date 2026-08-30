from __future__ import annotations

from fastapi import APIRouter, HTTPException

from backend.schemas.results import DrugResult
from backend.services.run_manager import run_manager

router = APIRouter(prefix="/api/v1/results", tags=["results"])


@router.get("/drug-ranking/{run_id}", response_model=list[DrugResult])
def drug_ranking(run_id: str) -> list[DrugResult]:
    state = run_manager.get(run_id)
    if state is None:
        raise HTTPException(status_code=404, detail=f"Unknown run_id: {run_id}")
    if state.status != "completed" or state.results is None:
        raise HTTPException(status_code=404, detail=f"Run {run_id} has no results yet (status={state.status}).")
    return [DrugResult(**row) for row in state.results]
