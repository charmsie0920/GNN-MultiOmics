from __future__ import annotations

from fastapi import APIRouter, HTTPException

from backend.model_backends.registry import get_backend
from backend.schemas.results import (
    DrugResult,
    EnrichmentResponse,
    EnrichmentTerm,
    GeneAttribution,
    GeneAttributionResponse,
    TargetRecovery,
    TrainingHistoryPoint,
)
from backend.services.run_manager import run_manager
from src.interpretation.enrichment import enrich
from src.interpretation.genes import load_symbol_map, target_recovery

router = APIRouter(prefix="/api/v1/results", tags=["results"])


@router.get("/drug-ranking/{run_id}", response_model=list[DrugResult])
def drug_ranking(run_id: str) -> list[DrugResult]:
    state = run_manager.get(run_id)
    if state is None:
        raise HTTPException(status_code=404, detail=f"Unknown run_id: {run_id}")
    if state.status != "completed" or state.results is None:
        raise HTTPException(status_code=404, detail=f"Run {run_id} has no results yet (status={state.status}).")
    return [DrugResult(**row) for row in state.results]


@router.get("/training-history/{run_id}", response_model=list[TrainingHistoryPoint])
def training_history(run_id: str) -> list[TrainingHistoryPoint]:
    state = run_manager.get(run_id)
    if state is None:
        raise HTTPException(status_code=404, detail=f"Unknown run_id: {run_id}")
    if state.status != "completed" or state.training_history is None:
        raise HTTPException(status_code=404, detail=f"Run {run_id} has no training history yet (status={state.status}).")
    return [TrainingHistoryPoint(**row) for row in state.training_history]


def _completed_state(run_id: str):
    """Fetch a run that is finished and has results, or raise the right 404."""
    state = run_manager.get(run_id)
    if state is None:
        raise HTTPException(status_code=404, detail=f"Unknown run_id: {run_id}")
    if state.status != "completed" or state.results is None:
        raise HTTPException(status_code=404, detail=f"Run {run_id} has no results yet (status={state.status}).")
    return state


def _attribute(run_id: str, drug_id: str, top_k: int) -> tuple[list[dict], dict, object]:
    """Shared front half of both interpretation endpoints.

    Dispatches through the registry rather than importing any particular
    backend, so a model that cannot explain is refused here by its own
    declaration instead of by a crash deeper down.
    """
    state = _completed_state(run_id)

    backend = get_backend(state.backend_name)
    if not backend.supports_interpretation:
        raise HTTPException(
            status_code=501,
            detail=f"Backend {state.backend_name!r} does not support gene-level interpretation.",
        )

    drug_row = next((row for row in state.results if str(row.get("drug_id")) == str(drug_id)), None)
    if drug_row is None:
        raise HTTPException(status_code=404, detail=f"Run {run_id} has no result for drug_id {drug_id!r}.")

    try:
        genes = backend.explain(run_id, state.target_cell_line, drug_id, top_k)
    except RuntimeError as exc:
        # Raised when the run's weights cannot be faithfully rebuilt.
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    return genes, drug_row, state


@router.get("/gene-attribution/{run_id}", response_model=GeneAttributionResponse)
def gene_attribution(run_id: str, drug_id: str, top_k: int = 50) -> GeneAttributionResponse:
    genes, drug_row, state = _attribute(run_id, drug_id, top_k)

    symbol_map = load_symbol_map()
    known_symbols = {symbol.upper() for symbol in symbol_map.values() if symbol}
    recovery = target_recovery(genes, drug_row.get("putative_target", ""), known_symbols)

    return GeneAttributionResponse(
        drug_id=str(drug_id),
        target_cell_line=state.target_cell_line,
        genes=[GeneAttribution(**gene) for gene in genes],
        target_recovery=TargetRecovery(**recovery),
    )


@router.get("/enrichment/{run_id}", response_model=EnrichmentResponse)
def enrichment(run_id: str, drug_id: str, top_k: int = 50) -> EnrichmentResponse:
    genes, _drug_row, _state = _attribute(run_id, drug_id, top_k)

    result = enrich([gene["gene_symbol"] for gene in genes])
    return EnrichmentResponse(
        status=result["status"],
        message=result["message"],
        terms=[EnrichmentTerm(**term) for term in result["terms"]],
    )
