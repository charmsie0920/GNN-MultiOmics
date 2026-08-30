from __future__ import annotations

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from backend.routers.dataset import router as dataset_router
from backend.routers.model_run import router as model_run_router
from backend.routers.results import router as results_router

app = FastAPI(title="MSC16 Placeholder API", version="0.1.0")

# Local desktop UI talking to a local backend — not public-facing.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(dataset_router)
app.include_router(model_run_router)
app.include_router(results_router)


class HealthResponse(BaseModel):
    status: str
    service: str


class DatasetInitResponse(BaseModel):
    run_id: str
    message: str


class ModelLogEntry(BaseModel):
    time: str
    level: str
    message: str


@app.get("/api/v1/health", response_model=HealthResponse)
def health() -> HealthResponse:
    # * Placeholder health check for local development and UI wiring.
    # * Replace with real dependency checks (database/model-store/queue) later.
    return HealthResponse(status="ok", service="msc16-fastapi")


@app.post("/api/v1/dataset/initialize", response_model=DatasetInitResponse)
def initialize_dataset() -> DatasetInitResponse:
    # * Placeholder dataset initialization response.
    # * Replace with actual ingestion + validation + job enqueue logic.
    return DatasetInitResponse(
        run_id="RUN-PLACEHOLDER-001",
        message="Dataset initialization accepted (placeholder).",
    )


@app.get("/api/v1/model/logs", response_model=list[ModelLogEntry])
def model_logs() -> list[ModelLogEntry]:
    # * Placeholder logs list.
    # * Replace with stream/log store retrieval from pipeline executor.
    return [
        ModelLogEntry(time="10:04:12", level="INFO", message="Initializing Preprocessing Engine v2.4.1"),
        ModelLogEntry(time="10:05:42", level="PROCESS", message="SEQ_001_A.fq.gz -> Normalizing Read Depth"),
        ModelLogEntry(time="10:12:30", level="SUCCESS", message="SEQ_001_A.fq.gz QC Passed"),
    ]


@app.get("/api/v1/model/visualization-state")
def model_visualization_state() -> dict[str, int | str]:
    # * Placeholder visualization status.
    # * Replace with graph execution state from running model job.
    return {
        "status": "computing",
        "progress_percent": 68,
        "elapsed_seconds": 252,
    }


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="127.0.0.1", port=8000)
