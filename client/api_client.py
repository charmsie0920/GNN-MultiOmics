"""Thin synchronous HTTP client for the FastAPI backend.

Blocking by design — callers running on the Qt UI thread should use
`client.workers.UploadWorker` instead of calling these functions directly.
"""

from __future__ import annotations

from pathlib import Path

import requests

BASE_URL = "http://127.0.0.1:8000"


class ApiError(Exception):
    """Raised when a backend request fails, with a user-facing message."""


def upload_dataset_csv(file_path: str) -> dict:
    """POST a CSV file to the dataset upload endpoint.

    Returns the parsed JSON response on success (2xx).

    Raises:
        ApiError: On a network failure or a non-2xx response, with the
            backend's error detail message when available.
    """
    path = Path(file_path)
    try:
        with path.open("rb") as handle:
            response = requests.post(
                f"{BASE_URL}/api/v1/dataset/upload",
                files={"file": (path.name, handle, "text/csv")},
                timeout=30,
            )
    except requests.RequestException as exc:
        raise ApiError(f"Could not reach the backend: {exc}") from exc

    return _unwrap(response)


def start_run(target_cell_line: str, backend_name: str | None = None) -> dict:
    """Start a model run for the most recently uploaded dataset.

    Args:
        target_cell_line: The `sanger_model_id` to predict drug rankings for
            once training finishes.
        backend_name: Optional backend override (defaults to the server's
            configured default).

    Returns the parsed JSON response (`run_id`, `device`,
    `expected_duration_seconds`, ...) on success.

    Raises:
        ApiError: On a network failure or a non-2xx response.
    """
    params = {"target_cell_line": target_cell_line}
    if backend_name:
        params["backend_name"] = backend_name
    try:
        response = requests.post(f"{BASE_URL}/api/v1/model/run/start", params=params, timeout=10)
    except requests.RequestException as exc:
        raise ApiError(f"Could not reach the backend: {exc}") from exc
    return _unwrap(response)


def get_run_status(run_id: str, since: int = 0) -> dict:
    """Fetch the latest status/progress/log lines for a run.

    Raises:
        ApiError: On a network failure or a non-2xx response.
    """
    try:
        response = requests.get(
            f"{BASE_URL}/api/v1/model/run/status/{run_id}",
            params={"since": since},
            timeout=10,
        )
    except requests.RequestException as exc:
        raise ApiError(f"Could not reach the backend: {exc}") from exc
    return _unwrap(response)


def get_drug_ranking(run_id: str) -> list[dict]:
    """Fetch the ranked drug predictions for a completed run.

    Raises:
        ApiError: On a network failure or a non-2xx response (including a
            run that hasn't produced results yet).
    """
    try:
        response = requests.get(f"{BASE_URL}/api/v1/results/drug-ranking/{run_id}", timeout=10)
    except requests.RequestException as exc:
        raise ApiError(f"Could not reach the backend: {exc}") from exc
    return _unwrap(response)


def get_training_history(run_id: str) -> list[dict]:
    """Fetch the per-epoch training curve (train_loss/val_rmse/val_pcc) for a completed run.

    Raises:
        ApiError: On a network failure or a non-2xx response (including a
            run that hasn't produced training history yet).
    """
    try:
        response = requests.get(f"{BASE_URL}/api/v1/results/training-history/{run_id}", timeout=10)
    except requests.RequestException as exc:
        raise ApiError(f"Could not reach the backend: {exc}") from exc
    return _unwrap(response)


def _unwrap(response: requests.Response) -> dict:
    if not response.ok:
        detail = response.text
        try:
            detail = response.json().get("detail", detail)
        except ValueError:
            pass
        raise ApiError(str(detail))
    return response.json()
