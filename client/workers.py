"""QThread workers that run blocking API calls off the Qt UI thread."""

from __future__ import annotations

import time

from PySide6.QtCore import QThread, Signal

from client.api_client import ApiError, get_drug_ranking, get_run_status, start_run, upload_dataset_csv

POLL_INTERVAL_SECONDS = 1.0


class UploadWorker(QThread):
    """Uploads a dataset CSV off the UI thread.

    Only uploads — starting a run needs a target cell line, chosen by the
    user from `cell_line_ids` in the result, so that's a separate step via
    `StartRunWorker`.
    """

    succeeded = Signal(dict)
    failed = Signal(str)

    def __init__(self, file_path: str, parent=None) -> None:
        super().__init__(parent)
        self._file_path = file_path

    def run(self) -> None:  # noqa: N802
        try:
            upload_result = upload_dataset_csv(self._file_path)
        except ApiError as exc:
            self.failed.emit(str(exc))
        else:
            self.succeeded.emit(upload_result)


class StartRunWorker(QThread):
    """Starts a model run for a chosen target cell line, off the UI thread."""

    succeeded = Signal(dict)
    failed = Signal(str)

    def __init__(self, target_cell_line: str, parent=None) -> None:
        super().__init__(parent)
        self._target_cell_line = target_cell_line

    def run(self) -> None:  # noqa: N802
        try:
            run_result = start_run(self._target_cell_line)
        except ApiError as exc:
            self.failed.emit(str(exc))
        else:
            self.succeeded.emit(run_result)


class ResultsWorker(QThread):
    """Fetches the ranked drug predictions for a completed run, off the UI thread."""

    succeeded = Signal(list)
    failed = Signal(str)

    def __init__(self, run_id: str, parent=None) -> None:
        super().__init__(parent)
        self._run_id = run_id

    def run(self) -> None:  # noqa: N802
        try:
            results = get_drug_ranking(self._run_id)
        except ApiError as exc:
            self.failed.emit(str(exc))
        else:
            self.succeeded.emit(results)


class RunStatusPoller(QThread):
    """Polls a model run's status/log lines until it finishes."""

    statusUpdate = Signal(dict)

    def __init__(self, run_id: str, parent=None) -> None:
        super().__init__(parent)
        self._run_id = run_id
        self._stop = False

    def stop(self) -> None:
        self._stop = True

    def run(self) -> None:  # noqa: N802
        since = 0
        while not self._stop:
            try:
                payload = get_run_status(self._run_id, since)
            except ApiError as exc:
                self.statusUpdate.emit({
                    "run_id": self._run_id,
                    "status": "failed",
                    "error_message": str(exc),
                    "new_log_lines": [],
                    "progress_percent": 0.0,
                    "elapsed_seconds": 0.0,
                    "expected_duration_seconds": 0.0,
                })
                return

            since = payload["next_since"]
            self.statusUpdate.emit(payload)
            if payload["status"] != "running":
                return
            time.sleep(POLL_INTERVAL_SECONDS)
