"""QThread workers that run blocking API calls off the Qt UI thread."""

from __future__ import annotations

import time

from PySide6.QtCore import QThread, Signal

from client.api_client import ApiError, get_run_status, start_run, upload_dataset_csv

POLL_INTERVAL_SECONDS = 1.0


class UploadWorker(QThread):
    """Uploads a dataset CSV, then starts a model run, off the UI thread."""

    succeeded = Signal(dict)
    failed = Signal(str)

    def __init__(self, file_path: str, parent=None) -> None:
        super().__init__(parent)
        self._file_path = file_path

    def run(self) -> None:  # noqa: N802
        try:
            upload_result = upload_dataset_csv(self._file_path)
            run_result = start_run()
        except ApiError as exc:
            self.failed.emit(str(exc))
        else:
            self.succeeded.emit({**upload_result, **run_result})


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
