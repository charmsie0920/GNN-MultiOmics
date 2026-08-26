"""Tracks background model-training runs started via the API."""

from __future__ import annotations

import random
import threading
import time
import traceback
import uuid
from dataclasses import dataclass, field

# Random range for each "still running, past the estimate" time extension,
# so the countdown doesn't visibly tick by the same amount every time.
_OVERDUE_EXTENSION_RANGE_SECONDS = (20.0, 45.0)

# Once a run is confirmed overdue, the displayed overall-progress percent is
# pinned here rather than recomputed (which would otherwise wobble slightly
# each time the expected duration is extended).
_OVERDUE_PROGRESS_PERCENT = 99.0

import torch

from backend.model_backends.registry import DEFAULT_BACKEND, get_backend


@dataclass
class RunState:
    run_id: str
    backend_name: str
    device: str
    expected_duration_seconds: float
    status: str = "running"  # "running" | "completed" | "failed"
    overdue: bool = False
    start_time: float = field(default_factory=time.monotonic)
    log_lines: list[str] = field(default_factory=list)
    error_message: str | None = None
    lock: threading.Lock = field(default_factory=threading.Lock)

    def append_log(self, line: str) -> None:
        with self.lock:
            self.log_lines.append(line)

    def log_lines_since(self, since: int) -> tuple[list[str], int]:
        with self.lock:
            return self.log_lines[since:], len(self.log_lines)

    def extend_if_overdue(self, elapsed: float) -> None:
        """Push the expected duration out if a run is taking longer than estimated.

        Once a run first crosses this point it's marked `overdue`, which
        pins its displayed progress percent (see `progress_percent`) so only
        the estimated time keeps moving from here — each extension is a
        random amount rather than a fixed one, so the countdown doesn't look
        mechanically identical every time it happens.
        """
        if elapsed >= self.expected_duration_seconds * 0.98:
            self.overdue = True
            extension = random.uniform(*_OVERDUE_EXTENSION_RANGE_SECONDS)
            self.expected_duration_seconds = elapsed + extension

    def progress_percent(self, elapsed: float) -> float:
        if self.status != "running":
            return 100.0
        if self.overdue:
            return _OVERDUE_PROGRESS_PERCENT
        return min(_OVERDUE_PROGRESS_PERCENT, (elapsed / self.expected_duration_seconds) * 100.0)


class RunManager:
    def __init__(self) -> None:
        self._runs: dict[str, RunState] = {}

    def start_run(self, backend_name: str = DEFAULT_BACKEND) -> RunState:
        backend = get_backend(backend_name)
        device = "cuda" if torch.cuda.is_available() else "cpu"
        state = RunState(
            run_id=uuid.uuid4().hex,
            backend_name=backend_name,
            device=device,
            expected_duration_seconds=backend.expected_duration_seconds[device],
        )
        self._runs[state.run_id] = state

        thread = threading.Thread(target=self._execute, args=(state, backend), daemon=True)
        thread.start()
        return state

    def get(self, run_id: str) -> RunState | None:
        return self._runs.get(run_id)

    @staticmethod
    def _execute(state: RunState, backend) -> None:
        try:
            backend.run(log=state.append_log)
        except Exception as exc:  # noqa: BLE001 - surfaced to the UI, full traceback stays in the terminal
            traceback.print_exc()
            state.error_message = str(exc)
            state.status = "failed"
        else:
            state.status = "completed"


run_manager = RunManager()
