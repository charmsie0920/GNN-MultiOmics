"""Interface every swappable model-training backend implements."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable


class ModelBackend(ABC):
    """A single experiment script wired up to run from the backend.

    Implementations should not modify the wrapped script's own CLI/training
    behavior; instead, load and invoke it (e.g. via `importlib`) from
    `run()`, adding thin post-training inference on top.
    """

    name: str
    expected_duration_seconds: dict[str, float]

    @abstractmethod
    def run(
        self,
        log: Callable[[str], None],
        target_cell_line: str,
        on_results: Callable[[list[dict]], None],
        on_progress: Callable[[int, int], None],
    ) -> None:
        """Execute the backend, forwarding relevant stdout lines to `log`.

        Call `on_progress(current_epoch, max_epochs)` as training progresses
        (best-effort — as often as the wrapped script's own logging makes
        available) so the UI can show real progress instead of a fixed
        time estimate. After training, run inference for `target_cell_line`
        across every drug the run was trained on and pass the ranked
        results to `on_results`.
        """
