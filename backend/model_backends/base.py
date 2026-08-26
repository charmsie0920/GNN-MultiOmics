"""Interface every swappable model-training backend implements."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable


class ModelBackend(ABC):
    """A single experiment script wired up to run from the backend.

    Implementations should not modify the wrapped script; instead, load and
    invoke it (e.g. via `importlib`) from `run()`.
    """

    name: str
    expected_duration_seconds: dict[str, float]

    @abstractmethod
    def run(self, log: Callable[[str], None]) -> None:
        """Execute the backend, forwarding relevant stdout lines to `log`."""
