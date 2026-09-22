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

    # Whether this backend can explain an individual prediction in terms of
    # genes. Not every model can: a backend trained on PCA-projected omics has
    # no per-gene identity left to attribute to, so the capability is declared
    # rather than assumed, and the UI hides its interpretation panels when a
    # run reports False.
    supports_interpretation: bool = False

    def explain(
        self,
        run_id: str,
        target_cell_line: str,
        drug_id: str,
        top_k: int = 50,
    ) -> list[dict]:
        """Per-gene attribution for one (cell line, drug) prediction.

        Only backends that set `supports_interpretation = True` need to
        implement this. Each returned dict describes one gene: `gene_symbol`,
        `protein_id`, signed `score`, `direction`, `is_driver_mutation` and
        `is_drug_target`, ordered most important first.
        """
        raise NotImplementedError(f"{self.name} does not support interpretation.")

    @abstractmethod
    def run(
        self,
        log: Callable[[str], None],
        target_cell_line: str,
        on_results: Callable[[list[dict]], None],
        on_progress: Callable[[int, int], None],
        on_training_history: Callable[[list[dict]], None],
    ) -> None:
        """Execute the backend, forwarding relevant stdout lines to `log`.

        Call `on_progress(current_epoch, max_epochs)` as training progresses
        (best-effort — as often as the wrapped script's own logging makes
        available) so the UI can show real progress instead of a fixed
        time estimate. After training, run inference for `target_cell_line`
        across every drug the run was trained on and pass the ranked
        results to `on_results`. If the wrapped script logs per-epoch
        metrics, call `on_training_history` once with the full list of
        `{"epoch": int, "train_loss": float, "val_rmse": float, "val_pcc": float}`
        entries collected during training (empty list if none were parsed).
        """
