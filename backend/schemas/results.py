from __future__ import annotations

from pydantic import BaseModel


class DrugResult(BaseModel):
    drug_name: str
    putative_target: str
    pathway_name: str
    predicted_ic50_um: float
    ranking: str
    confidence_percent: float


class TrainingHistoryPoint(BaseModel):
    """One point for the results page's training-curve panel.

    Two disjoint shapes share this list: a real per-epoch metric point (when
    the backend actually trained -- `epoch`/`train_loss`/`val_rmse`/`val_pcc`
    set, `actual_ln_ic50`/`predicted_ln_ic50` left None), or a held-out
    validation pair (when the backend loaded a pretrained checkpoint instead
    -- no epochs were run, so `actual_ln_ic50`/`predicted_ln_ic50` are set
    instead so the panel still has real data to plot). The client tells them
    apart by which fields are non-None.
    """

    epoch: int | None = None
    train_loss: float | None = None
    val_rmse: float | None = None
    val_pcc: float | None = None
    actual_ln_ic50: float | None = None
    predicted_ln_ic50: float | None = None
