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
    epoch: int
    train_loss: float
    val_rmse: float
    val_pcc: float
