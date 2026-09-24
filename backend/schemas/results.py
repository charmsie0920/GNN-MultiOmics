from __future__ import annotations

from pydantic import BaseModel


class DrugResult(BaseModel):
    # GDSC drug names are not unique, so the id is what identifies a row for
    # anything acting on a specific drug (interpretation, in particular).
    drug_id: str
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
    apart by `epoch` being set or not.

    Every point, either shape, also carries `eval_val_rmse`/`eval_val_pcc`:
    the backend's own validation metrics for the final model (see
    `model_backends._common.attach_eval_metrics`), for the UI to report as-is.
    """

    epoch: int | None = None
    train_loss: float | None = None
    val_rmse: float | None = None
    val_pcc: float | None = None
    actual_ln_ic50: float | None = None
    predicted_ln_ic50: float | None = None
    eval_val_rmse: float | None = None
    eval_val_pcc: float | None = None


class GeneAttribution(BaseModel):
    """One gene's contribution to a single prediction.

    `score` is signed against ln(IC50): negative pushed the prediction down,
    toward sensitivity. The two evidence flags say whether the gene was reached
    through real annotation (a cancer-driver mutation in this cell line, or this
    drug's curated target) rather than inferred through the PPI network.
    """

    gene_symbol: str
    protein_id: str
    score: float
    direction: str
    is_driver_mutation: bool
    is_drug_target: bool


class RecoveredTarget(BaseModel):
    gene_symbol: str
    rank: int


class TargetRecovery(BaseModel):
    """Whether the drug's GDSC-annotated target surfaced in the top genes.

    `unmatched_tokens` matters for honesty: `putative_target` mixes real symbols
    with free-text mechanisms ("Microtubule destabiliser") and informal names
    ("MEK1"), and those can never be recovered. Reporting them separately stops
    a naming mismatch being read as a model failure.
    """

    putative_target: str
    target_genes: list[str]
    unmatched_tokens: list[str]
    recovered: list[RecoveredTarget]
    checked: bool


class GeneAttributionResponse(BaseModel):
    drug_id: str
    target_cell_line: str
    genes: list[GeneAttribution]
    target_recovery: TargetRecovery


class EnrichmentTerm(BaseModel):
    term: str
    library: str
    p_value: float
    adjusted_p_value: float
    combined_score: float
    overlap: str
    genes: list[str]


class EnrichmentResponse(BaseModel):
    """Enriched terms plus why the list looks the way it does.

    The status is part of the contract: "nothing was significant" and "Enrichr
    was unreachable" are indistinguishable from an empty term list alone, and
    the UI needs to say which happened.
    """

    status: str
    message: str
    terms: list[EnrichmentTerm]
