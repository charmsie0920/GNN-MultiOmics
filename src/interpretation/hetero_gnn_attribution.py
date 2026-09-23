"""Per-protein attribution for `HeteroIC50GNN`.

The only module that knows this model's internals. The checkpoint's convs are
`SAGEConv`, so there are no attention weights to read; instead the predicted
ln(IC50) is differentiated with respect to the protein embedding table and each
protein scored by gradient x embedding -- a first-order estimate of how much
that protein's representation contributed to this particular prediction.

Protein features on the graph are all zeros (`03_graph_construction.py`
zero-fills them), so `protein_embedding.weight` is where every protein's
learned identity actually lives and is the only meaningful attribution target.

With two message-passing layers only proteins within two hops of the cell-line
or drug node can influence the output, so the gradient is zero everywhere else
and the candidate set stays sparse without any thresholding.
"""

from __future__ import annotations

import torch
from torch import nn

_MUTATION_EDGE = ("cell_line", "has_mutation", "protein")
_TARGET_EDGE = ("drug", "targets", "protein")


def attribute_pair(
    model: nn.Module,
    x_dict: dict,
    edge_index_dict: dict,
    cell_idx: int,
    drug_idx: int,
) -> torch.Tensor:
    """Signed per-protein contribution to one (cell line, drug) prediction.

    Returns a tensor of shape [num_proteins]. Negative means the protein pushed
    ln(IC50) down (toward sensitivity), positive means up (toward resistance).

    Runs in full eval mode -- unlike the MC-dropout prediction path, dropout
    must be off and BatchNorm must use its running statistics, or repeated calls
    on the same pair would disagree with each other and with the reported
    prediction.
    """
    nn.Module.eval(model)

    embedding = model.protein_embedding.weight
    embedding.requires_grad_(True)

    device = embedding.device
    cell_index = torch.tensor([cell_idx], dtype=torch.long, device=device)
    drug_index = torch.tensor([drug_idx], dtype=torch.long, device=device)

    prediction = model(x_dict, edge_index_dict, cell_index, drug_index).squeeze()
    gradient = torch.autograd.grad(prediction, embedding)[0]

    return (gradient * embedding).sum(dim=-1).detach()


def evidence_sets(
    edge_index_dict: dict,
    cell_idx: int,
    drug_idx: int,
    node_ids: list[str],
) -> tuple[set[str], set[str]]:
    """Protein ids this cell line mutates and this drug targets.

    Both edge sets are real biological annotation rather than anything learned:
    `has_mutation` is restricted to cancer-driver mutations and `targets` comes
    from GDSC's curated drug-target column. Surfacing them alongside the scores
    is what lets the UI distinguish a gene the model reached via the PPI network
    from one it was handed directly.
    """
    return (
        _neighbours(edge_index_dict, _MUTATION_EDGE, cell_idx, node_ids),
        _neighbours(edge_index_dict, _TARGET_EDGE, drug_idx, node_ids),
    )


def _neighbours(edge_index_dict: dict, edge_type: tuple, source_idx: int, node_ids: list[str]) -> set[str]:
    edge_index = edge_index_dict.get(edge_type)
    if edge_index is None:
        return set()
    source, destination = edge_index
    matched = destination[source == source_idx].tolist()
    return {node_ids[i] for i in matched}
