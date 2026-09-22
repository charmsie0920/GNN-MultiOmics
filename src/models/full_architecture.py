"""The full proposed architecture: cross-attention fusion + PPI graph + molecular graph.

Every prior model implements a *part* of the architecture in
`docs/Architecture Simplified.jpeg`:

- `experiments/06_full_matrix/cross_attention_matrix.py` fuses omics with
  cross-attention but has no graph at all.
- `experiments/07_gnn_ablation/gnn_baseline.py` (E11/E13) does message passing
  over the PPI graph, but feeds it raw concatenated PCA blocks as `cell_line`
  features and a 2048-bit fingerprint as `drug` features -- no fusion, no
  molecular structure.
- `experiments/11_molecular_graph/molecular_graph_matrix.py` adds the
  atom-level drug graph, but drops the PPI graph.

This module composes all three into the architecture as proposed:

    GE / Mut_CNV / Proteomics (PCA-128 each)
            |
            +-- cross-attention fusion ------> cell_line node features (256)
    SMILES  +-- molecular graph GCN ---------> drug node features (128)
            |
            +-- HeteroConv message passing over the PPI graph
            |     (protein<->protein, drug->protein, cell_line->protein)
            |
            +-- cell_line embedding ++ drug embedding --> MLP --> ln(IC50)

Both encoders stay differentiable end to end, and both follow the
encode-once-per-forward pattern: the 532 cell lines and 498 drugs are encoded
as whole node sets, then indexed per pair, rather than re-encoded per pair.

**Node ordering is the correctness-critical detail here.** The GNN indexes
`cell_line` and `drug` embeddings by the graph's own node indices, so the
fused omics rows and the drug embedding table must be built in exactly the
graph's `node_ids` order -- not the order the PCA CSVs or the SMILES cache
happen to use. `align_omics_to_graph` and `align_drug_graphs_to_graph` do
that alignment and assert coverage; nothing here assumes the orders agree.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

from src.data.drug_graphs import BatchedMolGraphs, MolGraph, collate_drug_graphs
from src.models.cross_attention_fusion import MultiOmicsCrossAttentionFusion
from src.models.drug_gcn import MolecularGraphEncoder
from src.models.hetero_gnn import HeteroGNN

# Pinned to the fusion settings used throughout the matrix, so this model's
# fusion block is the same one E10/E14 and the molecular-graph runs use.
D_MODEL = 128
NUM_HEADS = 4
FUSION_OUT_DIM = 256
FUSION_DROPOUT = 0.2


def align_omics_to_graph(
    omics: Dict[str, np.ndarray], cell_ids: pd.Index, node_ids: Sequence[str]
) -> Dict[str, torch.Tensor]:
    """Reorder per-modality PCA rows into the graph's `cell_line` node order."""
    missing = [c for c in node_ids if c not in set(cell_ids)]
    if missing:
        raise KeyError(
            f"{len(missing)} graph cell_line nodes have no omics row, e.g. {missing[:5]}"
        )
    row_of = pd.Series(np.arange(len(cell_ids)), index=cell_ids)
    rows = row_of.loc[list(node_ids)].to_numpy()
    return {key: torch.from_numpy(arr[rows]) for key, arr in omics.items()}


def align_drug_graphs_to_graph(
    graphs: Dict[str, MolGraph], node_ids: Sequence[str]
) -> BatchedMolGraphs:
    """Collate molecular graphs into the graph's `drug` node order."""
    missing = [d for d in node_ids if d not in graphs]
    if missing:
        raise KeyError(
            f"{len(missing)} graph drug nodes have no molecular graph, e.g. {missing[:5]}"
        )
    return collate_drug_graphs(graphs, list(node_ids))


class FullArchitectureRegressor(nn.Module):
    """Cross-attention fused omics + GCN-encoded drugs, message-passed over the PPI graph.

    `omics` must already be aligned to the graph's `cell_line` node order and
    `batched` to its `drug` node order (see the two helpers above). Both are
    held as non-persistent buffers: they are constants rebuilt from the CSVs
    and SMILES cache on load, so keeping them out of the state_dict avoids
    both checkpoint bloat and a stale copy overriding a rebuilt one.
    """

    def __init__(
        self,
        omics: Dict[str, torch.Tensor],
        batched: BatchedMolGraphs,
        metadata: Tuple[List[str], List[Tuple[str, str, str]]],
        num_proteins: int,
        variant: str = "gcn",
        hidden_dim: int = 128,
        num_layers: int = 2,
        heads: int = 4,
        dropout: float = 0.2,
        head_hidden_dims: Tuple[int, ...] = (256, 128),
        head_dropout: float = 0.3,
    ):
        super().__init__()
        self.modalities = list(omics)
        for modality, tensor in omics.items():
            self.register_buffer(f"omics_{modality}", tensor, persistent=False)

        self.fusion = MultiOmicsCrossAttentionFusion(
            d_model=D_MODEL,
            num_heads=NUM_HEADS,
            out_dim=FUSION_OUT_DIM,
            dropout=FUSION_DROPOUT,
            modalities=self.modalities,
        )
        self.drug_encoder = MolecularGraphEncoder(batched)
        self.gnn = HeteroGNN(
            metadata=metadata,
            cell_line_dim=FUSION_OUT_DIM,
            drug_dim=self.drug_encoder.out_dim,
            num_proteins=num_proteins,
            variant=variant,
            hidden_dim=hidden_dim,
            num_layers=num_layers,
            heads=heads,
            dropout=dropout,
            head_hidden_dims=head_hidden_dims,
            head_dropout=head_dropout,
        )

    def node_features(self) -> Dict[str, torch.Tensor]:
        """Encode every cell line and drug once: the GNN's `x_dict` input.

        `protein` is absent deliberately -- `HeteroGNN.encode` learns protein
        features through its own `nn.Embedding`, since 03_graph_construction.py
        zero-fills them as a placeholder.
        """
        fused = self.fusion(
            {m: getattr(self, f"omics_{m}") for m in self.modalities}
        )
        return {"cell_line": fused, "drug": self.drug_encoder()}

    def forward(
        self,
        edge_index_dict: Dict[Tuple[str, str, str], torch.Tensor],
        cell_index: torch.Tensor,
        drug_index: torch.Tensor,
    ) -> torch.Tensor:
        """Predict ln_ic50 for a batch of (cell_line, drug) node-index pairs."""
        return self.gnn(self.node_features(), edge_index_dict, cell_index, drug_index)


if __name__ == "__main__":
    from pathlib import Path

    from torch_geometric.data import HeteroData

    from src.data.drug_graphs import build_drug_graphs
    from src.data.experiment_utils import (
        GE_KEY,
        MUT_CNV_KEY,
        PROTEOMICS_KEY,
        load_omics_subset,
    )

    GRAPH_PATH = Path("src/graph/hetero_graph.pt")
    MODALITIES = [GE_KEY, MUT_CNV_KEY, PROTEOMICS_KEY]

    data: HeteroData = torch.load(GRAPH_PATH, weights_only=False)
    cell_node_ids = list(data["cell_line"].node_ids)
    drug_node_ids = list(data["drug"].node_ids)
    print(
        f"[graph] cell_line={len(cell_node_ids)} drug={len(drug_node_ids)} "
        f"protein={len(data['protein'].node_ids)}"
    )

    omics_raw, cell_ids = load_omics_subset(MODALITIES)
    omics = align_omics_to_graph(omics_raw, cell_ids, cell_node_ids)
    batched = align_drug_graphs_to_graph(build_drug_graphs(), drug_node_ids)

    torch.manual_seed(42)
    model = FullArchitectureRegressor(
        omics=omics,
        batched=batched,
        metadata=(list(data.node_types), list(data.edge_types)),
        num_proteins=len(data["protein"].node_ids),
        variant="gcn",
    )
    model.eval()

    edge_dict = {et: data[et].edge_index for et in data.edge_types}
    cell_index = torch.arange(8)
    drug_index = torch.arange(8)

    with torch.no_grad():
        features = model.node_features()
        out = model(edge_dict, cell_index, drug_index)

    print(f"\ncell_line features: {tuple(features['cell_line'].shape)}")
    print(f"drug features:      {tuple(features['drug'].shape)}")
    print(f"prediction:         {tuple(out.shape)}")

    assert features["cell_line"].shape == (len(cell_node_ids), FUSION_OUT_DIM)
    assert features["drug"].shape == (len(drug_node_ids), model.drug_encoder.out_dim)
    assert out.shape == (8,)
    assert torch.isfinite(out).all(), "non-finite predictions"

    # Alignment: the fused row for a cell line must be built from *that* cell
    # line's omics. Perturbing one graph node's input may only move its own row.
    probe = model.modalities[0]
    buffer = getattr(model, f"omics_{probe}")
    original = buffer[3].clone()
    buffer[3] = buffer[3] + 10.0
    with torch.no_grad():
        perturbed = model.node_features()["cell_line"]
    buffer[3] = original
    changed = (~torch.isclose(perturbed, features["cell_line"], atol=1e-6)).any(dim=1)
    assert changed[3] and changed.sum() == 1, (
        f"perturbing cell-line node 3 changed {int(changed.sum())} fused rows, expected 1"
    )
    print(f"omics row -> cell_line node alignment holds ({probe} probe)")

    total_params = sum(p.numel() for p in model.parameters())
    print(f"\nFullArchitectureRegressor parameter count: {total_params:,}")
    print("All full-architecture checks passed.")
