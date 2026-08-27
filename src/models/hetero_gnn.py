"""Heterogeneous GNN (GCN / GAT variants) over the cell_line-drug-protein graph.

The first graph-based model in the project: message passing over
`data/processed/hetero_graph.pt` (built by 03_graph_construction.py, linked
by 04_link_cell_lines.py), rather than the flat per-cell-line feature vectors
every prior baseline used.

Two conv variants are exposed so the ablation can test which suits noisy
biological graphs -- the MoGraphDRP paper found GCN beat GAT on their drug
molecular graphs (GAT overfit), which is worth checking on our very different
graph topology (a PPI network, not a molecule).

Protein nodes carry no real features (03_graph_construction.py zero-fills
them as an explicit placeholder), so this module gives them a learnable
`nn.Embedding` instead -- the step that doc flagged as the GNN stage's job.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import torch
import torch.nn as nn
from torch_geometric.nn import GATConv, HeteroConv, SAGEConv

NodeType = str
EdgeType = Tuple[str, str, str]


def _make_conv(variant: str, in_dim: int, out_dim: int, heads: int, dropout: float) -> nn.Module:
    """One conv layer for a single edge type.

    SAGEConv is used for the "GCN" variant rather than GCNConv because
    GCNConv requires a symmetric, same-node-type graph (it normalizes by node
    degree over a single node set); our edge types are bipartite
    (cell_line->protein, drug->protein), which SAGEConv handles natively via
    its (src, dst) tuple input dims.
    """
    if variant == "gcn":
        return SAGEConv((-1, -1), out_dim)
    if variant == "gat":
        return GATConv((-1, -1), out_dim // heads, heads=heads, dropout=dropout, add_self_loops=False)
    raise ValueError(f"Unknown conv variant: {variant!r}")


class HeteroGNN(nn.Module):
    """Message passing over the hetero graph -> (cell_line, drug) pair -> ln_ic50.

    Each edge type gets its own conv, and every relation is also added in
    reverse so information flows in both directions (e.g. a drug's target
    protein can inform the drug embedding, and a cell line's mutated proteins
    can inform the cell-line embedding).
    """

    def __init__(
        self,
        metadata: Tuple[List[NodeType], List[EdgeType]],
        cell_line_dim: int,
        drug_dim: int,
        num_proteins: int,
        variant: str = "gat",
        hidden_dim: int = 128,
        num_layers: int = 2,
        heads: int = 4,
        dropout: float = 0.2,
        head_hidden_dims: Tuple[int, ...] = (256, 128),
        head_dropout: float = 0.3,
    ):
        super().__init__()
        node_types, edge_types = metadata

        # Project the heterogeneous input feature dims onto one shared width.
        self.cell_line_proj = nn.Linear(cell_line_dim, hidden_dim)
        self.drug_proj = nn.Linear(drug_dim, hidden_dim)
        # Proteins have placeholder (zero) features -- learn them instead.
        self.protein_embedding = nn.Embedding(num_proteins, hidden_dim)

        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            conv_map: Dict[EdgeType, nn.Module] = {}
            for edge_type in edge_types:
                src, rel, dst = edge_type
                conv_map[edge_type] = _make_conv(variant, hidden_dim, hidden_dim, heads, dropout)
                if src != dst:
                    conv_map[(dst, f"rev_{rel}", src)] = _make_conv(
                        variant, hidden_dim, hidden_dim, heads, dropout
                    )
            self.convs.append(HeteroConv(conv_map, aggr="sum"))

        self.norm = nn.ModuleDict({nt: nn.LayerNorm(hidden_dim) for nt in node_types})
        self.dropout = nn.Dropout(dropout)

        layers: list[nn.Module] = []
        prev = hidden_dim * 2  # cell_line embedding ++ drug embedding
        for h in head_hidden_dims:
            layers += [nn.Linear(prev, h), nn.BatchNorm1d(h), nn.ReLU(), nn.Dropout(head_dropout)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.head = nn.Sequential(*layers)

    def encode(
        self, x_dict: Dict[NodeType, torch.Tensor], edge_index_dict: Dict[EdgeType, torch.Tensor]
    ) -> Dict[NodeType, torch.Tensor]:
        """Run message passing, returning one embedding per node of each type."""
        h_dict = {
            "cell_line": self.cell_line_proj(x_dict["cell_line"]),
            "drug": self.drug_proj(x_dict["drug"]),
            "protein": self.protein_embedding.weight,
        }

        # Mirror every relation so messages flow both ways.
        full_edges = dict(edge_index_dict)
        for (src, rel, dst), edge_index in edge_index_dict.items():
            if src != dst:
                full_edges[(dst, f"rev_{rel}", src)] = edge_index.flip(0)

        for conv in self.convs:
            out = conv(h_dict, full_edges)
            # HeteroConv only returns node types that received messages; keep
            # the previous embedding for any type that didn't (residual-style).
            h_dict = {
                nt: self.dropout(torch.relu(self.norm[nt](out[nt]))) if nt in out else h_dict[nt]
                for nt in h_dict
            }
        return h_dict

    def forward(
        self,
        x_dict: Dict[NodeType, torch.Tensor],
        edge_index_dict: Dict[EdgeType, torch.Tensor],
        cell_index: torch.Tensor,
        drug_index: torch.Tensor,
    ) -> torch.Tensor:
        """Predict ln_ic50 for a batch of (cell_line, drug) index pairs."""
        h_dict = self.encode(x_dict, edge_index_dict)
        pair = torch.cat([h_dict["cell_line"][cell_index], h_dict["drug"][drug_index]], dim=-1)
        return self.head(pair).squeeze(-1)
