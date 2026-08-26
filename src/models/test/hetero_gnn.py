import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import HeteroConv, SAGEConv


class HeteroIC50GNN(nn.Module):

    def __init__(self, hidden_dim=128):
        super().__init__()

        # 1. Project mismatched raw feature dimensions into unified hidden space D
        self.cell_proj = nn.Linear(384, hidden_dim)
        self.drug_proj = nn.Linear(2048, hidden_dim)
        self.prot_proj = nn.Linear(128, hidden_dim)

        # 2. Layer 1 Message Passing
        self.conv1 = HeteroConv(
            {
                ("protein", "interacts_with", "protein"): SAGEConv(
                    hidden_dim, hidden_dim
                ),
                ("drug", "targets", "protein"): SAGEConv(
                    hidden_dim, hidden_dim
                ),
                ("protein", "rev_targets", "drug"): SAGEConv(hidden_dim, hidden_dim),
            },
            aggr="sum",
        )

        # 3. Layer 2 Message Passing
        self.conv2 = HeteroConv(
            {
                ("protein", "interacts_with", "protein"): SAGEConv(
                    hidden_dim, hidden_dim
                ),
                ("drug", "targets", "protein"): SAGEConv(
                    hidden_dim, hidden_dim
                ),
                ("protein", "rev_targets", "drug"): SAGEConv(hidden_dim, hidden_dim),
            },
            aggr="sum",
        )

        # 4. IC50 Regression Head
        self.predictor = nn.Sequential(
            nn.Linear(hidden_dim * 2, 128),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
        )

    def forward(self, x_dict, edge_index_dict, cell_idx, drug_idx):
        # Initial projection to hidden_dim
        h_dict = {
            "cell_line": F.relu(self.cell_proj(x_dict["cell_line"])),
            "drug": F.relu(self.drug_proj(x_dict["drug"])),
            "protein": F.relu(self.prot_proj(x_dict["protein"])),
        }

        out = self.conv1(h_dict, edge_index_dict)
        h_dict = {**h_dict, **{k: F.relu(v) for k, v in out.items()}}

        out = self.conv2(h_dict, edge_index_dict)
        h_dict = {**h_dict, **{k: F.relu(v) for k, v in out.items()}}

        # Extract target pairs for this training batch
        batch_cell_emb = h_dict["cell_line"][cell_idx]
        batch_drug_emb = h_dict["drug"][drug_idx]

        # Concatenate and pass through MLP predictor
        pair_representation = torch.cat(
            [batch_cell_emb, batch_drug_emb], dim=-1
        )
        return self.predictor(pair_representation).squeeze(-1)