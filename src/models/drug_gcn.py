"""GCN encoder over atom-level molecular graphs, for the drug branch.

The third drug representation in the matrix, alongside the one-hot identity
and the 2048-bit Morgan fingerprint that `FingerprintEncoder`
(`cross_attention_fusion.py`) compresses. Where a fingerprint enumerates which
substructures are present, this encodes how the atoms actually connect --
the representation the MoGraphDRP benchmark uses (their section 2.2.2), and the
one their Table 6 attributes the largest single component effect to.

Hyperparameters follow their Table 1 where stated: 3 GCN layers, ReLU,
dropout 0.4, 128-dim drug embedding, mean pooling over atoms.

**Encode-once design.** `forward()` takes no arguments and returns the whole
(n_drugs, out_dim) embedding table. The ~498 molecular graphs are collated
into one disjoint graph at construction time (see
`src/data/drug_graphs.py`) and held as non-persistent buffers, so a training
step encodes 498 small molecules once and the regressor indexes the result by
drug code -- rather than re-encoding the same molecule once per pair across
~111,799 pairs. The graph is block-diagonal, so this is mathematically
identical to encoding each molecule separately (asserted in the smoke test
below), just far cheaper.
"""

from __future__ import annotations

import torch
import torch.nn as nn
from torch_geometric.nn import GCNConv, global_mean_pool

from src.data.drug_graphs import ATOM_FEATURE_DIM, BatchedMolGraphs

# MoGraphDRP Table 1: 3 GCN layers, dropout 0.4, 128-dim drug embedding.
GCN_HIDDEN_DIM = 128
GCN_OUT_DIM = 128
GCN_NUM_LAYERS = 3
GCN_DROPOUT = 0.4


class MolecularGraphEncoder(nn.Module):
    """Batched molecular graphs -> one pooled embedding per drug.

    `GCNConv` adds self-loops, so a bond-free atom (a bare counter-ion
    fragment, e.g. the sodium in a salt) still contributes its own features
    to the pooled vector instead of vanishing.
    """

    def __init__(
        self,
        batched: BatchedMolGraphs,
        in_dim: int = ATOM_FEATURE_DIM,
        hidden_dim: int = GCN_HIDDEN_DIM,
        out_dim: int = GCN_OUT_DIM,
        num_layers: int = GCN_NUM_LAYERS,
        dropout: float = GCN_DROPOUT,
    ):
        super().__init__()
        if num_layers < 1:
            raise ValueError(f"MolecularGraphEncoder needs >=1 layer, got {num_layers}")

        # Non-persistent: the graph is rebuilt from SMILES on load, so keeping
        # it out of the state_dict avoids bloating every checkpoint with a
        # constant, and avoids a stale copy silently overriding a rebuilt one.
        self.register_buffer("x", batched.x, persistent=False)
        self.register_buffer("edge_index", batched.edge_index, persistent=False)
        self.register_buffer("batch", batched.batch, persistent=False)
        self.n_drugs = batched.n_graphs
        self.out_dim = out_dim

        dims = [in_dim] + [hidden_dim] * (num_layers - 1) + [out_dim]
        self.convs = nn.ModuleList(
            GCNConv(dims[i], dims[i + 1]) for i in range(num_layers)
        )
        self.activation = nn.ReLU()
        self.dropout = nn.Dropout(dropout)

    def forward(self) -> torch.Tensor:
        """-> (n_drugs, out_dim), row `i` being the drug with code `i`."""
        h = self.x
        for conv in self.convs:
            h = self.dropout(self.activation(conv(h, self.edge_index)))
        return global_mean_pool(h, self.batch, size=self.n_drugs)


if __name__ == "__main__":
    from rdkit import Chem

    from src.data.drug_graphs import collate_drug_graphs, mol_to_graph

    SAMPLES = {
        "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "imatinib_like": "CN1CCN(CC1)Cc1ccc(cc1)C(=O)Nc1ccc(C)c(Nc2nccc(n2)-c2cccnc2)c1",
        "benzene": "c1ccccc1",
        "lone_atom": "[Na+]",
    }
    graphs = {
        name: mol_to_graph(Chem.MolFromSmiles(smiles)) for name, smiles in SAMPLES.items()
    }
    order = list(SAMPLES)

    torch.manual_seed(42)
    batched = collate_drug_graphs(graphs, order)
    encoder = MolecularGraphEncoder(batched)
    encoder.eval()  # dropout off, so the isolation check below is exact

    with torch.no_grad():
        table = encoder()

    print(f"\nDrug embedding table: {tuple(table.shape)}")
    assert table.shape == (len(order), encoder.out_dim)
    assert torch.isfinite(table).all(), "non-finite drug embeddings"

    # The pooled embedding of each drug must not depend on which other drugs
    # were in the batch. If this fails, the collated graph is not truly
    # disjoint and molecules are leaking structure into each other.
    for code, name in enumerate(order):
        solo = MolecularGraphEncoder(collate_drug_graphs(graphs, [name]))
        solo.load_state_dict(encoder.state_dict())
        solo.eval()
        with torch.no_grad():
            solo_embedding = solo()
        assert torch.allclose(table[code], solo_embedding[0], atol=1e-6), (
            f"{name}: batched embedding differs from solo encoding"
        )

    total_params = sum(p.numel() for p in encoder.parameters())
    print(f"MolecularGraphEncoder parameter count: {total_params:,}")
    print("\nPer-drug encoding is independent of batch composition. Checks passed.")
