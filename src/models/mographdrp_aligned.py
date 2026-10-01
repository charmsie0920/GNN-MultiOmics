"""A MoGraphDRP-aligned baseline, so comparisons against the benchmark are fair.

Every model in this project so far differs from MoGraphDRP on several axes at
once -- fusion, drug representation, prediction head and hyperparameters --
which makes the RMSE gap uninterpretable: it cannot be attributed to any one
choice. This module reconstructs their architecture as closely as the
available data allows, so that:

1. the benchmark can be reproduced under **their** protocol (random pair split)
   and re-measured under **ours** (grouped by cell line), and
2. this project's own contributions -- cross-attention fusion, the PPI graph --
   become single-component deltas on a fair base rather than one large
   unattributable difference.

Their architecture, as implemented in their released code (not only the paper
text -- the two disagree on the number of bilinear heads and on how the two
drug views are combined, and the code is what produced their numbers):

    per-omics 2-layer MLP branches -> concat -> cell vector
    Morgan/PubChem/ESPF encoders + 3-layer molecular-graph GCN -> gated sum -> drug vector
    cell x drug -> bilinear attention (BAN layer) -> fusion vector -> MLP -> ln(IC50)

**What this does not match, and why the residual gap is not purely
architectural:**

- They use 4 omics types (expression, mutation, methylation, GSVA pathway
  scores); this project has 3 (GE, Mut_CNV, proteomics). Proteomics is an
  addition of ours; methylation and pathway activity would need new data.
- They fuse 3 fingerprint types (Morgan 2048, PubChem 881, ESPF 2586); only
  Morgan is available here. PubChem and ESPF are not derivable with RDKit
  alone.
- They filter to ~600-700 COSMIC cancer genes per branch; this project applies
  PCA to 128 components (retaining 68-75% of variance).
- Their atoms carry a 78-dim sum-normalised feature vector; this project keeps
  its own atom vocabulary (`src/data/drug_graphs.py`). The GCN layer widths
  follow their in/2x/4x pattern on our input width.
- Their training code sets no random seed, so their published figure is one
  unreproducible draw of the split and the initialisation.

Switches are provided for each component so the ablation is a flag rather than
a fork: `fusion`, `drug_mode` and `head`.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence

import torch
import torch.nn as nn
from torch.nn.utils import weight_norm
from torch_geometric.nn import GCNConv, global_add_pool

from src.data.drug_graphs import ATOM_FEATURE_DIM, BatchedMolGraphs
from src.models.cross_attention_fusion import MultiOmicsCrossAttentionFusion

# MoGraphDRP's `config.py` and `models/drug_gnn.py`. These differ from this
# project's inherited defaults (lr 1e-3, batch 256, dropout 0.2/0.3) on every
# line, which is itself one of the unattributed differences this module exists
# to remove.
LATENT_DIM = 128          # their "cell embedding" per branch, and drug embedding size
BILINEAR_HEADS = 3        # `ban_heads=3` in their code (their Table 1 says 4)
BAN_K = 3                 # BAN pooling factor
BAN_DROPOUT = 0.5
GCN_FC_DIM = 1024
DROPOUT = 0.4
LR = 1e-4
BATCH_SIZE = 128
MAX_EPOCHS = 200
HEAD_HIDDEN_DIMS = (512, 128)


class OmicsBranchEncoder(nn.Module):
    """Their §2.1: one independent 2-layer branch per omics type, then concat.

    Each branch is Linear -> BatchNorm -> ReLU -> Dropout -> Linear -> ReLU,
    projecting into a shared `latent_dim` space. The branches are combined by
    plain concatenation -- the paper reports that learned fusion weights were
    *less* stable than concatenation in their experiments, which is precisely
    the claim this project's cross-attention arm tests.
    """

    def __init__(
        self,
        modalities: Sequence[str],
        in_dim: int = 128,
        latent_dim: int = LATENT_DIM,
        dropout: float = DROPOUT,
    ):
        super().__init__()
        self.modalities = list(modalities)
        self.branches = nn.ModuleDict({
            m: nn.Sequential(
                nn.Linear(in_dim, latent_dim * 2),
                nn.BatchNorm1d(latent_dim * 2),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(latent_dim * 2, latent_dim),
                nn.ReLU(),
            )
            for m in self.modalities
        })
        self.out_dim = latent_dim * len(self.modalities)

    def forward(self, omics: Dict[str, torch.Tensor]) -> torch.Tensor:
        return torch.cat([self.branches[m](omics[m]) for m in self.modalities], dim=-1)


class MultiHeadBilinearAttention(nn.Module):
    """Their §2.3, ported from the `BANLayer` in their `models/fusion.py`.

    Both vectors are projected into a shared `hidden_dim * k` space by a
    Dropout -> weight-normed Linear -> ReLU net. Each head scores the pair with
    a learned bilinear form over that space, the Hadamard product is weighted
    by the summed head scores, sum-pooled in groups of `k` down to
    `hidden_dim`, and batch-normalised.

    Their layer is written for token sequences, but they call it with exactly
    one drug token and one cell token, so each head's attention map is a single
    scalar. This is that single-token case written out directly.

    Unlike concatenation, every output unit depends on a *product* of a cell
    feature and a drug feature, so the head can express "this omics pattern
    matters only for this kind of molecule" -- which a concatenation head
    followed by an MLP can only approximate.
    """

    def __init__(
        self,
        cell_dim: int,
        drug_dim: int,
        hidden_dim: int = LATENT_DIM,
        heads: int = BILINEAR_HEADS,
        dropout: float = BAN_DROPOUT,
        k: int = BAN_K,
    ):
        super().__init__()
        self.hidden_dim, self.k = hidden_dim, k

        def proj(in_dim: int) -> nn.Sequential:
            return nn.Sequential(
                nn.Dropout(dropout),
                weight_norm(nn.Linear(in_dim, hidden_dim * k), dim=None),
                nn.ReLU(),
            )

        self.drug_net = proj(drug_dim)
        self.cell_net = proj(cell_dim)
        self.h_mat = nn.Parameter(torch.randn(heads, hidden_dim * k))
        self.h_bias = nn.Parameter(torch.randn(heads))
        self.bn = nn.BatchNorm1d(hidden_dim)
        self.out_dim = hidden_dim

    def forward(self, cell: torch.Tensor, drug: torch.Tensor) -> torch.Tensor:
        joint = self.drug_net(drug) * self.cell_net(cell)        # (B, hidden*k)
        att = joint @ self.h_mat.t() + self.h_bias               # (B, heads)
        fused = att.sum(dim=1, keepdim=True) * joint
        fused = fused.view(-1, self.hidden_dim, self.k).sum(dim=-1)
        return self.bn(fused)


class MorganEncoder(nn.Module):
    """Their `build_dynamic_encoder`: halve the width until the bottleneck.

    For 2048 -> 128 that is 2048 -> 1024 -> 512 -> 256 -> 128, each step
    Linear -> ReLU -> Dropout. They sum three such encoders (Morgan, ESPF,
    PubChem); only Morgan is available here.
    """

    def __init__(self, fp_size: int = 2048, out_dim: int = LATENT_DIM, dropout: float = DROPOUT):
        super().__init__()
        n_layers = max(1, int(math.log2(fp_size / out_dim)))
        ratio = (out_dim / fp_size) ** (1.0 / n_layers)
        sizes = [fp_size]
        for _ in range(1, n_layers):
            sizes.append(int(sizes[-1] * ratio))
        sizes.append(out_dim)

        layers: List[nn.Module] = []
        for a, b in zip(sizes[:-1], sizes[1:]):
            layers += [nn.Linear(a, b), nn.ReLU(), nn.Dropout(dropout)]
        self.net = nn.Sequential(*layers)
        self.out_dim = out_dim

    def forward(self, fingerprint: torch.Tensor) -> torch.Tensor:
        return self.net(fingerprint)


class MolGraphGCN(nn.Module):
    """Their `GCNNetMultiOmics` graph sub-branch.

    Three GCN layers widening in -> in -> 2*in -> 4*in with ReLU and no
    dropout, **sum** pooling over atoms, then 4*in -> 1024 -> out. Encodes all
    drugs once per step and is indexed by drug code, like
    `MolecularGraphEncoder` (`src/models/drug_gcn.py`), which is kept unchanged
    for the other experiments.
    """

    def __init__(
        self,
        batched: BatchedMolGraphs,
        in_dim: int = ATOM_FEATURE_DIM,
        out_dim: int = LATENT_DIM,
        fc_dim: int = GCN_FC_DIM,
    ):
        super().__init__()
        self.register_buffer("x", batched.x, persistent=False)
        self.register_buffer("edge_index", batched.edge_index, persistent=False)
        self.register_buffer("batch", batched.batch, persistent=False)
        self.n_drugs = batched.n_graphs
        self.out_dim = out_dim

        self.convs = nn.ModuleList([
            GCNConv(in_dim, in_dim),
            GCNConv(in_dim, in_dim * 2),
            GCNConv(in_dim * 2, in_dim * 4),
        ])
        self.fc1 = nn.Linear(in_dim * 4, fc_dim)
        self.fc2 = nn.Linear(fc_dim, out_dim)
        self.activation = nn.ReLU()

    def forward(self) -> torch.Tensor:
        """-> (n_drugs, out_dim), row `i` being the drug with code `i`."""
        h = self.x
        for conv in self.convs:
            h = self.activation(conv(h, self.edge_index))
        h = global_add_pool(h, self.batch, size=self.n_drugs)
        return self.fc2(self.activation(self.fc1(h)))


class DualDrugEncoder(nn.Module):
    """Their §2.2.3: fingerprint and molecular graph used *together*, not as alternatives.

    Every experiment in this project so far picked one or the other; the
    benchmark combines both, on the argument that fingerprints capture local
    substructure while the graph captures global topology. `drug_mode` selects
    `fingerprint`, `graph`, or `both` so that claim is testable here.

    `both` follows their `DrugFusionModule` ("simple" mode): a learned
    per-dimension sigmoid gate on each view, summed, then one Linear -- so the
    drug vector stays `latent_dim` wide rather than doubling. A single view is
    passed through as-is.
    """

    def __init__(
        self,
        drug_mode: str,
        batched: Optional[BatchedMolGraphs] = None,
        fp_size: int = 2048,
        latent_dim: int = LATENT_DIM,
        dropout: float = DROPOUT,
    ):
        super().__init__()
        if drug_mode not in {"fingerprint", "graph", "both"}:
            raise ValueError(f"drug_mode must be fingerprint/graph/both, got {drug_mode!r}")
        if drug_mode in {"graph", "both"} and batched is None:
            raise ValueError(f"drug_mode={drug_mode!r} needs the batched molecular graphs")
        self.drug_mode = drug_mode

        self.fingerprint_encoder = (
            MorganEncoder(fp_size=fp_size, out_dim=latent_dim, dropout=dropout)
            if drug_mode in {"fingerprint", "both"} else None
        )
        self.graph_encoder = (
            MolGraphGCN(batched, out_dim=latent_dim)
            if drug_mode in {"graph", "both"} else None
        )
        if drug_mode == "both":
            self.gate = nn.Parameter(torch.randn(latent_dim * 2))
            self.fuse = nn.Linear(latent_dim, latent_dim)
        self.out_dim = latent_dim

    def forward(
        self,
        fingerprints: Optional[torch.Tensor] = None,
        drug_codes: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        if self.drug_mode == "fingerprint":
            return self.fingerprint_encoder(fingerprints)
        graph = self.graph_encoder()[drug_codes]
        if self.drug_mode == "graph":
            return graph
        w_graph, w_fp = torch.sigmoid(self.gate).chunk(2)
        return self.fuse(graph * w_graph + self.fingerprint_encoder(fingerprints) * w_fp)


class MoGraphDRPAligned(nn.Module):
    """The aligned baseline, with each divergence from the benchmark as a switch.

    `fusion="concat"` + `head="bilinear"` + `drug_mode="both"` is the
    MoGraphDRP-aligned configuration. Changing one switch at a time measures
    one component.
    """

    def __init__(
        self,
        modalities: Sequence[str],
        drug_mode: str = "both",
        fusion: str = "concat",
        head: str = "bilinear",
        batched: Optional[BatchedMolGraphs] = None,
        omics_in_dim: int = 128,
        latent_dim: int = LATENT_DIM,
        dropout: float = DROPOUT,
        head_hidden_dims: Sequence[int] = HEAD_HIDDEN_DIMS,
    ):
        super().__init__()
        if fusion not in {"concat", "cross_attention"}:
            raise ValueError(f"fusion must be concat/cross_attention, got {fusion!r}")
        if head not in {"bilinear", "mlp"}:
            raise ValueError(f"head must be bilinear/mlp, got {head!r}")
        self.modalities, self.fusion_mode, self.head_mode = list(modalities), fusion, head

        if fusion == "concat":
            self.fusion = OmicsBranchEncoder(modalities, omics_in_dim, latent_dim, dropout)
        else:
            self.fusion = MultiOmicsCrossAttentionFusion(
                d_model=omics_in_dim, num_heads=4, out_dim=latent_dim * len(modalities),
                dropout=dropout, modalities=list(modalities),
            )
            self.fusion.out_dim = latent_dim * len(modalities)

        self.drug_encoder = DualDrugEncoder(drug_mode, batched, latent_dim=latent_dim, dropout=dropout)
        cell_dim = self.fusion.out_dim

        if head == "bilinear":
            self.interaction = MultiHeadBilinearAttention(
                cell_dim, self.drug_encoder.out_dim, latent_dim,
            )
            prev = self.interaction.out_dim
        else:
            self.interaction = None
            prev = cell_dim + self.drug_encoder.out_dim

        # Their `models/mlp.py`: Linear -> ReLU -> BatchNorm per layer, no dropout.
        layers: List[nn.Module] = []
        for h in head_hidden_dims:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.BatchNorm1d(h)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.predictor = nn.Sequential(*layers)

    def interaction_vector(
        self,
        omics: Dict[str, torch.Tensor],
        fingerprints: Optional[torch.Tensor] = None,
        drug_codes: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """The drug-cell representation fed to the predictor.

        Exposed separately because the benchmark's XGBoost stage consumes this
        vector alongside the scalar prediction (their §2.4).
        """
        cell = self.fusion(omics)
        drug = self.drug_encoder(fingerprints, drug_codes)
        return self.interaction(cell, drug) if self.interaction is not None \
            else torch.cat([cell, drug], dim=-1)

    def forward(
        self,
        omics: Dict[str, torch.Tensor],
        fingerprints: Optional[torch.Tensor] = None,
        drug_codes: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        return self.predictor(
            self.interaction_vector(omics, fingerprints, drug_codes)
        ).squeeze(-1)


if __name__ == "__main__":
    from rdkit import Chem

    from src.data.drug_graphs import collate_drug_graphs, mol_to_graph

    SMILES = {"a": "CC(=O)Oc1ccccc1C(=O)O", "b": "c1ccccc1", "c": "[Na+]",
              "d": "COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC"}
    graphs = {k: mol_to_graph(Chem.MolFromSmiles(v)) for k, v in SMILES.items()}
    batched = collate_drug_graphs(graphs, list(SMILES))

    B, MODS = 16, ["GE", "Proteomics"]
    omics = {m: torch.randn(B, 128) for m in MODS}
    fps = torch.randint(0, 2, (B, 2048)).float()
    codes = torch.randint(0, batched.n_graphs, (B,))

    print(f"{'configuration':<46}{'out':>8}{'params':>12}")
    print("-" * 66)
    for fusion in ("concat", "cross_attention"):
        for head in ("bilinear", "mlp"):
            for drug_mode in ("both", "fingerprint", "graph"):
                torch.manual_seed(42)
                m = MoGraphDRPAligned(MODS, drug_mode, fusion, head, batched)
                m.eval()
                with torch.no_grad():
                    out = m(omics, fps, codes)
                assert out.shape == (B,) and torch.isfinite(out).all()
                label = f"{fusion} + {head} + {drug_mode}"
                print(f"{label:<46}{str(tuple(out.shape)):>8}"
                      f"{sum(p.numel() for p in m.parameters()):>12,}")

    # The bilinear head must be genuinely *interactive*, not merely a fancier
    # way of adding two vectors. The test is the mixed difference
    #     d = f(c1,d1) - f(c1,d2) - f(c2,d1) + f(c2,d2)
    # which cancels to exactly zero for any head that is additive in its two
    # arguments (concat -> linear), and does not for a bilinear one.
    torch.manual_seed(0)
    bil = MultiHeadBilinearAttention(256, 128)
    bil.eval()
    c1, c2 = torch.randn(4, 256), torch.randn(4, 256)
    d1, d2 = torch.randn(4, 128), torch.randn(4, 128)

    additive = nn.Linear(384, 128)

    def add_f(c, d):
        return additive(torch.cat([c, d], dim=-1))

    with torch.no_grad():
        f11, f12, f21, f22 = bil(c1, d1), bil(c1, d2), bil(c2, d1), bil(c2, d2)
        mixed_bilinear = (f11 - f12 - f21 + f22).abs().mean()
        mixed_additive = (add_f(c1, d1) - add_f(c1, d2)
                          - add_f(c2, d1) + add_f(c2, d2)).abs().mean()

    assert not torch.allclose(f11, f12), "output ignores the drug"
    assert not torch.allclose(f11, f21), "output ignores the cell line"
    assert mixed_additive < 1e-5, "additive control should cancel exactly"
    assert mixed_bilinear > 1e-2, "bilinear head is behaving additively"
    print(f"\nmixed difference — bilinear {mixed_bilinear:.4f} vs additive "
          f"{mixed_additive:.2e}: the head models a genuine interaction")

    # Gradients must reach every branch of the aligned configuration.
    torch.manual_seed(42)
    m = MoGraphDRPAligned(MODS, "both", "concat", "bilinear", batched)
    # Squared, not summed: the predictor now ends BatchNorm -> Linear, and a
    # plain sum over the batch is constant through BatchNorm (zero gradient).
    m(omics, fps, codes).pow(2).sum().backward()
    for part in ("fusion", "drug_encoder.fingerprint_encoder",
                 "drug_encoder.graph_encoder", "interaction", "predictor"):
        norm = sum(float(p.grad.norm()) for n, p in m.named_parameters()
                   if n.startswith(part) and p.grad is not None)
        assert norm > 0, f"no gradient reached {part}"
        print(f"  grad norm {part:<38} {norm:.4f}")

    print("\nAll MoGraphDRP-aligned model checks passed.")
