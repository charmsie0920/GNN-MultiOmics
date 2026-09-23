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

Their architecture (§2.1-2.3):

    per-omics 2-layer MLP branches -> concat -> cell vector
    Morgan/PubChem/ESPF encoders + 3-layer molecular-graph GCN -> concat -> drug vector
    cell x drug -> multi-head bilinear attention -> fusion vector -> MLP -> ln(IC50)

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

Switches are provided for each component so the ablation is a flag rather than
a fork: `fusion`, `drug_mode` and `head`.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import torch
import torch.nn as nn

from src.data.drug_graphs import BatchedMolGraphs
from src.models.cross_attention_fusion import (
    FingerprintEncoder,
    MultiOmicsCrossAttentionFusion,
)
from src.models.drug_gcn import MolecularGraphEncoder

# MoGraphDRP Table 1. These differ from this project's inherited defaults
# (lr 1e-3, batch 256, dropout 0.2/0.3) on every line, which is itself one of
# the unattributed differences this module exists to remove.
LATENT_DIM = 128          # their "cell embedding" per branch, and drug embedding size
BILINEAR_HEADS = 4        # "Number of Attention Heads in BAN"
DROPOUT = 0.4
LR = 1e-4
BATCH_SIZE = 128
MAX_EPOCHS = 200
HEAD_HIDDEN_DIMS = (256, 128)


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
    """Their §2.3: feature-level interaction between the cell and drug vectors.

    A full bilinear form needs a `cell_dim x drug_dim` matrix per head, which is
    both large and prone to overfitting, so this uses the standard **low-rank
    factorization** used by bilinear attention networks: project both vectors
    into a shared per-head space and take their Hadamard product, which is
    equivalent to a rank-constrained bilinear map. The paper describes the same
    three stages -- map into a shared space, compute multiple bilinear
    interactions, compress into one representation.

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
        out_dim: int = LATENT_DIM * 2,
        dropout: float = DROPOUT,
    ):
        super().__init__()
        self.heads = heads
        self.hidden_dim = hidden_dim
        self.cell_proj = nn.Linear(cell_dim, hidden_dim * heads)
        self.drug_proj = nn.Linear(drug_dim, hidden_dim * heads)
        self.dropout = nn.Dropout(dropout)
        self.out = nn.Sequential(
            nn.Linear(hidden_dim * heads, out_dim),
            nn.BatchNorm1d(out_dim),
            nn.ReLU(),
        )
        self.out_dim = out_dim

    def forward(self, cell: torch.Tensor, drug: torch.Tensor) -> torch.Tensor:
        batch = cell.shape[0]
        c = self.cell_proj(cell).view(batch, self.heads, self.hidden_dim)
        d = self.drug_proj(drug).view(batch, self.heads, self.hidden_dim)
        interaction = self.dropout(c * d)  # low-rank bilinear, per head
        return self.out(interaction.reshape(batch, -1))


class DualDrugEncoder(nn.Module):
    """Their §2.2.3: fingerprint and molecular graph used *together*, not as alternatives.

    Every experiment in this project so far picked one or the other; the
    benchmark concatenates both, on the argument that fingerprints capture
    local substructure while the graph captures global topology. `drug_mode`
    selects `fingerprint`, `graph`, or `both` so that claim is testable here.
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
            FingerprintEncoder(fp_size=fp_size, hidden_dim=latent_dim,
                               out_dim=latent_dim, dropout=dropout)
            if drug_mode in {"fingerprint", "both"} else None
        )
        self.graph_encoder = (
            MolecularGraphEncoder(batched, out_dim=latent_dim, dropout=dropout)
            if drug_mode in {"graph", "both"} else None
        )
        self.out_dim = latent_dim * (2 if drug_mode == "both" else 1)

    def forward(
        self,
        fingerprints: Optional[torch.Tensor] = None,
        drug_codes: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        parts: List[torch.Tensor] = []
        if self.fingerprint_encoder is not None:
            parts.append(self.fingerprint_encoder(fingerprints))
        if self.graph_encoder is not None:
            parts.append(self.graph_encoder()[drug_codes])
        return torch.cat(parts, dim=-1) if len(parts) > 1 else parts[0]


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
                cell_dim, self.drug_encoder.out_dim, latent_dim, BILINEAR_HEADS,
                latent_dim * 2, dropout,
            )
            prev = self.interaction.out_dim
        else:
            self.interaction = None
            prev = cell_dim + self.drug_encoder.out_dim

        layers: List[nn.Module] = []
        for h in head_hidden_dims:
            layers += [nn.Linear(prev, h), nn.BatchNorm1d(h), nn.ReLU(), nn.Dropout(dropout)]
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
    bil = MultiHeadBilinearAttention(256, 256)
    bil.eval()
    c1, c2 = torch.randn(4, 256), torch.randn(4, 256)
    d1, d2 = torch.randn(4, 256), torch.randn(4, 256)

    additive = nn.Linear(512, 256)

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
    m(omics, fps, codes).sum().backward()
    for part in ("fusion", "drug_encoder.fingerprint_encoder",
                 "drug_encoder.graph_encoder", "interaction", "predictor"):
        norm = sum(float(p.grad.norm()) for n, p in m.named_parameters()
                   if n.startswith(part) and p.grad is not None)
        assert norm > 0, f"no gradient reached {part}"
        print(f"  grad norm {part:<38} {norm:.4f}")

    print("\nAll MoGraphDRP-aligned model checks passed.")
