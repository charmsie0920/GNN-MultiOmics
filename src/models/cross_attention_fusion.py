"""Directional pairwise cross-attention fusion of GE / Mut_CNV / Proteomics embeddings.

Each omics modality is assumed to already be PCA-compressed to a shared
dimension `d_model` (see `src/data/omics_preprocessing.py`). This module
computes directional multi-head cross-attention across all 6 ordered
modality pairs and fuses them into a single cell-line embedding suitable for
binding into PyTorch Geometric node features.

Note on attention weights: with `n_tokens=1` (the default, and what every
result in docs/results.md was produced with) each modality is a single pooled
vector per sample rather than a sequence, so each cross-attention call has
exactly one query and one key and softmax reduces to a weight of 1.0 by
construction. The learning capacity of that configuration comes from the
Q/K/V projections and the residual FFN, not from a non-trivial attention
distribution -- so "cross-attention fusion" is, at `n_tokens=1`, a claim the
architecture does not actually support.

`n_tokens > 1` fixes that by splitting each modality's `d_model` vector into
`n_tokens` contiguous chunks of `d_model // n_tokens` before attending, giving
softmax something to distribute over. The chunks are meaningful rather than
arbitrary: the PCA components feeding this module are ordered by explained
variance, so the first chunk holds the dominant directions and later chunks
progressively finer structure -- attention can then learn to weight coarse
against fine signal, per modality pair.

`n_tokens=1` is kept as the default so existing callers and every recorded
result stay bit-identical; the token count is an ablation axis, not a
migration.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import torch
import torch.nn as nn

from src.data.omics_preprocessing import MODALITIES


class PairwiseCrossAttention(nn.Module):
    """Single directional multi-head cross-attention block with residual + LayerNorm.

    Each modality vector is viewed as `n_tokens` chunks of `d_model //
    n_tokens` and attended token-to-token. At `n_tokens=1` this is the
    original single-token block, where softmax is 1.0 by construction; above
    1, the attention distribution is genuinely learned.
    """

    def __init__(
        self,
        d_model: int = 128,
        num_heads: int = 4,
        dropout: float = 0.1,
        n_tokens: int = 1,
    ):
        super().__init__()
        d_token, remainder = divmod(d_model, n_tokens)
        if n_tokens < 1 or remainder:
            raise ValueError(
                f"n_tokens must be >=1 and divide d_model ({d_model}), got {n_tokens}"
            )
        if d_token % num_heads:
            raise ValueError(
                f"num_heads ({num_heads}) must divide the per-token width "
                f"d_model // n_tokens = {d_token}"
            )
        self.n_tokens = n_tokens
        self.d_token = d_token
        self.mha = nn.MultiheadAttention(embed_dim=d_token, num_heads=num_heads, dropout=dropout, batch_first=True)
        self.norm = nn.LayerNorm(d_token)
        # Off during training: `need_weights=True` costs an extra materialized
        # (B, n_tokens, n_tokens) tensor per block per step.
        self.record_attention = False
        self.last_attention: Optional[torch.Tensor] = None

    def forward(self, query_modality: torch.Tensor, key_value_modality: torch.Tensor) -> torch.Tensor:
        """query_modality, key_value_modality: (B, d_model) -> (B, d_model)."""
        batch = query_modality.shape[0]
        q = query_modality.view(batch, self.n_tokens, self.d_token)
        kv = key_value_modality.view(batch, self.n_tokens, self.d_token)
        attn_out, weights = self.mha(q, kv, kv, need_weights=self.record_attention)
        if self.record_attention:
            self.last_attention = weights.detach()
        return self.norm(q + attn_out).reshape(batch, -1)


class MultiOmicsCrossAttentionFusion(nn.Module):
    """Fuses an arbitrary subset (>=2) of omics modality embeddings into one cell-line vector.

    Defaults to all 3 modalities (GE, Mut_CNV, Proteomics) for backward
    compatibility with existing callers. Pass `modalities` explicitly to
    fuse a 2-of-3 subset (needed by the omics-ablation experiment matrix,
    docs/plan/experiment_matrix_plan.md) -- `self.pairs` and `output_proj`'s
    input width are both derived from the given list, so a 2-modality call
    builds 2 pairwise-attention blocks instead of 6.
    """

    def __init__(
        self,
        d_model: int = 128,
        num_heads: int = 4,
        out_dim: int = 256,
        dropout: float = 0.2,
        modalities: Optional[Sequence[str]] = None,
        n_tokens: int = 1,
    ):
        super().__init__()
        self.n_tokens = n_tokens
        self.modalities = list(modalities) if modalities is not None else list(MODALITIES)
        if len(self.modalities) < 2:
            raise ValueError(
                f"MultiOmicsCrossAttentionFusion needs >=2 modalities to cross-attend "
                f"across, got {self.modalities!r}"
            )
        self.pairs: List[Tuple[str, str]] = [
            (i, j) for i in self.modalities for j in self.modalities if i != j
        ]

        self.cross_attn = nn.ModuleDict(
            {
                f"{i}->{j}": PairwiseCrossAttention(d_model, num_heads, dropout, n_tokens)
                for i, j in self.pairs
            }
        )
        self.ffn = nn.Sequential(
            nn.Linear(d_model, d_model * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model * 4, d_model),
        )
        self.ffn_norm = nn.LayerNorm(d_model)

        concat_dim = d_model * len(self.pairs)
        self.output_proj = nn.Sequential(
            nn.Linear(concat_dim, out_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.LayerNorm(out_dim),
        )

    def set_record_attention(self, enabled: bool = True) -> None:
        """Toggle attention-weight capture on every pairwise block.

        Left off during training; turn it on for an evaluation pass to inspect
        what the fusion actually learned. At `n_tokens=1` the captured weights
        are all exactly 1.0 -- which is the point being demonstrated.
        """
        for block in self.cross_attn.values():
            block.record_attention = enabled
            if not enabled:
                block.last_attention = None

    def attention_weights(self) -> Dict[str, torch.Tensor]:
        """Per-pair attention from the last forward, as {'GE->Proteomics': (B, T, T)}."""
        return {
            name: block.last_attention
            for name, block in self.cross_attn.items()
            if block.last_attention is not None
        }

    def forward(self, omics: Dict[str, torch.Tensor]) -> torch.Tensor:
        """omics: {modality: (B, d_model)} for each modality in `self.modalities` -> (B, out_dim)."""
        pair_outputs = []
        for i, j in self.pairs:
            attn_out = self.cross_attn[f"{i}->{j}"](omics[i], omics[j])
            ffn_out = self.ffn_norm(attn_out + self.ffn(attn_out))
            pair_outputs.append(ffn_out)
        fused = torch.cat(pair_outputs, dim=-1)
        return self.output_proj(fused)


class FingerprintEncoder(nn.Module):
    """Compresses a sparse 2048-bit Morgan fingerprint to a dense low-dim vector.

    Concatenating raw 2048 sparse bits directly onto a 256-dim continuous
    fused-omics embedding would let the drug block dominate the head's input
    purely by width, so fingerprint-mode experiments project it down first.
    Not needed for tree models (RF handles the raw 2048-dim input fine).
    """

    def __init__(self, fp_size: int = 2048, hidden_dim: int = 128, out_dim: int = 128, dropout: float = 0.2):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(fp_size, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, out_dim),
        )
        self.out_dim = out_dim

    def forward(self, fingerprint: torch.Tensor) -> torch.Tensor:
        """fingerprint: (B, fp_size) -> (B, out_dim)."""
        return self.net(fingerprint)


if __name__ == "__main__":
    import numpy as np
    import pandas as pd

    from src.data.omics_preprocessing import GE_KEY, MUT_CNV_KEY, PROTEOMICS_KEY, OmicsPreprocessingPipeline

    N_SAMPLES = 200  # must exceed d_target for PCA(d_target) to be achievable
    D_TARGET = 128
    OUT_DIM = 256
    rng = np.random.default_rng(seed=42)

    cell_line_ids = [f"CELL_{i:04d}" for i in range(N_SAMPLES)]
    ge_df = pd.DataFrame(rng.lognormal(mean=1.0, sigma=1.0, size=(N_SAMPLES, 10_000)), index=cell_line_ids)
    mut_cnv_df = pd.DataFrame(rng.normal(size=(N_SAMPLES, 5_000)), index=cell_line_ids)
    proteomics_df = pd.DataFrame(rng.normal(size=(N_SAMPLES, 2_000)), index=cell_line_ids)
    proteomics_df.iloc[::10, ::7] = np.nan  # inject missing values to exercise median imputation

    print("Input shapes:")
    print(f"  {GE_KEY}: {ge_df.shape}")
    print(f"  {MUT_CNV_KEY}: {mut_cnv_df.shape}")
    print(f"  {PROTEOMICS_KEY}: {proteomics_df.shape}")

    pipeline = OmicsPreprocessingPipeline(d_target=D_TARGET, model_dir="models/omics_pca")
    reduced = pipeline.fit_transform({GE_KEY: ge_df, MUT_CNV_KEY: mut_cnv_df, PROTEOMICS_KEY: proteomics_df})
    pipeline.save()

    print("\nPost-PCA (intermediate pair) shapes:")
    for modality, arr in reduced.items():
        print(f"  {modality}: {arr.shape}")

    omics_tensors = {modality: torch.from_numpy(arr) for modality, arr in reduced.items()}

    fusion_model = MultiOmicsCrossAttentionFusion(d_model=D_TARGET, num_heads=4, out_dim=OUT_DIM, dropout=0.2)
    fusion_model.eval()
    with torch.no_grad():
        h_cell = fusion_model(omics_tensors)

    print(f"\nFinal fused cell-line embedding shape: {tuple(h_cell.shape)}")

    total_params = sum(p.numel() for p in fusion_model.parameters())
    trainable_params = sum(p.numel() for p in fusion_model.parameters() if p.requires_grad)
    approx_vram_mb = total_params * 4 / (1024 ** 2)  # fp32 parameters only, excludes activations/optimizer state
    print(f"\nMultiOmicsCrossAttentionFusion parameter count: {total_params:,} (trainable: {trainable_params:,})")
    print(f"Approx. parameter memory (fp32): {approx_vram_mb:.2f} MB")
