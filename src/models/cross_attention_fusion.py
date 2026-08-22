"""Directional pairwise cross-attention fusion of GE / Mut_CNV / Proteomics embeddings.

Each omics modality is assumed to already be PCA-compressed to a shared
dimension `d_model` (see `src/data/omics_preprocessing.py`). This module
computes directional multi-head cross-attention across all 6 ordered
modality pairs and fuses them into a single cell-line embedding suitable for
binding into PyTorch Geometric node features.

Note on attention weights: because each modality is a single pooled vector
per sample (not a sequence of tokens), each cross-attention call has exactly
one query and one key, so softmax reduces to a weight of 1.0 by construction.
The learning capacity of this block therefore comes from the Q/K/V
projections and the residual FFN, not from a non-trivial attention
distribution.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import torch
import torch.nn as nn

from src.data.omics_preprocessing import MODALITIES


class PairwiseCrossAttention(nn.Module):
    """Single directional multi-head cross-attention block with residual + LayerNorm."""

    def __init__(self, d_model: int = 128, num_heads: int = 4, dropout: float = 0.1):
        super().__init__()
        self.mha = nn.MultiheadAttention(embed_dim=d_model, num_heads=num_heads, dropout=dropout, batch_first=True)
        self.norm = nn.LayerNorm(d_model)

    def forward(self, query_modality: torch.Tensor, key_value_modality: torch.Tensor) -> torch.Tensor:
        """query_modality, key_value_modality: (B, d_model) -> (B, d_model)."""
        q = query_modality.unsqueeze(1)
        kv = key_value_modality.unsqueeze(1)
        attn_out, _ = self.mha(q, kv, kv)
        attn_out = attn_out.squeeze(1)
        return self.norm(query_modality + attn_out)


class MultiOmicsCrossAttentionFusion(nn.Module):
    """Fuses GE, Mut_CNV, and Proteomics embeddings into one cell-line vector."""

    def __init__(self, d_model: int = 128, num_heads: int = 4, out_dim: int = 256, dropout: float = 0.2):
        super().__init__()
        self.modalities = list(MODALITIES)
        self.pairs: List[Tuple[str, str]] = [
            (i, j) for i in self.modalities for j in self.modalities if i != j
        ]

        self.cross_attn = nn.ModuleDict(
            {f"{i}->{j}": PairwiseCrossAttention(d_model, num_heads, dropout) for i, j in self.pairs}
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

    def forward(self, omics: Dict[str, torch.Tensor]) -> torch.Tensor:
        """omics: {'GE': (B, d_model), 'Mut_CNV': (B, d_model), 'Proteomics': (B, d_model)} -> (B, out_dim)."""
        pair_outputs = []
        for i, j in self.pairs:
            attn_out = self.cross_attn[f"{i}->{j}"](omics[i], omics[j])
            ffn_out = self.ffn_norm(attn_out + self.ffn(attn_out))
            pair_outputs.append(ffn_out)
        fused = torch.cat(pair_outputs, dim=-1)
        return self.output_proj(fused)


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
