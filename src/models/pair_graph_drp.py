"""PairGraphDRP: the frozen MoGraphDRP-aligned base plus a pair-specific input.

The base (`src/models/mographdrp_aligned.py`) encodes the cell line and the
drug separately and combines the two vectors. Nothing in it sees how *this*
cell line's mutated proteins relate to *this* drug's target proteins. This
model adds a per-pair vector carrying that relation and concatenates it with
the base's interaction vector before the predictor.

Phase 3 (the gate) supplies that vector as hand-built features
(`src/data/pair_features.py`). Phase 4 replaces it with learned attention over
the PPI graph, in this file.

The base is subclassed, not copied: encoders, fusion, head and predictor
layout all come from `MoGraphDRPAligned`, whose file and defaults are frozen.
Only the first predictor layer is widened to accept the extra columns.
"""

from __future__ import annotations

from typing import Dict, Optional

import torch
import torch.nn as nn

from src.models.mographdrp_aligned import MoGraphDRPAligned


class PairGraphDRP(MoGraphDRPAligned):
    """`MoGraphDRPAligned` with `pair_dim` extra inputs to the predictor.

    Takes every argument of the base, by keyword, plus `pair_dim`.
    """

    def __init__(self, *args, pair_dim: int, **kwargs):
        super().__init__(*args, **kwargs)
        if pair_dim < 1:
            raise ValueError(f"pair_dim must be positive, got {pair_dim}")
        self.pair_dim = pair_dim
        first = self.predictor[0]
        self.predictor[0] = nn.Linear(first.in_features + pair_dim, first.out_features)

    def forward(
        self,
        omics: Dict[str, torch.Tensor],
        fingerprints: Optional[torch.Tensor] = None,
        drug_codes: Optional[torch.Tensor] = None,
        pair: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        if pair is None:
            raise ValueError("PairGraphDRP needs the per-pair vector")
        joint = self.interaction_vector(omics, fingerprints, drug_codes)
        return self.predictor(torch.cat([joint, pair], dim=-1)).squeeze(-1)


if __name__ == "__main__":
    B, MODS, PAIR_DIM = 16, ["GE", "Mut_CNV"], 5
    omics = {m: torch.randn(B, 128) for m in MODS}
    fps = torch.randint(0, 2, (B, 2048)).float()
    pair = torch.rand(B, PAIR_DIM)

    print(f"{'configuration':<30}{'out':>8}{'params':>12}{'vs base':>10}")
    print("-" * 60)
    for head in ("mlp", "bilinear"):
        torch.manual_seed(42)
        base = MoGraphDRPAligned(MODS, "fingerprint", "concat", head)
        m = PairGraphDRP(MODS, "fingerprint", "concat", head, pair_dim=PAIR_DIM)
        m.eval()
        with torch.no_grad():
            out = m(omics, fps, pair=pair)
        assert out.shape == (B,) and torch.isfinite(out).all()
        n_base = sum(p.numel() for p in base.parameters())
        n_pair = sum(p.numel() for p in m.parameters())
        # The only new parameters are the extra columns of the first Linear.
        assert n_pair - n_base == PAIR_DIM * m.predictor[0].out_features
        print(f"{'concat + ' + head + ' + pair':<30}{str(tuple(out.shape)):>8}"
              f"{n_pair:>12,}{n_pair - n_base:>+10,}")

    # The prediction must depend on the pair vector with the cell line and the
    # drug held fixed; otherwise the new input is decorative.
    torch.manual_seed(42)
    m = PairGraphDRP(MODS, "fingerprint", "concat", "mlp", pair_dim=PAIR_DIM)
    m.eval()
    with torch.no_grad():
        changed = (m(omics, fps, pair=pair) - m(omics, fps, pair=1 - pair)).abs().mean()
    assert changed > 1e-4, "output ignores the pair vector"
    print(f"\nmean |change| when only the pair vector changes: {changed:.4f}")

    # Squared, not summed: a plain sum over the batch is constant through the
    # predictor's BatchNorm (see the base module's own gradient check).
    m.train()
    m(omics, fps, pair=pair).pow(2).sum().backward()
    pair_columns = m.predictor[0].weight.grad[:, -PAIR_DIM:]
    assert pair_columns.abs().sum() > 0, "no gradient reached the pair columns"
    for part in ("fusion", "drug_encoder", "predictor"):
        norm = sum(float(p.grad.norm()) for n, p in m.named_parameters()
                   if n.startswith(part) and p.grad is not None)
        assert norm > 0, f"no gradient reached {part}"
        print(f"  grad norm {part:<14} {norm:.4f}")
    print(f"  grad norm {'pair columns':<14} {float(pair_columns.norm()):.4f}")

    try:
        m(omics, fps)
    except ValueError:
        pass
    else:
        raise AssertionError("a missing pair vector must fail, not be silently skipped")

    print("\nAll PairGraphDRP checks passed.")
