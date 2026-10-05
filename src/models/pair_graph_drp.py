"""PairGraphDRP: the frozen MoGraphDRP-aligned base plus a pair-specific input.

The base (`src/models/mographdrp_aligned.py`) encodes the cell line and the
drug separately and combines the two vectors. Nothing in it sees how *this*
cell line's mutated proteins relate to *this* drug's target proteins. This
model adds a per-pair vector carrying that relation and concatenates it with
the base's interaction vector before the predictor.

The pair vector comes from one of two sources:

    features    hand-built target-vs-mutation features (Phase 3 gate,
                `src/data/pair_features.py`)
    attention   `PairAttention` (Phase 4 stage A): the drug's target proteins
                attend over the cell line's mutated proteins, each protein a
                learned embedding. No message passing, so not yet a GNN; the
                PPI layers of stage B will replace `protein_states`.

The base is subclassed, not copied: encoders, fusion, head and predictor
layout all come from `MoGraphDRPAligned`, whose file and defaults are frozen.
Only the first predictor layer is widened to accept the extra columns.
"""

from __future__ import annotations

from typing import Dict, NamedTuple, Optional, Tuple, Union

import torch
import torch.nn as nn
import torch.nn.functional as F

from src.data.pair_features import PairSets
from src.models.mographdrp_aligned import MoGraphDRPAligned

TOKEN_DIM = 64          # protein embedding width; small for the 4 GB GPU
ATTENTION_HEADS = 4


class PairAttentionWeights(NamedTuple):
    """Attention of one batch, with what each query and key position is."""

    weights: torch.Tensor      # (B, heads, targets, 1 + mutations), each row sums to 1
    query_index: torch.Tensor  # (B, targets), protein index of each query
    query_mask: torch.Tensor   # (B, targets), False on padding
    key_index: torch.Tensor    # (B, 1 + mutations), protein index of each key
    key_mask: torch.Tensor     # (B, 1 + mutations), False on padding


class PairAttention(nn.Module):
    """A drug's target proteins attend over a cell line's mutated proteins.

    Every protein is one learned token. Queries are the drug's targets, keys
    and values the cell line's mutated proteins; the attended targets are
    mean-pooled to one vector per pair.

    Two learned "none" tokens, stored after the proteins in the embedding
    table, keep every row defined:

    - `no_mutation` is a key for **every** cell line, not only those without
      a mutation. Softmax weights always sum to 1, so without it a target
      with no relevant mutation would have to attend to an irrelevant one.
    - `no_target` is the single query of a drug with no known target.
    """

    def __init__(self, sets: PairSets, dim: int = TOKEN_DIM, heads: int = ATTENTION_HEADS):
        super().__init__()
        self.no_target, self.no_mutation = sets.n_proteins, sets.n_proteins + 1
        self.embedding = nn.Embedding(sets.n_proteins + 2, dim)
        self.attention = nn.MultiheadAttention(dim, heads, batch_first=True)
        self.out_dim = dim

        target_index = torch.from_numpy(sets.target_index).clone()
        target_mask = torch.from_numpy(sets.target_mask).clone()
        no_target = ~target_mask.any(dim=1)
        target_index[no_target, 0] = self.no_target
        target_mask[no_target, 0] = True

        n_cells = len(sets.mutated_index)
        mutated_index = torch.cat([
            torch.full((n_cells, 1), self.no_mutation, dtype=torch.long),
            torch.from_numpy(sets.mutated_index),
        ], dim=1)
        mutated_mask = torch.cat([
            torch.ones(n_cells, 1, dtype=torch.bool),
            torch.from_numpy(sets.mutated_mask),
        ], dim=1)

        self.register_buffer("target_index", target_index, persistent=False)
        self.register_buffer("target_mask", target_mask, persistent=False)
        self.register_buffer("mutated_index", mutated_index, persistent=False)
        self.register_buffer("mutated_mask", mutated_mask, persistent=False)

    def protein_states(self) -> torch.Tensor:
        """-> (n_proteins + 2, dim), the token of every protein and "none" token.

        Stage A: the raw embeddings. Stage B replaces this with the embeddings
        after message passing over the PPI edges; nothing else changes.
        """
        return self.embedding.weight

    def forward(
        self,
        cell_codes: torch.Tensor,
        drug_codes: torch.Tensor,
        return_weights: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, PairAttentionWeights]]:
        """-> (B, dim), and the attention weights when `return_weights`."""
        states = self.protein_states()
        query_index, query_mask = self.target_index[drug_codes], self.target_mask[drug_codes]
        key_index, key_mask = self.mutated_index[cell_codes], self.mutated_mask[cell_codes]
        # `F.embedding`, not `states[index]`: the same lookup, but the indexing
        # form's backward pass is about 9x slower on CUDA (7.6 vs 0.8 ms a step).
        keys = F.embedding(key_index, states)
        attended, weights = self.attention(
            F.embedding(query_index, states), keys, keys, key_padding_mask=~key_mask,
            need_weights=return_weights, average_attn_weights=False,
        )
        real = query_mask.unsqueeze(-1).to(attended.dtype)
        pair = (attended * real).sum(dim=1) / real.sum(dim=1)
        if not return_weights:
            return pair
        return pair, PairAttentionWeights(weights, query_index, query_mask, key_index, key_mask)


class PairGraphDRP(MoGraphDRPAligned):
    """`MoGraphDRPAligned` with a per-pair vector as extra input to the predictor.

    Takes every argument of the base, by keyword, plus exactly one of:

    - `pair_dim`: the pair vector is given to `forward` as `pair_dim` features.
    - `pair_sets`: the pair vector is computed by `PairAttention`; `forward`
      is given each row's cell-line code instead, and uses `drug_codes` for
      the targets.
    """

    def __init__(self, *args, pair_dim: Optional[int] = None,
                 pair_sets: Optional[PairSets] = None, **kwargs):
        super().__init__(*args, **kwargs)
        if (pair_dim is None) == (pair_sets is None):
            raise ValueError("give exactly one of pair_dim (features) or pair_sets (attention)")
        if pair_sets is None and pair_dim < 1:
            raise ValueError(f"pair_dim must be positive, got {pair_dim}")
        self.pair_dim = pair_dim if pair_sets is None else TOKEN_DIM
        first = self.predictor[0]
        self.predictor[0] = nn.Linear(first.in_features + self.pair_dim, first.out_features)
        self.pair_attention = None if pair_sets is None else PairAttention(pair_sets, self.pair_dim)

    def forward(
        self,
        omics: Dict[str, torch.Tensor],
        fingerprints: Optional[torch.Tensor] = None,
        drug_codes: Optional[torch.Tensor] = None,
        pair: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """`pair` is the feature matrix, or the cell-line codes with attention."""
        if pair is None:
            raise ValueError("PairGraphDRP needs the per-pair vector")
        if self.pair_attention is not None:
            pair = self.pair_attention(pair, drug_codes)
        joint = self.interaction_vector(omics, fingerprints, drug_codes)
        return self.predictor(torch.cat([joint, pair], dim=-1)).squeeze(-1)

    @torch.no_grad()
    def attention_weights(self, cell_codes: torch.Tensor,
                          drug_codes: torch.Tensor) -> PairAttentionWeights:
        """Which mutated proteins each target attends to, for the given pairs."""
        if self.pair_attention is None:
            raise ValueError("this model has no attention module (pair features only)")
        return self.pair_attention(cell_codes, drug_codes, return_weights=True)[1]


if __name__ == "__main__":
    import numpy as np

    from src.data.pair_features import pad_index_sets

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

    # --- attention ---------------------------------------------------------
    # A toy graph of 40 proteins. Drug 1 has no target and cell line 0 has no
    # mutation, so every combination of empty and non-empty sets is in a batch.
    N_PROTEINS = 40
    targets = [[3, 7], [], [11], [3, 7, 11, 12, 13, 14, 15]]
    mutated = [[], [3, 20], [11], [21, 22, 23, 24, 25, 26], [7, 11, 30]]
    as_sets = lambda groups: pad_index_sets([np.array(g, dtype=np.int64) for g in groups])  # noqa: E731
    sets = PairSets(np.zeros(0, np.int64), np.zeros(0, np.int64),
                    *as_sets(targets), *as_sets(mutated), n_proteins=N_PROTEINS)
    drugs = torch.arange(B) % len(targets)
    cells = torch.arange(B) // len(targets) % len(mutated)

    print(f"\n{'configuration':<30}{'out':>8}{'params':>12}{'vs base':>10}")
    print("-" * 60)
    for head in ("mlp", "bilinear"):
        torch.manual_seed(42)
        base = MoGraphDRPAligned(MODS, "fingerprint", "concat", head)
        m = PairGraphDRP(MODS, "fingerprint", "concat", head, pair_sets=sets)
        n_base = sum(p.numel() for p in base.parameters())
        n_pair = sum(p.numel() for p in m.parameters())
        n_attention = sum(p.numel() for p in m.pair_attention.parameters())
        assert n_pair - n_base == n_attention + TOKEN_DIM * m.predictor[0].out_features
        for mode in ("train", "eval"):
            getattr(m, mode)()
            with torch.no_grad():
                out = m(omics, fps, drugs, cells)
            assert out.shape == (B,) and torch.isfinite(out).all(), f"{head} {mode}: bad output"
        print(f"{'concat + ' + head + ' + attention':<30}{str(tuple(out.shape)):>8}"
              f"{n_pair:>12,}{n_pair - n_base:>+10,}")

    torch.manual_seed(42)
    m = PairGraphDRP(MODS, "fingerprint", "concat", "mlp", pair_sets=sets)
    m.eval()

    # Empty sets: a drug with no target, a cell line with no mutation, and both.
    empty_drug, empty_cell = torch.tensor([1, 0, 1]), torch.tensor([1, 0, 0])
    with torch.no_grad():
        vectors = m.pair_attention(empty_cell, empty_drug)
    assert torch.isfinite(vectors).all(), "an empty target or mutation set gave NaN"

    w = m.attention_weights(cells, drugs)
    assert w.weights.shape == (B, ATTENTION_HEADS, 7, 1 + 6)
    assert torch.allclose(w.weights.sum(dim=-1), torch.ones(B, ATTENTION_HEADS, 7), atol=1e-5)
    assert (w.weights.masked_select(~w.key_mask[:, None, None, :]) == 0).all(), \
        "attention on a padded key"
    assert (w.key_index[:, 0] == m.pair_attention.no_mutation).all() and w.key_mask[:, 0].all()
    assert (w.query_index[drugs == 1, 0] == m.pair_attention.no_target).all()
    assert (w.query_mask.sum(dim=1) == torch.tensor([2, 1, 1, 7])[drugs]).all()
    print(f"\nattention weights {tuple(w.weights.shape)}: rows sum to 1, zero on padding")

    # The targets must matter beyond the fingerprint, and the mutations beyond
    # the omics: change only the drug code (the fingerprint path ignores it),
    # then only the cell-line code, with every other input held fixed.
    with torch.no_grad():
        ref = m(omics, fps, drugs, cells)
        by_target = (m(omics, fps, (drugs + 2) % len(targets), cells) - ref).abs().mean()
        by_mutation = (m(omics, fps, drugs, (cells + 1) % len(mutated)) - ref).abs().mean()
    assert by_target > 1e-4, "output ignores the target set"
    assert by_mutation > 1e-4, "output ignores the mutated set"
    print(f"mean |change| when only the target set changes:  {by_target:.4f}")
    print(f"mean |change| when only the mutated set changes: {by_mutation:.4f}")

    # Gradients reach the proteins of the batch's pairs and both "none" tokens,
    # and no others. (The batch holds every drug but not cell line 4.)
    m.train()
    m(omics, fps, drugs, cells).pow(2).sum().backward()
    grad = m.pair_attention.embedding.weight.grad.abs().sum(dim=1)
    batch_sets = [targets[d] for d in set(drugs.tolist())] + [mutated[c] for c in set(cells.tolist())]
    in_use = sorted({p for group in batch_sets for p in group}
                    | {m.pair_attention.no_target, m.pair_attention.no_mutation})
    assert 30 not in in_use, "cell line 4 was meant to be left out of the batch"
    assert (grad[in_use] > 0).all(), "no gradient reached a protein that is in use"
    unused = torch.ones(N_PROTEINS + 2, dtype=torch.bool)
    unused[in_use] = False
    assert (grad[unused] == 0).all(), "gradient reached a protein no pair refers to"
    for part in ("pair_attention.embedding", "pair_attention.attention", "predictor"):
        norm = sum(float(p.grad.norm()) for n, p in m.named_parameters()
                   if n.startswith(part) and p.grad is not None)
        assert norm > 0, f"no gradient reached {part}"
        print(f"  grad norm {part:<26} {norm:.4f}")
    print(f"  embedding rows with a gradient: {int((grad > 0).sum())} of {N_PROTEINS + 2}")

    try:
        PairGraphDRP(MODS, "fingerprint", "concat", "mlp", pair_dim=PAIR_DIM, pair_sets=sets)
    except ValueError:
        pass
    else:
        raise AssertionError("features and attention together must be refused")

    print("\nAll PairGraphDRP checks passed.")
