"""Hand-built pair features: how a cell line's mutations sit relative to a drug's targets.

The Phase 3 gate for the pair-specific idea. Every model so far describes the
cell line and the drug separately; nothing tells it whether *this* cell line's
mutated proteins are the proteins *this* drug acts on. These features state
that directly, from the heterogeneous graph (`src/graph/hetero_graph.pt`):

    has_target    the drug has at least one known target protein
    has_mutation  the cell line has at least one driver-mutation edge
    direct        a target protein is itself mutated (min hops == 0)
    proximity     1 / (1 + min PPI hops between the mutated set and the target
                  set); 0 when either set is empty or no path exists
    near_count    log1p(number of mutated proteins within 1 hop of a target)

They are not a GNN: no parameters, no message passing. If they carry no
signal the base model lacks, a learned module over the same edges is unlikely
to find one, which is what the gate decides before Phase 4 is built.

Phase 4 reads the same two edge types as sets instead of features:
`build_pair_sets` gives each drug's target proteins and each cell line's
mutated proteins as padded index tables, for the attention module in
`src/models/pair_graph_drp.py`.

**Not leakage.** Mutation and target edges are inputs, not labels. A held-out
cell line's mutations are part of its description in the same way its
expression profile is, so its edges being in the graph tells the model nothing
about its measured ln(IC50).

Run from the repository root to print the coverage table:
    python src/data/pair_features.py
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from scipy.sparse.csgraph import shortest_path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.experiment_utils import (  # noqa: E402
    COL_CELL_LINE,
    COL_DRUG,
    DTYPE,
    GE_KEY,
    MUT_CNV_KEY,
    PROTEOMICS_KEY,
    build_drug_features,
    grouped_split,
    load_omics_subset,
    load_targets,
)

GRAPH_PATH = Path("src/graph/hetero_graph.pt")

PPI = ("protein", "interacts_with", "protein")
TARGETS = ("drug", "targets", "protein")
MUTATION = ("cell_line", "has_mutation", "protein")

FEATURE_NAMES = ["has_target", "has_mutation", "direct", "proximity", "near_count"]


class PairFeatures(NamedTuple):
    """Per-pair features and the row subsets the gate is reported on."""

    matrix: np.ndarray               # (n_pairs, len(FEATURE_NAMES)), float32
    names: List[str]
    min_hops: np.ndarray             # (n_pairs,), inf when undefined or unreachable
    subsets: Dict[str, np.ndarray]   # name -> boolean row mask


class PairGraph(NamedTuple):
    """The graph, with the pair rows' cell lines and drugs matched to its nodes."""

    data: object                         # the HeteroData graph
    cell_codes: np.ndarray               # (n_pairs,), index into `cell_levels`
    cell_levels: pd.Index                # sorted cell-line IDs of the pairs
    drug_codes: np.ndarray               # (n_pairs,), index into `drug_levels`
    drug_levels: pd.Index                # sorted drug IDs of the pairs
    targets_of: Dict[str, np.ndarray]    # drug ID -> target protein indices
    mutated_of: Dict[str, np.ndarray]    # cell-line ID -> mutated protein indices


class PairSets(NamedTuple):
    """Each drug's target proteins and each cell line's mutated proteins, padded.

    What the Phase 4 attention module reads. A pair row picks its drug's row of
    the target table and its cell line's row of the mutation table.
    """

    cell_codes: np.ndarray       # (n_pairs,), int64, row of `mutated_index`
    drug_codes: np.ndarray       # (n_pairs,), int64, row of `target_index`
    target_index: np.ndarray     # (n_drugs, most targets), protein node index
    target_mask: np.ndarray      # same shape, True on real entries
    mutated_index: np.ndarray    # (n_cell_lines, most mutations), protein node index
    mutated_mask: np.ndarray     # same shape, True on real entries
    n_proteins: int


def _edges_by_source(edge_index: torch.Tensor, source_ids: List[str]) -> Dict[str, np.ndarray]:
    """`{source node id: protein indices}` for one (source -> protein) edge type."""
    grouped: Dict[str, List[int]] = {}
    for source, protein in edge_index.t().tolist():
        grouped.setdefault(source_ids[source], []).append(protein)
    return {k: np.unique(v) for k, v in grouped.items()}


def _load_pair_graph(y_used: pd.DataFrame, graph_path: Path) -> PairGraph:
    """Load the graph and match the pair rows' cell lines and drugs to its nodes.

    Matched by node ID, never by position: the graph holds 498 drugs and the
    pairs use 240 of them.
    """
    data = torch.load(graph_path, weights_only=False)
    cell_ids = [str(c) for c in data["cell_line"].node_ids]
    drug_ids = [str(d) for d in data["drug"].node_ids]

    cell_codes, cell_levels = pd.factorize(y_used[COL_CELL_LINE].astype(str), sort=True)
    drug_codes, drug_levels = pd.factorize(y_used[COL_DRUG].astype(str), sort=True)
    missing_cells = set(cell_levels) - set(cell_ids)
    missing_drugs = set(drug_levels) - set(drug_ids)
    if missing_cells or missing_drugs:
        raise AssertionError(
            f"{len(missing_cells)} cell lines and {len(missing_drugs)} drugs in the pairs "
            f"are not graph nodes, e.g. {sorted(missing_cells)[:3]} {sorted(missing_drugs)[:3]}"
        )
    return PairGraph(
        data, cell_codes, cell_levels, drug_codes, drug_levels,
        targets_of=_edges_by_source(data[TARGETS].edge_index, drug_ids),
        mutated_of=_edges_by_source(data[MUTATION].edge_index, cell_ids),
    )


def pad_index_sets(sets: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    """Variable-length index sets -> (index, mask), both (len(sets), longest set).

    Entries are left-aligned; padded slots hold index 0 and are False in the mask.
    """
    width = max(1, max(len(s) for s in sets))
    index = np.zeros((len(sets), width), dtype=np.int64)
    mask = np.zeros((len(sets), width), dtype=bool)
    for i, members in enumerate(sets):
        index[i, :len(members)] = members
        mask[i, :len(members)] = True
    return index, mask


def build_pair_sets(y_used: pd.DataFrame, graph_path: Path = GRAPH_PATH) -> PairSets:
    """The protein sets the Phase 4 attention module reads, for the rows of `y_used`.

    `drug_codes` follow the same sorted drug order as the ladder's own drug
    codes; the ladder asserts that before using them.
    """
    graph = _load_pair_graph(y_used, graph_path)
    empty = np.empty(0, dtype=np.int64)
    target_index, target_mask = pad_index_sets(
        [graph.targets_of.get(d, empty) for d in graph.drug_levels])
    mutated_index, mutated_mask = pad_index_sets(
        [graph.mutated_of.get(c, empty) for c in graph.cell_levels])
    used = np.union1d(target_index[target_mask], mutated_index[mutated_mask])
    print(
        f"[pair]    targets per drug: up to {target_mask.shape[1]}, "
        f"{int((~target_mask.any(axis=1)).sum())} of {len(target_mask)} drugs have none  |  "
        f"mutated proteins per cell line: up to {mutated_mask.shape[1]}, "
        f"{int((~mutated_mask.any(axis=1)).sum())} of {len(mutated_mask)} have none  |  "
        f"{len(used)} distinct proteins in either"
    )
    return PairSets(
        graph.cell_codes.astype(np.int64), graph.drug_codes.astype(np.int64),
        target_index, target_mask, mutated_index, mutated_mask,
        n_proteins=graph.data["protein"].x.shape[0],
    )


def build_pair_features(y_used: pd.DataFrame, graph_path: Path = GRAPH_PATH) -> PairFeatures:
    """Features for every row of `y_used`, in its row order.

    Cell lines and drugs are matched to graph nodes by their IDs, never by
    position: the graph holds 498 drugs and the pairs use 240 of them.
    """
    data, cell_codes, cell_levels, drug_codes, drug_levels, targets_of, mutated_of = \
        _load_pair_graph(y_used, graph_path)
    n_proteins = data["protein"].x.shape[0]

    # Hop counts from every target protein the pairs can reference. STRING
    # lists each interaction in both directions, so the graph is undirected.
    ppi = data[PPI].edge_index.numpy()
    adjacency = sp.csr_matrix(
        (np.ones(ppi.shape[1], dtype=np.int8), (ppi[0], ppi[1])), shape=(n_proteins, n_proteins)
    )
    target_proteins = np.unique(
        np.concatenate([targets_of[d] for d in drug_levels if d in targets_of])
    )
    row_of_target = {p: i for i, p in enumerate(target_proteins)}
    hops = shortest_path(adjacency, unweighted=True, indices=target_proteins)

    # One value per (cell line, drug); the pair rows are gathered from this grid.
    min_hops = np.full((len(cell_levels), len(drug_levels)), np.inf)
    near_count = np.zeros((len(cell_levels), len(drug_levels)))
    for j, drug in enumerate(drug_levels):
        if drug not in targets_of:
            continue
        to_target = hops[[row_of_target[p] for p in targets_of[drug]]].min(axis=0)
        for i, cell in enumerate(cell_levels):
            mutated = mutated_of.get(cell)
            if mutated is None:
                continue
            distances = to_target[mutated]
            min_hops[i, j] = distances.min()
            near_count[i, j] = (distances <= 1).sum()

    has_target = np.array([d in targets_of for d in drug_levels])[drug_codes]
    has_mutation = np.array([c in mutated_of for c in cell_levels])[cell_codes]
    pair_hops = min_hops[cell_codes, drug_codes]
    direct = pair_hops == 0

    matrix = np.stack(
        [
            has_target,
            has_mutation,
            direct,
            1.0 / (1.0 + pair_hops),          # inf -> 0
            np.log1p(near_count[cell_codes, drug_codes]),
        ],
        axis=1,
    ).astype(DTYPE)
    assert matrix.shape == (len(y_used), len(FEATURE_NAMES)) and np.isfinite(matrix).all()

    subsets = {
        "all": np.ones(len(y_used), dtype=bool),
        "has_target": has_target,
        "no_target": ~has_target,
        "direct_hit": direct,
        "hop_1": pair_hops == 1,
        "hop_2plus": np.isfinite(pair_hops) & (pair_hops >= 2),
    }
    print(
        f"[pair]    {len(y_used)} pairs x {matrix.shape[1]} features  |  "
        f"drug has a target: {has_target.mean():.2%}  |  target directly mutated: "
        f"{direct.mean():.2%} ({int(direct.sum())} pairs)"
    )
    return PairFeatures(matrix, list(FEATURE_NAMES), pair_hops, subsets)


def load_pair_rows() -> pd.DataFrame:
    """The ladder's pair rows (`y_used`), without the omics tensors or molecular graphs.

    Same filters and row order as `load_pairs` in `src/final_model/run_ablation.py`,
    for scripts that only need to know which pair each row is.
    """
    _, cell_ids = load_omics_subset((GE_KEY, MUT_CNV_KEY, PROTEOMICS_KEY))
    y_used, _, _, _ = build_drug_features(load_targets(cell_ids), "fingerprint")
    return y_used


if __name__ == "__main__":
    y_used = load_pair_rows()
    pair = build_pair_features(y_used)
    train_idx, val_idx, test_idx = grouped_split(y_used[COL_CELL_LINE].to_numpy())

    print(f"\n{'subset':<14}{'all pairs':>12}{'share':>9}{'train':>9}{'val':>9}{'test':>9}")
    print("-" * 62)
    for name, mask in pair.subsets.items():
        print(f"{name:<14}{int(mask.sum()):>12}{mask.mean():>9.2%}"
              f"{int(mask[train_idx].sum()):>9}{int(mask[val_idx].sum()):>9}"
              f"{int(mask[test_idx].sum()):>9}")

    print(f"\n{'feature':<14}{'min':>9}{'mean':>9}{'max':>9}")
    print("-" * 41)
    for k, name in enumerate(pair.names):
        column = pair.matrix[:, k]
        print(f"{name:<14}{column.min():>9.3f}{column.mean():>9.3f}{column.max():>9.3f}")

    s, m = pair.subsets, dict(zip(pair.names, pair.matrix.T))
    assert not (s["direct_hit"] & ~s["has_target"]).any(), "direct hit without a target"
    assert (m["near_count"][s["direct_hit"]] > 0).all(), "direct hit must count as within 1 hop"
    assert (m["proximity"][s["direct_hit"]] == 1).all()
    assert (m["proximity"][s["no_target"]] == 0).all(), "no-target pairs must carry no proximity"
    assert s["has_target"].sum() + s["no_target"].sum() == len(y_used)
    # A drug's row features must differ between cell lines, or they add nothing
    # to the drug fingerprint the model already has.
    one_drug = y_used[COL_DRUG] == y_used[COL_DRUG][s["direct_hit"]].iloc[0]
    assert len(np.unique(pair.matrix[one_drug.to_numpy()], axis=0)) > 1

    # The padded sets must describe the same graph as the features: a pair is a
    # direct hit exactly when its target set and its mutated set share a protein.
    print()
    sets = build_pair_sets(y_used)
    shared = np.zeros((len(sets.mutated_index), len(sets.target_index)), dtype=bool)
    for j, (targets, real) in enumerate(zip(sets.target_index, sets.target_mask)):
        shared[:, j] = (np.isin(sets.mutated_index, targets[real]) & sets.mutated_mask).any(axis=1)
    assert np.array_equal(shared[sets.cell_codes, sets.drug_codes], s["direct_hit"]), \
        "padded sets disagree with the direct-hit feature"
    drug_has_target = sets.target_mask.any(axis=1)[sets.drug_codes]
    assert np.array_equal(drug_has_target, s["has_target"])
    assert np.array_equal(sets.mutated_mask.any(axis=1)[sets.cell_codes],
                          pair.matrix[:, pair.names.index("has_mutation")] == 1)
    assert sets.target_index.max() < sets.n_proteins and sets.mutated_index.max() < sets.n_proteins
    print("\nAll pair-feature checks passed.")
