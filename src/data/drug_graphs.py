"""Atom-level molecular graph featurization for the drug branch.

The experiment matrix currently represents a drug either as a one-hot identity
or as a 2048-bit Morgan fingerprint (`build_morgan_fingerprints` in
`experiment_utils.py`). Both discard molecular topology: the fingerprint
enumerates which substructures are present but not how they connect. This
module adds the third representation used by the MoGraphDRP benchmark
(their section 2.2.2) -- each drug as a graph whose nodes are atoms and whose
edges are bonds -- so a GCN can learn over the structure directly.

Deliberately depends on `rdkit` and `torch` only, not on
`torch_geometric`. The 498 resolvable drugs are collated into a single
batched graph once (see the `build_drug_graph_codes`/collate step) rather than
streamed through a PyG `DataLoader`, so nothing here needs PyG's `Data`
container; only the model stage does, for `GCNConv`/`global_mean_pool`.

**Population parity matters more than the featurization details.** This module
reads the same `DRUG_SMILES_PATH` cache and applies the same
drop-the-unresolvable rule as `build_morgan_fingerprints`, so the molecular
graph arm is evaluated on exactly the pairs the fingerprint arm sees
(111,799). If the two drug sets ever diverge, the graph-vs-fingerprint
comparison is confounded by population and the result is worthless -- the
next stage asserts they match rather than assuming it.
"""

from __future__ import annotations

import csv
from typing import Dict, List, NamedTuple, Sequence, Tuple

import numpy as np
import pandas as pd
import torch

from src.data.experiment_utils import (
    COL_CELL_LINE,
    COL_DRUG,
    COL_TARGET,
    DRUG_SMILES_PATH,
    DTYPE,
    cell_row_indices,
)

# --- atom feature vocabulary --------------------------------------------------
# Fixed vocabularies (rather than ones inferred from whichever drugs happen to
# resolve) keep the feature width identical across runs and populations, so a
# checkpoint trained on one drug set stays loadable against another. Anything
# outside a vocabulary falls into that vocabulary's trailing "other" bucket.
ATOM_ELEMENTS: Tuple[str, ...] = (
    "C", "N", "O", "S", "F", "Cl", "Br", "I",
    "P", "B", "Si", "Se", "Na", "Pt", "As", "Zn",
)
ATOM_DEGREES: Tuple[int, ...] = (0, 1, 2, 3, 4, 5)
ATOM_NUM_HS: Tuple[int, ...] = (0, 1, 2, 3, 4)
ATOM_HYBRIDIZATIONS: Tuple[str, ...] = ("SP", "SP2", "SP3", "SP3D", "SP3D2")

# element + degree + n_hydrogens + hybridization one-hots (each with an
# "other" bucket), then formal charge, aromaticity, ring membership.
ATOM_FEATURE_DIM = (
    (len(ATOM_ELEMENTS) + 1)
    + (len(ATOM_DEGREES) + 1)
    + (len(ATOM_NUM_HS) + 1)
    + (len(ATOM_HYBRIDIZATIONS) + 1)
    + 3
)


class MolGraph(NamedTuple):
    """One drug's molecular graph.

    `x` is (n_atoms, ATOM_FEATURE_DIM) float32; `edge_index` is (2, n_bonds*2)
    int64 in COO form, with both directions of every bond present since
    chemical bonds are undirected. A single-atom molecule (e.g. a bare metal
    salt fragment) yields a well-formed (2, 0) `edge_index` rather than an
    error -- message passing over it is a no-op, which is correct.
    """

    x: torch.Tensor
    edge_index: torch.Tensor


def _one_hot(value, vocabulary: Sequence) -> List[float]:
    """One-hot `value` over `vocabulary`, with a trailing catch-all bucket."""
    encoding = [0.0] * (len(vocabulary) + 1)
    try:
        encoding[vocabulary.index(value)] = 1.0
    except ValueError:
        encoding[-1] = 1.0
    return encoding


def atom_features(atom) -> List[float]:
    """Featurize one RDKit atom into a fixed-width `ATOM_FEATURE_DIM` vector.

    Mirrors the descriptor set the benchmark names (element type, bond count,
    hydrogen count, aromaticity, "and other descriptive properties"), plus
    hybridization and ring membership, which are standard in this family of
    models and cost nothing to compute.
    """
    return (
        _one_hot(atom.GetSymbol(), ATOM_ELEMENTS)
        + _one_hot(atom.GetDegree(), ATOM_DEGREES)
        + _one_hot(atom.GetTotalNumHs(), ATOM_NUM_HS)
        + _one_hot(str(atom.GetHybridization()), ATOM_HYBRIDIZATIONS)
        + [
            float(atom.GetFormalCharge()),
            float(atom.GetIsAromatic()),
            float(atom.IsInRing()),
        ]
    )


def mol_to_graph(mol) -> MolGraph:
    """Convert an RDKit molecule to node features and a COO bond edge index."""
    x = torch.tensor(
        [atom_features(atom) for atom in mol.GetAtoms()], dtype=torch.float32
    )

    sources: List[int] = []
    targets: List[int] = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        sources += [i, j]
        targets += [j, i]

    edge_index = torch.tensor([sources, targets], dtype=torch.long)
    if edge_index.numel() == 0:  # bond-free molecule: keep the (2, 0) shape
        edge_index = torch.zeros((2, 0), dtype=torch.long)

    return MolGraph(x=x, edge_index=edge_index)


def build_drug_graphs() -> Dict[str, MolGraph]:
    """Parse cached SMILES into per-drug molecular graphs, keyed by GDSC drug_id.

    Drugs with missing/unparseable SMILES are absent from the returned dict
    (not zero-filled), matching `build_morgan_fingerprints` so that callers
    drop the same pairs. Molecules that parse to zero atoms are dropped too
    and counted separately: RDKit returns an empty molecule rather than `None`
    for some degenerate inputs, which the fingerprint path would silently
    featurize as an all-zero vector.
    """
    from rdkit import Chem

    if not DRUG_SMILES_PATH.exists():
        raise FileNotFoundError(
            f"Missing {DRUG_SMILES_PATH}. Run src/data/fetch_drug_smiles.py first."
        )

    with DRUG_SMILES_PATH.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    graphs: Dict[str, MolGraph] = {}
    unparseable = 0
    empty = 0
    for row in rows:
        smiles = row["canonical_smiles"]
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            unparseable += 1
            continue
        if mol.GetNumAtoms() == 0:
            empty += 1
            continue
        graphs[row["drug_id"]] = mol_to_graph(mol)

    total_atoms = sum(g.x.shape[0] for g in graphs.values())
    total_edges = sum(g.edge_index.shape[1] for g in graphs.values())
    print(
        f"[molgraphs] {len(graphs)} / {len(rows)} drug_ids resolved to a molecular "
        f"graph ({unparseable} missing/unparseable SMILES dropped, {empty} zero-atom)"
    )
    print(
        f"[molgraphs] {total_atoms} atoms, {total_edges} directed bond edges, "
        f"{ATOM_FEATURE_DIM}-dim atom features"
    )
    return graphs


# --- population parity with the fingerprint arm -------------------------------
def assert_fingerprint_population_parity(graphs: Dict[str, MolGraph]) -> None:
    """Fail loudly if the molecular graph and fingerprint arms cover different drugs.

    The whole point of the molecular graph arm is a controlled comparison
    against the fingerprint arm (E10 in docs/results.md) with only the drug
    encoder changed. That holds only if both arms are fit and scored on the
    same rows, which in turn requires the same resolvable drug set. RDKit can
    in principle resolve a SMILES that `build_morgan_fingerprints` keeps but
    this module drops (a molecule that parses to zero atoms), so this is
    checked rather than assumed.
    """
    from src.data.experiment_utils import build_morgan_fingerprints

    fingerprint_ids = set(build_morgan_fingerprints())
    graph_ids = set(graphs)
    if graph_ids != fingerprint_ids:
        graph_only = sorted(graph_ids - fingerprint_ids)
        fingerprint_only = sorted(fingerprint_ids - graph_ids)
        raise AssertionError(
            "Molecular graph and fingerprint arms cover different drugs, so their "
            "results would not be comparable.\n"
            f"  graph-only ({len(graph_only)}): {graph_only[:10]}\n"
            f"  fingerprint-only ({len(fingerprint_only)}): {fingerprint_only[:10]}"
        )
    print(f"[molgraphs] population parity OK: {len(graph_ids)} drugs, identical to fingerprint arm")


# --- batched graph over the whole drug vocabulary ------------------------------
class BatchedMolGraphs(NamedTuple):
    """All drugs' molecular graphs concatenated into one disjoint graph.

    `x` is (total_atoms, ATOM_FEATURE_DIM), `edge_index` is (2, total_edges)
    with per-drug atom offsets already applied, and `batch` is
    (total_atoms,) mapping each atom to its drug's code, which is what
    `global_mean_pool` needs to pool back to one vector per drug.

    There are only ~498 distinct drugs but ~111,799 pairs, so the model
    encodes this batch **once per forward pass** and indexes the resulting
    (n_drugs, out_dim) table by drug code, instead of re-encoding the same
    molecule once per pair. That keeps the existing dense `TensorDataset`
    loader intact -- the drug column becomes an int64 code rather than a
    2048-bit row -- and avoids PyG's variable-size batching entirely.
    """

    x: torch.Tensor
    edge_index: torch.Tensor
    batch: torch.Tensor
    n_graphs: int


def collate_drug_graphs(
    graphs: Dict[str, MolGraph], drug_order: Sequence[str]
) -> BatchedMolGraphs:
    """Concatenate `graphs` into one disjoint graph, ordered by `drug_order`.

    `drug_order` must be the same drug-code ordering the pair tensors use, so
    that `batch == code` holds and the pooled table can be indexed directly.
    """
    xs: List[torch.Tensor] = []
    edge_indices: List[torch.Tensor] = []
    batches: List[torch.Tensor] = []

    atom_offset = 0
    for code, drug_id in enumerate(drug_order):
        graph = graphs[drug_id]
        n_atoms = graph.x.shape[0]
        xs.append(graph.x)
        edge_indices.append(graph.edge_index + atom_offset)
        batches.append(torch.full((n_atoms,), code, dtype=torch.long))
        atom_offset += n_atoms

    batched = BatchedMolGraphs(
        x=torch.cat(xs, dim=0),
        edge_index=torch.cat(edge_indices, dim=1),
        batch=torch.cat(batches, dim=0),
        n_graphs=len(drug_order),
    )
    print(
        f"[molgraphs] batched {batched.n_graphs} drugs -> "
        f"x={tuple(batched.x.shape)}, edge_index={tuple(batched.edge_index.shape)}"
    )
    return batched


# --- pair tensors (molecular graph drug arm) -----------------------------------
def build_graph_pair_tensors(
    omics: Dict[str, np.ndarray],
    cell_ids: pd.Index,
    y: pd.DataFrame,
    graphs: Dict[str, MolGraph],
) -> Tuple[Dict[str, np.ndarray], np.ndarray, np.ndarray, np.ndarray, BatchedMolGraphs, pd.DataFrame]:
    """Molecular-graph counterpart of `experiment_utils.build_pair_tensors`.

    Differs in exactly one respect: the drug block is an int64 *code* per pair
    rather than a dense feature row, with the structure itself carried once in
    the returned `BatchedMolGraphs`. Pairs whose drug has no resolvable graph
    are dropped, matching the fingerprint arm's rule.

    Returns (gathered_omics, drug_codes, target, groups, batched_graphs, y_used).
    """
    n_before = len(y)
    y_used = y[y[COL_DRUG].isin(graphs.keys())].reset_index(drop=True)
    print(
        f"[drug_features] graph mode: {n_before - len(y_used)} pairs dropped "
        f"(drug has no resolvable molecular graph) -> {len(y_used)} pairs remain"
    )

    # `sort=True` makes the code assignment deterministic (alphabetical by
    # drug_id) rather than dependent on row order. The code is only a lookup
    # key into the batched graph table -- it never reaches the model as a
    # feature, so the drug is represented by its structure alone.
    drug_codes, drug_levels = pd.factorize(y_used[COL_DRUG], sort=True)
    batched = collate_drug_graphs(graphs, list(drug_levels))

    rows = cell_row_indices(cell_ids, y_used)
    gathered = {key: arr[rows] for key, arr in omics.items()}
    target = y_used[COL_TARGET].to_numpy(dtype=DTYPE)
    groups = y_used[COL_CELL_LINE].to_numpy()

    omics_bytes = sum(a.nbytes for a in gathered.values())
    print(
        f"[design]  {len(y_used)} pairs  |  omics: "
        f"{' + '.join(f'{k}({v.shape[1]})' for k, v in gathered.items())}"
        f"  |  drug/graph: {batched.n_graphs} molecules  =  "
        f"{(omics_bytes + batched.x.nbytes) / 1024**2:.1f} MB"
    )
    return gathered, drug_codes.astype(np.int64), target, groups, batched, y_used


if __name__ == "__main__":
    from rdkit import Chem

    # Hand-checkable molecules spanning the edge cases: a ring + substituents,
    # a bare aromatic ring, a platinum complex whose fragments are disconnected,
    # and a single bond-free atom.
    SAMPLES = {
        "aspirin": ("CC(=O)Oc1ccccc1C(=O)O", 13, 13),
        "benzene": ("c1ccccc1", 6, 6),
        "cisplatin": ("N.N.Cl[Pt]Cl", 5, 2),
        "lone_atom": ("[Na+]", 1, 0),
    }

    print(f"ATOM_FEATURE_DIM = {ATOM_FEATURE_DIM}\n")
    for name, (smiles, expected_atoms, expected_bonds) in SAMPLES.items():
        mol = Chem.MolFromSmiles(smiles)
        graph = mol_to_graph(mol)
        n_atoms, n_edges = graph.x.shape[0], graph.edge_index.shape[1]
        print(
            f"  {name:<10} {smiles:<22} x={tuple(graph.x.shape)} "
            f"edge_index={tuple(graph.edge_index.shape)}"
        )
        assert n_atoms == expected_atoms, f"{name}: {n_atoms} atoms != {expected_atoms}"
        assert n_edges == expected_bonds * 2, f"{name}: {n_edges} != {expected_bonds * 2}"
        assert graph.x.shape[1] == ATOM_FEATURE_DIM
        assert graph.edge_index.dtype == torch.long
        if n_edges:  # every edge must index a real atom
            assert int(graph.edge_index.max()) < n_atoms

    # Collation: the batched graph must stay block-diagonal -- every edge has
    # to stay inside its own molecule, or message passing would leak structure
    # between unrelated drugs.
    graphs = {
        name: mol_to_graph(Chem.MolFromSmiles(smiles))
        for name, (smiles, _, _) in SAMPLES.items()
    }
    order = list(SAMPLES)
    batched = collate_drug_graphs(graphs, order)

    total_atoms = sum(graphs[d].x.shape[0] for d in order)
    total_edges = sum(graphs[d].edge_index.shape[1] for d in order)
    assert batched.x.shape == (total_atoms, ATOM_FEATURE_DIM)
    assert batched.edge_index.shape == (2, total_edges)
    assert batched.batch.shape == (total_atoms,)
    assert batched.n_graphs == len(order)
    # Atoms appear in drug-code order, so `batch` is non-decreasing and covers
    # every code exactly as many times as that molecule has atoms.
    assert torch.equal(batched.batch, batched.batch.sort().values)
    assert batched.batch.bincount().tolist() == [graphs[d].x.shape[0] for d in order]
    # No edge may cross a molecule boundary.
    assert torch.equal(
        batched.batch[batched.edge_index[0]], batched.batch[batched.edge_index[1]]
    )

    print("\nAll molecular-graph shape checks passed.")
