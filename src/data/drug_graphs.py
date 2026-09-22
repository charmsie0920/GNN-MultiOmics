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

import torch

from src.data.experiment_utils import DRUG_SMILES_PATH

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

    print("\nAll molecular-graph shape checks passed.")
