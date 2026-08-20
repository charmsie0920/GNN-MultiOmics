"""Build the heterogeneous `cell_line` / `drug` / `protein` graph and save it to disk.

Reads only local files produced by earlier pipeline stages (PCA-compressed
omics, `fetch_drug_smiles.py`, `fetch_string_aliases.py`, and the raw STRING
PPI links) and assembles a `torch_geometric.data.HeteroData` object:

- `cell_line` nodes: PCA-compressed transcriptomics/genomics/proteomics,
  concatenated per cell line.
- `drug` nodes: 2048-bit Morgan fingerprints computed from cached SMILES.
- `protein` nodes: zero-initialized placeholder embeddings (learned later by
  an `nn.Embedding` in the GNN stage) for every protein that survives the
  STRING confidence filter or is targeted by a resolved drug.

Memory-efficiency techniques used, per the graph construction plan:
- STRING's 13.7M-row PPI file is streamed in 1M-row chunks with a confidence
  filter, never materialized whole in memory.
- Node indices are built incrementally in plain Python dicts/lists; only the
  final edge tensors are converted to `torch`.
- The STRING alias file (protein -> gene symbol) is filtered to only the
  ~380 gene symbols actually referenced by GDSC drug targets, not loaded
  wholesale (the full file has ~40M rows).
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from torch_geometric.data import HeteroData

TRANSCRIPTOMICS_PCA_PATH = Path("data/processed/transcriptomics_pca.csv")
GENOMICS_PCA_PATH = Path("data/processed/genomics_pca.csv")
PROTEOMICS_PCA_PATH = Path("data/processed/proteomics_pca.csv")
DRUG_SMILES_PATH = Path("data/raw/pubchem/gdsc_drug_smiles.csv")
GDSC_COMPOUNDS_PATH = Path("data/raw/gdsc/screened_compounds_rel_8.5.csv")
STRING_LINKS_PATH = Path("data/raw/string/9606.protein.links.v12.0.txt.gz")
STRING_ALIASES_PATH = Path("data/raw/string/9606.protein.aliases.v12.0.txt.gz")
OUTPUT_PATH = Path("data/processed/hetero_graph.pt")

# STRING's standard "high confidence" cutoff on combined_score (0-999 scale).
PPI_SCORE_THRESHOLD = 700
MORGAN_RADIUS = 2
MORGAN_FP_SIZE = 2048
PROTEIN_EMBED_DIM = 128
PPI_CHUNKSIZE = 1_000_000


def build_cell_line_nodes() -> Tuple[torch.Tensor, List[str]]:
    """Inner-join the three PCA CSVs on `standard_model_id` and concatenate features."""

    transcriptomics = pd.read_csv(TRANSCRIPTOMICS_PCA_PATH, index_col="standard_model_id")
    genomics = pd.read_csv(GENOMICS_PCA_PATH, index_col="standard_model_id")
    proteomics = pd.read_csv(PROTEOMICS_PCA_PATH, index_col="standard_model_id")

    common_ids = sorted(set(transcriptomics.index) & set(genomics.index) & set(proteomics.index))
    dropped = (set(transcriptomics.index) | set(genomics.index) | set(proteomics.index)) - set(common_ids)
    if dropped:
        print(f"[cell_line] {len(dropped)} cell line IDs did not match across all 3 modalities, dropped: {sorted(dropped)}")

    features = np.concatenate(
        [
            transcriptomics.loc[common_ids].to_numpy(dtype=np.float32),
            genomics.loc[common_ids].to_numpy(dtype=np.float32),
            proteomics.loc[common_ids].to_numpy(dtype=np.float32),
        ],
        axis=1,
    )
    return torch.from_numpy(features), common_ids


def build_drug_nodes() -> Tuple[torch.Tensor, List[str], Dict[str, int]]:
    """Parse cached SMILES into Morgan fingerprints; drop unparseable rows."""

    with DRUG_SMILES_PATH.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=MORGAN_RADIUS, fpSize=MORGAN_FP_SIZE)

    fingerprints: List[np.ndarray] = []
    drug_ids: List[str] = []
    unparseable = 0
    for row in rows:
        smiles = row["canonical_smiles"]
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            unparseable += 1
            continue
        fp = mfpgen.GetFingerprint(mol)
        arr = np.zeros((MORGAN_FP_SIZE,), dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fingerprints.append(arr)
        drug_ids.append(row["drug_id"])

    print(f"[drug] {unparseable} / {len(rows)} rows dropped (missing/unparseable SMILES)")

    features = torch.from_numpy(np.stack(fingerprints, axis=0)) if fingerprints else torch.empty((0, MORGAN_FP_SIZE))
    drug_id_to_index = {drug_id: index for index, drug_id in enumerate(drug_ids)}
    return features, drug_ids, drug_id_to_index


def stream_ppi_edges(score_threshold: int) -> Tuple[torch.Tensor, List[str], Dict[str, int]]:
    """Stream STRING PPI links in chunks, filter by score, and index proteins incrementally.

    STRING's protein.links file already lists both directions of every
    interaction (verified: row count is 13,715,404, matching the documented
    total, and a spot check confirmed A-B and B-A are both present with the
    same score) -- so edges are used as-is, with no manual mirroring needed.
    """

    protein_to_index: Dict[str, int] = {}
    protein_ids: List[str] = []
    src_indices: List[int] = []
    dst_indices: List[int] = []

    def get_index(protein_id: str) -> int:
        index = protein_to_index.get(protein_id)
        if index is None:
            index = len(protein_ids)
            protein_to_index[protein_id] = index
            protein_ids.append(protein_id)
        return index

    total_rows = 0
    kept_rows = 0
    for chunk in pd.read_csv(STRING_LINKS_PATH, sep=" ", chunksize=PPI_CHUNKSIZE):
        total_rows += len(chunk)
        filtered = chunk[chunk["combined_score"] >= score_threshold]
        kept_rows += len(filtered)
        for protein1, protein2 in zip(filtered["protein1"], filtered["protein2"]):
            src_indices.append(get_index(protein1))
            dst_indices.append(get_index(protein2))

    print(f"[interacts_with] {kept_rows} / {total_rows} STRING rows kept at score >= {score_threshold}")

    edge_index = torch.tensor([src_indices, dst_indices], dtype=torch.long)
    return edge_index, protein_ids, protein_to_index


def load_drug_targets() -> Dict[str, List[str]]:
    """Load GDSC's `TARGET` column, keyed by `DRUG_ID`, split into gene-symbol tokens."""

    with GDSC_COMPOUNDS_PATH.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    targets_by_drug_id: Dict[str, List[str]] = {}
    for row in rows:
        raw = row.get("TARGET", "").strip()
        if not raw or raw.lower() == "not defined":
            continue
        symbols = [token.strip() for token in raw.split(",") if token.strip()]
        targets_by_drug_id[row["DRUG_ID"]] = symbols
    return targets_by_drug_id


def build_gene_symbol_to_ensp(referenced_symbols: Set[str]) -> Dict[str, str]:
    """Build a gene_symbol -> ENSP map from the STRING alias file, filtered to referenced symbols.

    The full alias file has ~40M rows; loading it wholesale into memory is
    unnecessary when only ~380 distinct GDSC target symbols need resolving, so
    this does a single streamed pass and only retains rows whose alias is in
    `referenced_symbols`. The first ENSP seen for a symbol is kept (STRING
    lists the same gene symbol from multiple alias sources for the same
    canonical protein in practice; no attempt is made to adjudicate
    conflicting sources beyond first-match).
    """

    symbol_to_ensp: Dict[str, str] = {}
    with gzip.open(STRING_ALIASES_PATH, "rt", encoding="utf-8") as handle:
        next(handle)  # header: #string_protein_id, alias, source
        for line in handle:
            protein_id, alias, _source = line.rstrip("\n").split("\t")
            if alias in referenced_symbols and alias not in symbol_to_ensp:
                symbol_to_ensp[alias] = protein_id
    return symbol_to_ensp


def build_drug_target_edges(
    targets_by_drug_id: Dict[str, List[str]],
    drug_id_to_index: Dict[str, int],
    protein_to_index: Dict[str, int],
    protein_ids: List[str],
) -> torch.Tensor:
    """Map resolved drug targets to protein nodes, adding new protein nodes as needed."""

    referenced_symbols = {symbol for symbols in targets_by_drug_id.values() for symbol in symbols}
    symbol_to_ensp = build_gene_symbol_to_ensp(referenced_symbols)

    resolved_symbols = set(symbol_to_ensp.keys())
    unresolved_symbols = referenced_symbols - resolved_symbols
    print(f"[targets] {len(referenced_symbols)} distinct target symbols referenced by GDSC drugs")
    print(f"[targets] {len(resolved_symbols)} resolved to an ENSP, {len(unresolved_symbols)} unresolved")
    if unresolved_symbols:
        print(f"[targets] unresolved symbols: {sorted(unresolved_symbols)}")

    src_indices: List[int] = []
    dst_indices: List[int] = []
    for drug_id, symbols in targets_by_drug_id.items():
        drug_index = drug_id_to_index.get(drug_id)
        if drug_index is None:
            continue  # drug node was dropped (unparseable/missing SMILES)
        for symbol in symbols:
            ensp = symbol_to_ensp.get(symbol)
            if ensp is None:
                continue
            protein_index = protein_to_index.get(ensp)
            if protein_index is None:
                protein_index = len(protein_ids)
                protein_to_index[ensp] = protein_index
                protein_ids.append(ensp)
            src_indices.append(drug_index)
            dst_indices.append(protein_index)

    return torch.tensor([src_indices, dst_indices], dtype=torch.long)


def tensor_memory_mb(data: HeteroData) -> float:
    """Sum tensor.numel() * tensor.element_size() across every node/edge store, in MB."""

    total_bytes = 0
    for store in list(data.node_stores) + list(data.edge_stores):
        for value in store.values():
            if isinstance(value, torch.Tensor):
                total_bytes += value.numel() * value.element_size()
    return total_bytes / (1024 * 1024)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Build the HeteroData graph from local files and save it to disk."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-threshold", type=int, default=PPI_SCORE_THRESHOLD)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args(argv)

    required = [
        TRANSCRIPTOMICS_PCA_PATH,
        GENOMICS_PCA_PATH,
        PROTEOMICS_PCA_PATH,
        DRUG_SMILES_PATH,
        GDSC_COMPOUNDS_PATH,
        STRING_LINKS_PATH,
        STRING_ALIASES_PATH,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required input files: " + ", ".join(missing))

    data = HeteroData()

    cell_line_x, cell_line_ids = build_cell_line_nodes()
    data["cell_line"].x = cell_line_x
    data["cell_line"].node_ids = cell_line_ids

    drug_x, drug_ids, drug_id_to_index = build_drug_nodes()
    data["drug"].x = drug_x
    data["drug"].node_ids = drug_ids

    ppi_edge_index, protein_ids, protein_to_index = stream_ppi_edges(args.score_threshold)

    targets_by_drug_id = load_drug_targets()
    target_edge_index = build_drug_target_edges(targets_by_drug_id, drug_id_to_index, protein_to_index, protein_ids)

    data["protein"].x = torch.zeros((len(protein_ids), PROTEIN_EMBED_DIM), dtype=torch.float32)
    data["protein"].node_ids = protein_ids

    data["protein", "interacts_with", "protein"].edge_index = ppi_edge_index
    data["drug", "targets", "protein"].edge_index = target_edge_index

    # Integrity check: every edge index must point at a valid node in its endpoint type.
    for edge_type in data.edge_types:
        src_type, _, dst_type = edge_type
        edge_index = data[edge_type].edge_index
        num_src = data[src_type].x.shape[0]
        num_dst = data[dst_type].x.shape[0]
        if edge_index.numel() > 0:
            assert edge_index[0].max().item() < num_src, f"{edge_type} source index out of range"
            assert edge_index[1].max().item() < num_dst, f"{edge_type} destination index out of range"

    print("\n=== HeteroData summary ===")
    for node_type in data.node_types:
        print(f"  node '{node_type}': num_nodes={data[node_type].x.shape[0]}, feature_dim={data[node_type].x.shape[1]}")
    for edge_type in data.edge_types:
        print(f"  edge {edge_type}: num_edges={data[edge_type].edge_index.shape[1]}")
    print(f"  total tensor memory: {tensor_memory_mb(data):.2f} MB")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    torch.save(data, args.output)
    print(f"\nSaved HeteroData to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
