"""Shared setup for the interpretability checks: graph, trained model, and untrained controls.

The untrained controls are `HeteroIC50GNN`s with the checkpoint's exact
architecture and the exact same graph, but freshly initialised weights that
never see an IC50 label. Anything their attributions show is therefore a
product of the graph's wiring and the architecture alone, which is the null
the trained model's attributions have to beat.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import torch
from torch import nn

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.interpretation.genes import load_symbol_map  # noqa: E402
from src.models.train.hetero_gnn import HeteroIC50GNN  # noqa: E402

GRAPH_PATH = REPO_ROOT / "src" / "graph" / "hetero_graph.pt"
CHECKPOINT_PATH = REPO_ROOT / "models" / "checkpoints" / "best_hetero_gnn.pt"
TARGET_CSV = REPO_ROOT / "data" / "processed" / "aligned" / "gdsc2_response_master.csv"
MUTATIONS_PATH = REPO_ROOT / "data" / "processed" / "aligned" / "mutations_aligned.csv"
OUT_DIR = REPO_ROOT / "experiments" / "15_interpretability_validation"

HIDDEN_DIM = 128
N_RANDOM_MODELS = 50
# Well clear of the training seed (42) so no control shares the checkpoint's init.
RANDOM_SEED_START = 1000

BRAF_SYMBOL = "BRAF"
BRAF_INHIBITORS = {"1036": "PLX-4720", "1061": "SB590885", "1373": "Dabrafenib"}


def load_graph_and_maps(device: torch.device) -> dict:
    """Same construction as backend/model_backends/hetero_gnn.py::_load_graph."""
    graph = torch.load(GRAPH_PATH, weights_only=False)
    x_dict = {nt: graph[nt].x.to(device, dtype=torch.float32) for nt in graph.node_types}
    edge_index_dict = {et: graph[et].edge_index.to(device, dtype=torch.long) for et in graph.edge_types}

    src, dst = edge_index_dict[("drug", "targets", "protein")]
    edge_index_dict[("protein", "rev_targets", "drug")] = torch.stack([dst, src])
    src, dst = edge_index_dict[("cell_line", "has_mutation", "protein")]
    edge_index_dict[("protein", "rev_has_mutation", "cell_line")] = torch.stack([dst, src])

    return {
        "x_dict": x_dict,
        "edge_index_dict": edge_index_dict,
        "cell_to_idx": {str(c): i for i, c in enumerate(graph["cell_line"].node_ids)},
        "drug_to_idx": {str(d): i for i, d in enumerate(graph["drug"].node_ids)},
        "num_proteins": int(x_dict["protein"].shape[0]),
        "protein_node_ids": [str(p) for p in graph["protein"].node_ids],
    }


def load_trained_model(graph_maps: dict, device: torch.device) -> nn.Module:
    model = HeteroIC50GNN(num_proteins=graph_maps["num_proteins"], hidden_dim=HIDDEN_DIM).to(device)
    model.load_state_dict(torch.load(CHECKPOINT_PATH, map_location=device))
    return model.eval()


def build_random_models(
    graph_maps: dict, device: torch.device, n: int = N_RANDOM_MODELS, seed_start: int = RANDOM_SEED_START
) -> list[nn.Module]:
    models = []
    for seed in range(seed_start, seed_start + n):
        torch.manual_seed(seed)
        model = HeteroIC50GNN(num_proteins=graph_maps["num_proteins"], hidden_dim=HIDDEN_DIM).to(device)
        models.append(model.eval())
    return models


def score_pairs(model: nn.Module, graph_maps: dict, cell_idx: list[int], drug_idx: list[int]) -> torch.Tensor:
    """Gradient x embedding per protein for each pair -> [n_pairs, n_proteins].

    Identical to calling `attribute_pair` once per pair: in eval mode every
    output row depends only on its own pair, so one forward pass followed by a
    per-row gradient gives the same numbers without re-running message passing
    for every pair.
    """
    model.eval()
    embedding = model.protein_embedding.weight
    embedding.requires_grad_(True)
    device = embedding.device

    cell_index = torch.tensor(cell_idx, dtype=torch.long, device=device)
    drug_index = torch.tensor(drug_idx, dtype=torch.long, device=device)
    predictions = model(graph_maps["x_dict"], graph_maps["edge_index_dict"], cell_index, drug_index)

    rows = []
    for i in range(len(cell_idx)):
        (gradient,) = torch.autograd.grad(predictions[i], embedding, retain_graph=i < len(cell_idx) - 1)
        rows.append((gradient * embedding).sum(dim=-1).detach())
    return torch.stack(rows)


def protein_index_for_symbol(graph_maps: dict, symbol: str) -> int:
    symbol_map = load_symbol_map(REPO_ROOT / "data" / "processed" / "protein_symbol_map.csv")
    matches = [i for i, pid in enumerate(graph_maps["protein_node_ids"]) if symbol_map.get(pid) == symbol]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one protein node for {symbol}, found {len(matches)}")
    return matches[0]


def driver_mutant_cell_lines(symbol: str) -> set[str]:
    """Cell lines with a cancer-driver mutation in `symbol` -- the same filter 04_link_cell_lines.py uses."""
    cells: set[str] = set()
    reader = pd.read_csv(
        MUTATIONS_PATH, usecols=["standard_model_id", "gene_symbol", "cancer_driver"], dtype=str, chunksize=1_000_000
    )
    for chunk in reader:
        hit = chunk[(chunk["gene_symbol"] == symbol) & (chunk["cancer_driver"].str.lower() == "t")]
        cells.update(hit["standard_model_id"])
    return cells


def braf_case_study_pairs(graph_maps: dict) -> pd.DataFrame:
    """Real GDSC2 response rows for BRAF-driver-mutant cell lines on BRAF inhibitors."""
    frame = pd.read_csv(
        TARGET_CSV, usecols=["sanger_model_id", "drug_id", "drug_name", "cell_line_name", "ln_ic50"], dtype=str
    )
    mutants = driver_mutant_cell_lines(BRAF_SYMBOL)
    frame = frame[frame["drug_id"].isin(BRAF_INHIBITORS) & frame["sanger_model_id"].isin(mutants)].copy()
    frame["cell_idx"] = frame["sanger_model_id"].map(graph_maps["cell_to_idx"])
    frame["drug_idx"] = frame["drug_id"].map(graph_maps["drug_to_idx"])
    frame = frame.dropna(subset=["cell_idx", "drug_idx"])
    frame["cell_idx"] = frame["cell_idx"].astype(int)
    frame["drug_idx"] = frame["drug_idx"].astype(int)
    frame["ln_ic50"] = frame["ln_ic50"].astype(float)
    return frame.drop_duplicates(["sanger_model_id", "drug_id"]).reset_index(drop=True)
