"""full_architecture_matrix.py — the proposed architecture, end to end.

Runs `src/models/full_architecture.py`: cross-attention fusion over all three
omics modalities, an atom-level GCN over drug molecular graphs, message
passing over the STRING PPI graph, and an MLP head -- every box in
`docs/Architecture Simplified.jpeg` in one model.

The controlled comparison is **E13** (GNN-GCN, tri-omics, fingerprint,
+mutation edges, RMSE 1.3513 in docs/07_gnn_ablation_results.md), which uses
this same `HeteroGNN` over this same graph. Two things change:

- `cell_line` node features: raw concatenated PCA blocks (384) -> cross-attention
  fused embedding (256).
- `drug` node features: 2048-bit Morgan fingerprint -> molecular graph GCN
  embedding (128).

Everything else -- graph, edge types, conv variant, hidden dim, layer count,
optimizer, split, seed -- is pinned to `experiments/07_gnn_ablation/gnn_baseline.py`,
so the delta is attributable to those two encoder swaps.

Both conv variants are swept because the benchmark's own Table 3 and this
project's E57/E59 disagree with the proposal's GAT choice; see
docs/07_gnn_ablation_results.md.

Run from the repository root:
    python "experiments/12_full_architecture/full_architecture_matrix.py"
    python "experiments/12_full_architecture/full_architecture_matrix.py" --variant gcn
"""

from __future__ import annotations

import argparse
import copy
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset
from torch_geometric.data import HeteroData

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.drug_graphs import build_drug_graphs  # noqa: E402
from src.data.experiment_utils import (  # noqa: E402
    COL_CELL_LINE,
    COL_DRUG,
    COL_TARGET,
    GE_KEY,
    MUT_CNV_KEY,
    PROTEOMICS_KEY,
    compute_shared_threshold,
    evaluate,
    grouped_split,
    load_omics_subset,
    load_targets,
    mean_only_floor,
    peak_rss_gb,
    print_metric_block,
)
from src.models.full_architecture import (  # noqa: E402
    FullArchitectureRegressor,
    align_drug_graphs_to_graph,
    align_omics_to_graph,
)

GRAPH_PATH = Path("src/graph/hetero_graph.pt")
RESULTS_CSV = Path("experiments/12_full_architecture/full_architecture_results.csv")

ALL_MODALITIES = [GE_KEY, MUT_CNV_KEY, PROTEOMICS_KEY]
MUTATION_EDGE = ("cell_line", "has_mutation", "protein")

# Pinned to experiments/07_gnn_ablation/gnn_baseline.py so the comparison
# against E13 isolates the two encoder swaps.
TORCH_SEED = 42
HIDDEN_DIM = 128
NUM_LAYERS = 2
HEADS = 4
DROPOUT = 0.2
HEAD_HIDDEN_DIMS = (256, 128)
HEAD_DROPOUT = 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-5
BATCH_SIZE = 1024  # larger than the flat models': every step re-runs full-graph message passing
MAX_EPOCHS = 200
PATIENCE = 15
LR_PATIENCE = 5

# E13's population, for a comparability check rather than a hard requirement.
EXPECTED_PAIRS = 111_799


def load_graph_and_pairs():
    """Load the hetero graph and join GDSC IC50 labels onto its node index space."""
    if not GRAPH_PATH.exists():
        raise FileNotFoundError(
            f"Missing {GRAPH_PATH}. Run src/data/03_graph_construction.py "
            f"then src/data/04_link_cell_lines.py first."
        )

    data: HeteroData = torch.load(GRAPH_PATH, weights_only=False)
    cell_ids = list(data["cell_line"].node_ids)
    drug_ids = list(data["drug"].node_ids)
    print(f"[graph] cell_line={len(cell_ids)}  drug={len(drug_ids)}  protein={len(data['protein'].node_ids)}")
    for et in data.edge_types:
        print(f"[graph] edge {et}: {data[et].edge_index.shape[1]}")

    y_df = load_targets(pd.Index(cell_ids))

    drug_to_index = {d: i for i, d in enumerate(drug_ids)}
    cell_to_index = {c: i for i, c in enumerate(cell_ids)}
    before = len(y_df)
    y_df = y_df[y_df[COL_DRUG].isin(drug_to_index)].reset_index(drop=True)
    print(f"[pairs] {before - len(y_df)} pairs dropped (drug absent from graph) -> {len(y_df)} remain")
    if len(y_df) != EXPECTED_PAIRS:
        print(
            f"[warn] {len(y_df)} pairs != E13's {EXPECTED_PAIRS}; the comparison "
            f"against docs/results.md is no longer like-for-like."
        )

    cell_index = y_df[COL_CELL_LINE].map(cell_to_index).to_numpy(dtype=np.int64)
    drug_index = y_df[COL_DRUG].map(drug_to_index).to_numpy(dtype=np.int64)
    target = y_df[COL_TARGET].to_numpy(dtype=np.float32)
    groups = y_df[COL_CELL_LINE].to_numpy()

    return data, cell_index, drug_index, target, groups


@torch.no_grad()
def predict(model, edge_dict, cell_index: np.ndarray, drug_index: np.ndarray, device) -> np.ndarray:
    model.eval()
    preds = []
    for start in range(0, len(cell_index), BATCH_SIZE):
        end = start + BATCH_SIZE
        ci = torch.from_numpy(cell_index[start:end]).to(device)
        di = torch.from_numpy(drug_index[start:end]).to(device)
        preds.append(model(edge_dict, ci, di).cpu().numpy())
    return np.concatenate(preds)


def train(model, edge_dict, loader, cell_va, drug_va, y_va, device) -> tuple[nn.Module, int]:
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=LR_PATIENCE
    )
    best_state = copy.deepcopy(model.state_dict())
    best_rmse, best_epoch, stale = float("inf"), 0, 0

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        for ci, di, yb in loader:
            ci, di, yb = ci.to(device), di.to(device), yb.to(device)
            optimizer.zero_grad()
            loss = criterion(model(edge_dict, ci, di), yb)
            loss.backward()
            optimizer.step()

        val_rmse = float(np.sqrt(np.mean((y_va - predict(model, edge_dict, cell_va, drug_va, device)) ** 2)))
        scheduler.step(val_rmse)

        if val_rmse < best_rmse - 1e-4:
            best_rmse, best_epoch, stale = val_rmse, epoch, 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            stale += 1

        if epoch % 5 == 0:
            print(f"  [epoch {epoch:>3}] val_rmse={val_rmse:.4f}  best={best_rmse:.4f} @ {best_epoch}")

        if stale >= PATIENCE:
            print(f"  [early stop] epoch {epoch}, best epoch {best_epoch}, val_rmse={best_rmse:.4f}")
            break

    model.load_state_dict(best_state)
    return model, best_epoch


def run_one(data, omics, batched, cell_index, drug_index, y, groups,
            modalities: list[str], variant: str, threshold: float, device) -> dict:
    label = f"{variant.upper()} | omics={'+'.join(modalities)}"
    print("\n" + "#" * 82)
    print(f"# FULL ARCHITECTURE | {label}")
    print("#" * 82)

    torch.manual_seed(TORCH_SEED)
    t0 = time.perf_counter()

    edge_dict = {et: data[et].edge_index.to(device) for et in data.edge_types}
    print(f"[edges] using {len(edge_dict)} edge types: {list(edge_dict.keys())}")

    train_idx, val_idx, test_idx = grouped_split(groups)
    t_prep = time.perf_counter() - t0

    ds = TensorDataset(
        torch.from_numpy(cell_index[train_idx]),
        torch.from_numpy(drug_index[train_idx]),
        torch.from_numpy(y[train_idx]),
    )
    loader = DataLoader(ds, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)

    model = FullArchitectureRegressor(
        omics={k: v.to(device) for k, v in omics.items()},
        batched=batched,
        metadata=(list(data.node_types), list(edge_dict.keys())),
        num_proteins=data["protein"].x.shape[0],
        variant=variant,
        hidden_dim=HIDDEN_DIM,
        num_layers=NUM_LAYERS,
        heads=HEADS,
        dropout=DROPOUT,
        head_hidden_dims=HEAD_HIDDEN_DIMS,
        head_dropout=HEAD_DROPOUT,
    ).to(device)

    # SAGEConv/GATConv are built with (-1, -1) input dims, so their weights are
    # lazily created on the first forward pass. Run one dummy batch to
    # materialize them before counting parameters or handing them to Adam.
    model.eval()  # BatchNorm in the head rejects a 1-sample batch in train mode
    with torch.no_grad():
        model(
            edge_dict,
            torch.from_numpy(cell_index[train_idx[:2]]).to(device),
            torch.from_numpy(drug_index[train_idx[:2]]).to(device),
        )
    n_params = sum(p.numel() for p in model.parameters())

    t1 = time.perf_counter()
    model, best_epoch = train(
        model, edge_dict, loader, cell_index[val_idx], drug_index[val_idx], y[val_idx], device
    )
    t_fit = time.perf_counter() - t1

    val = evaluate(y[val_idx], predict(model, edge_dict, cell_index[val_idx], drug_index[val_idx], device), threshold)
    test = evaluate(y[test_idx], predict(model, edge_dict, cell_index[test_idx], drug_index[test_idx], device), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])

    print_metric_block(f"FULL | {label}", val, test, floor)
    print(f"n_pairs={len(y)}  n_drugs={batched.n_graphs}  n_pairs_attn={len(model.fusion.pairs)}  "
          f"params={n_params:,}  best_epoch={best_epoch}  prep={t_prep:.1f}s  fit={t_fit:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {
        "model": f"FullArchitecture-{variant.upper()}",
        "omics": "+".join(modalities),
        "n_modalities": len(modalities),
        "drug_rep": "molecular_graph",
        "fusion": "cross_attention",
        "graph": "PPI +mutation edges",
        "n_pairs": len(y),
        "n_attention_pairs": len(model.fusion.pairs),
        "n_edge_types": len(edge_dict),
        "mean_only_rmse": floor,
        **{f"val_{k}": v for k, v in val.items()},
        **{f"test_{k}": v for k, v in test.items()},
        "params": n_params,
        "best_epoch": best_epoch,
        "fit_seconds": t_fit,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """`argv` is explicit so a notebook can call `main([...])` without argparse
    picking up the kernel's own `-f kernel.json` argument."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--variant", choices=("gcn", "gat"), default=None,
        help="Run a single conv variant. Defaults to sweeping both.",
    )
    parser.add_argument(
        "--omics", nargs="+", choices=ALL_MODALITIES, default=None,
        help="Omics subset (>=2 modalities). Defaults to all three, as proposed.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    modalities = list(dict.fromkeys(args.omics)) if args.omics else list(ALL_MODALITIES)
    if len(modalities) < 2:
        raise SystemExit("--omics needs >=2 modalities: cross-attention attends between modalities.")
    variants = [args.variant] if args.variant else ["gcn", "gat"]

    t_start = time.perf_counter()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    threshold = compute_shared_threshold()

    data, cell_index, drug_index, y, groups = load_graph_and_pairs()
    if MUTATION_EDGE not in data.edge_types:
        raise SystemExit(
            f"{MUTATION_EDGE} absent from the graph — cell_line nodes would be "
            f"disconnected and message passing could not reach them. Run "
            f"src/data/04_link_cell_lines.py first."
        )

    # Aligned once to the graph's node order, then shared across variants.
    omics_raw, cell_ids = load_omics_subset(modalities)
    omics = align_omics_to_graph(omics_raw, cell_ids, data["cell_line"].node_ids)
    batched = align_drug_graphs_to_graph(build_drug_graphs(), data["drug"].node_ids)

    results = [
        run_one(data, omics, batched, cell_index, drug_index, y, groups,
                modalities, variant, threshold, device)
        for variant in variants
    ]

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"FULL ARCHITECTURE COMPLETE — {len(df)} runs in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    summary = df[["model", "omics", "n_pairs", "test_rmse", "test_pcc", "test_r2", "test_auc"]]
    print(summary.sort_values("test_rmse").to_string(index=False))
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
