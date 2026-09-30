"""Retrain HeteroIC50GNN on the original or target-edge-ablated graph with 05_train.py's recipe.

Both arms use seed 42 and the same split file, so an `original` vs `ablated`
pair differs only in the 7 removed probe-drug target edges. The production
checkpoint can't serve as the control: 05_train.py never seeds, so comparing
against it would fold unknown seed noise into the ablation effect.

The probe drugs' response rows stay in training with their labels -- only
their edges are gone. This is not leave-drugs-out.

Run from the repository root:
    python experiments/16_target_edge_ablation/train_ablated.py --graph original
    python experiments/16_target_edge_ablation/train_ablated.py --graph ablated
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.models.train.hetero_gnn import HeteroIC50GNN  # noqa: E402

OUT_DIR = REPO_ROOT / "experiments" / "16_target_edge_ablation"
GRAPHS = {
    "original": REPO_ROOT / "src" / "graph" / "hetero_graph.pt",
    "ablated": OUT_DIR / "hetero_graph_ablated.pt",
}
SPLIT_CSV = REPO_ROOT / "data" / "raw" / "aligned_ic50_pairs.csv"

# 05_train.py's recipe.
HIDDEN_DIM = 128
LR = 1e-3
WEIGHT_DECAY = 1e-5
MAX_EPOCHS = 200
PATIENCE = 15
TRAIN_BATCH = 1024
EVAL_BATCH = 2048
TORCH_SEED = 42


def load_graph(path: Path, device: torch.device) -> tuple[dict, dict]:
    graph = torch.load(path, weights_only=False)
    x_dict = {nt: graph[nt].x.to(device, dtype=torch.float32) for nt in graph.node_types}
    edge_index_dict = {et: graph[et].edge_index.to(device, dtype=torch.long) for et in graph.edge_types}
    src, dst = edge_index_dict[("drug", "targets", "protein")]
    edge_index_dict[("protein", "rev_targets", "drug")] = torch.stack([dst, src])
    src, dst = edge_index_dict[("cell_line", "has_mutation", "protein")]
    edge_index_dict[("protein", "rev_has_mutation", "cell_line")] = torch.stack([dst, src])
    return x_dict, edge_index_dict


def make_loader(frame: pd.DataFrame, batch_size: int, shuffle: bool = False, drop_last: bool = False) -> DataLoader:
    dataset = TensorDataset(
        torch.tensor(frame["cell_idx"].to_numpy(), dtype=torch.long),
        torch.tensor(frame["drug_idx"].to_numpy(), dtype=torch.long),
        torch.tensor(frame["ln_ic50"].to_numpy(), dtype=torch.float32),
    )
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, drop_last=drop_last)


@torch.no_grad()
def evaluate(model: nn.Module, loader: DataLoader, x_dict: dict, edge_index_dict: dict, device) -> tuple[float, float]:
    model.eval()
    preds, targets = [], []
    for cell_idx, drug_idx, y in loader:
        preds.append(model(x_dict, edge_index_dict, cell_idx.to(device), drug_idx.to(device)).cpu().numpy())
        targets.append(y.numpy())
    pred, true = np.concatenate(preds), np.concatenate(targets)
    return float(np.sqrt(np.mean((pred - true) ** 2))), float(pearsonr(true, pred)[0])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", choices=sorted(GRAPHS), required=True)
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[setup] device={device}  graph={args.graph}  seed={TORCH_SEED}")

    x_dict, edge_index_dict = load_graph(GRAPHS[args.graph], device)
    frame = pd.read_csv(SPLIT_CSV, usecols=["cell_idx", "drug_idx", "ln_ic50", "split"])
    train_loader = make_loader(frame[frame["split"] == "train"], TRAIN_BATCH, shuffle=True, drop_last=True)
    val_loader = make_loader(frame[frame["split"] == "val"], EVAL_BATCH)
    test_loader = make_loader(frame[frame["split"] == "test"], EVAL_BATCH)

    torch.manual_seed(TORCH_SEED)
    model = HeteroIC50GNN(num_proteins=x_dict["protein"].shape[0], hidden_dim=HIDDEN_DIM).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=5)
    criterion = nn.MSELoss()

    checkpoint = OUT_DIR / "checkpoints" / f"hetero_gnn_{args.graph}_seed{TORCH_SEED}.pt"
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    history = []
    best_val, best_epoch, stale = float("inf"), 0, 0
    start = time.time()

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        running = 0.0
        for cell_idx, drug_idx, y in train_loader:
            y = y.to(device)
            optimizer.zero_grad()
            loss = criterion(model(x_dict, edge_index_dict, cell_idx.to(device), drug_idx.to(device)), y)
            loss.backward()
            optimizer.step()
            running += loss.item() * len(y)
        train_loss = running / len(train_loader.dataset)
        val_rmse, val_pcc = evaluate(model, val_loader, x_dict, edge_index_dict, device)
        scheduler.step(val_rmse)

        history.append({"epoch": epoch, "train_loss": train_loss, "val_rmse": val_rmse, "val_pcc": val_pcc,
                        "lr": optimizer.param_groups[0]["lr"], "seconds": time.time() - start})
        print(f"epoch {epoch:3d} | train {train_loss:.4f} | val rmse {val_rmse:.4f} pcc {val_pcc:.4f} "
              f"| {time.time() - start:.0f}s")

        # Same selection rule as 05_train.py.
        if val_rmse < best_val - 1e-4:
            best_val, best_epoch, stale = val_rmse, epoch, 0
            torch.save(model.state_dict(), checkpoint)
        else:
            stale += 1
            if stale >= PATIENCE:
                print(f"[early stop] epoch {epoch}, best epoch {best_epoch}, val_rmse={best_val:.4f}")
                break

    pd.DataFrame(history).to_csv(OUT_DIR / f"train_history_{args.graph}.csv", index=False)
    model.load_state_dict(torch.load(checkpoint, map_location=device))
    test_rmse, test_pcc = evaluate(model, test_loader, x_dict, edge_index_dict, device)
    print(f"\n[{args.graph}] best epoch {best_epoch} | val rmse {best_val:.4f} | test rmse {test_rmse:.4f} "
          f"pcc {test_pcc:.4f} | {time.time() - start:.0f}s")
    print(f"Saved -> {checkpoint.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
