import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))
from src.models.test.hetero_gnn import HeteroIC50GNN


def compute_metrics(y_true, y_pred):
    mse = float(np.mean((y_true - y_pred) ** 2))
    rmse = float(np.sqrt(mse))
    r_score, _ = pearsonr(y_true, y_pred)
    # rho, _ = spearmanr(y_true, y_pred)
    return mse, rmse, r_score


def train_epoch(model, dataloader, x_dict, edge_index_dict, optimizer, criterion, device):
    model.train()
    total_loss = 0.0

    for cell_idx, drug_idx, targets in dataloader:
        cell_idx = cell_idx.to(device)
        drug_idx = drug_idx.to(device)
        targets = targets.to(device)

        optimizer.zero_grad()
        preds = model(x_dict, edge_index_dict, cell_idx, drug_idx)
        loss = criterion(preds, targets)
        loss.backward()
        optimizer.step()

        total_loss += loss.item() * len(targets)

    return total_loss / len(dataloader.dataset)


@torch.no_grad()
def evaluate(model, dataloader, x_dict, edge_index_dict, device):
    model.eval()
    all_preds, all_targets = [], []

    for cell_idx, drug_idx, targets in dataloader:
        cell_idx = cell_idx.to(device)
        drug_idx = drug_idx.to(device)

        preds = model(x_dict, edge_index_dict, cell_idx, drug_idx)
        all_preds.extend(preds.cpu().numpy())
        all_targets.extend(targets.numpy())

    return compute_metrics(np.array(all_targets), np.array(all_preds))


def main():
    # 1. Setup Device & Directory
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    os.makedirs("models/checkpoints", exist_ok=True)
    print(f"Using device: {device}")

    # 2. Load Graph & Move Data Tensors to Device
    graph_path = (
        "src/graph/hetero_graph.pt"
        if os.path.exists("data/processed/hetero_graph.pt")
        else "src/graph/hetero_graph.pt"
    )
    graph = torch.load(graph_path, weights_only=False)

    x_dict = {
        k: graph[k].x.to(device, dtype=torch.float32)
        for k in graph.node_types
    }
    edge_index_dict = {
        k: graph[k].edge_index.to(device, dtype=torch.long)
        for k in graph.edge_types
    }

    src, dst = edge_index_dict[("drug", "targets", "protein")]
    edge_index_dict[("protein", "rev_targets", "drug")] = torch.stack([dst, src])

    src, dst = edge_index_dict[("cell_line", "has_mutation", "protein")]
    edge_index_dict[("protein", "rev_has_mutation", "cell_line")] = torch.stack([dst, src])
    # 3. Load Supervisory Pairs
    df = pd.read_csv("data/raw/aligned_ic50_pairs.csv")
    target_col = "ln_ic50"
    assert target_col in df.columns, f"{target_col} missing; have {df.columns.tolist()}"
    
    # print("target_col:", target_col)
    # print("columns:", df.columns.tolist())
    # print(df[target_col].describe())
    # print("skew:", df[target_col].skew())
    # import sys; sys.exit()

    train_df = df[df["split"] == "train"]
    val_df = df[df["split"] == "val"]
    test_df = df[df["split"] == "test"]

    # 4. Create PyTorch DataLoaders
    def make_loader(data_frame, shuffle=False, batch_size=1024, drop_last=False):
        dataset = TensorDataset(
            torch.tensor(data_frame["cell_idx"].values, dtype=torch.long),
            torch.tensor(data_frame["drug_idx"].values, dtype=torch.long),
            torch.tensor(data_frame[target_col].values, dtype=torch.float32),
        )
        return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, drop_last=drop_last)

    # drop_last so BatchNorm in the head never sees a 1-sample final batch.
    train_loader = make_loader(train_df, shuffle=True, batch_size=1024, drop_last=True)
    val_loader = make_loader(val_df, shuffle=False, batch_size=2048)
    test_loader = make_loader(test_df, shuffle=False, batch_size=2048)

    # 5. Initialize Model, Loss, and Optimizer
    num_proteins = x_dict["protein"].shape[0]
    model = HeteroIC50GNN(num_proteins=num_proteins, hidden_dim=128).to(device)
    # weight_decay overridable from the CLI so the 1e-4 vs 1e-5 arm can be
    # compared without editing the file; 1e-5 matches the tracked HeteroGNN.
    weight_decay = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-5
    print(f"[config] weight_decay={weight_decay:g}")
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=weight_decay)
    # Halve the LR when val RMSE plateaus, as the tracked gnn_baseline.py does.
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=5
    )
    criterion = nn.MSELoss()

    # 6. Training Loop
    max_epochs = 200
    patience = 15
    best_val_rmse = float("inf")
    best_epoch = 0
    stale = 0

    print(f"\nStarting GNN Training for up to {max_epochs} epochs (patience={patience})...\n" + "-" * 55)
    for epoch in range(1, max_epochs + 1):
        train_loss = train_epoch(
            model, train_loader, x_dict, edge_index_dict, optimizer, criterion, device
        )
        val_mse, val_rmse, val_r = evaluate(
            model, val_loader, x_dict, edge_index_dict, device
        )
        scheduler.step(val_rmse)

        print(
            f"Epoch {epoch:02d}/{max_epochs:02d} | "
            f"Train Loss: {train_loss:.4f} | "
            f"Val MSE: {val_mse:.4f} | "
            f"Val RMSE: {val_rmse:.4f} | "
            f"Val Pearson r: {val_r:.4f}"
        )

        # Save Best Checkpoint (selected on val RMSE, the metric we actually
        # care about minimizing -- not Pearson r, which kept climbing well
        # past the point the model started overfitting on RMSE).
        if val_rmse < best_val_rmse - 1e-4:
            best_val_rmse, best_epoch, stale = val_rmse, epoch, 0
            torch.save(
                model.state_dict(), "models/checkpoints/best_hetero_gnn.pt"
            )
        else:
            stale += 1
            if stale >= patience:
                print(f"[early stop] epoch {epoch}, best epoch {best_epoch}, val_rmse={best_val_rmse:.4f}")
                break

    # 7. Final Test Set Evaluation
    print("-" * 55)
    model.load_state_dict(torch.load("models/checkpoints/best_hetero_gnn.pt"))
    test_mse, test_rmse, test_r = evaluate(
        model, test_loader, x_dict, edge_index_dict, device
    )
    print(f"BEST MODEL TEST METRICS:")
    print(f"Test MSE: {test_mse:.4f} | Test RMSE: {test_rmse:.4f} | Test Pearson r: {test_r:.4f}")


if __name__ == "__main__":
    main()