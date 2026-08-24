"""cross_attention_matrix.py — Cross-attention fusion across omics subsets x drug reps.

Phase 2 of the experiment matrix (docs/plan/experiment_matrix_plan.md).
Unlike `rf_matrix.py`/`mlp_matrix.py`, which flatten every modality into one
concatenated vector, this uses `MultiOmicsCrossAttentionFusion` to learn
directional attention between modality pairs before the prediction head --
so the RF/MLP-vs-cross-attention delta at a fixed omics subset isolates the
value of *learned fusion* specifically.

Only the 4 omics subsets with >=2 modalities are run (8 cells): a single
modality has nothing to cross-attend against, and that case is already
covered by the MLP matrix.

Run from the repository root:
    python "experiments/Full Matrix/cross_attention_matrix.py"
"""

from __future__ import annotations

import copy
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.experiment_utils import (  # noqa: E402
    DRUG_ARMS,
    GE_KEY,
    MUT_CNV_KEY,
    PROTEOMICS_KEY,
    build_pair_tensors,
    compute_shared_threshold,
    evaluate,
    grouped_split,
    load_omics_subset,
    load_targets,
    mean_only_floor,
    peak_rss_gb,
    print_metric_block,
)
from src.models.cross_attention_fusion import (  # noqa: E402
    FingerprintEncoder,
    MultiOmicsCrossAttentionFusion,
)

ALL_MODALITIES = [GE_KEY, MUT_CNV_KEY, PROTEOMICS_KEY]

TORCH_SEED = 42
D_MODEL = 128
NUM_HEADS = 4
FUSION_OUT_DIM = 256
FUSION_DROPOUT = 0.2
FP_ENCODED_DIM = 128
HEAD_HIDDEN_DIMS = [256, 128]
HEAD_DROPOUT = 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-5
BATCH_SIZE = 256
MAX_EPOCHS = 200
PATIENCE = 15
LR_PATIENCE = 5

RESULTS_CSV = Path("experiments/Full Matrix/cross_attention_matrix_results.csv")


def multi_modality_subsets() -> list[list[str]]:
    """The 4 subsets with >=2 modalities (cross-attention needs a pair to attend across)."""
    subsets: list[list[str]] = []
    for size in (2, 3):
        for combo in combinations(ALL_MODALITIES, size):
            subsets.append(list(combo))
    return subsets


class CrossAttentionRegressor(nn.Module):
    """Fusion over an omics subset -> concat drug block -> MLP head -> scalar.

    For fingerprint mode the 2048-bit block is encoded down to FP_ENCODED_DIM
    first; for one-hot mode the (already modest) drug vector is concatenated
    directly, matching the original cross_attention_baseline.py behavior.
    """

    def __init__(self, modalities: list[str], n_drug_features: int, drug_mode: str):
        super().__init__()
        self.fusion = MultiOmicsCrossAttentionFusion(
            d_model=D_MODEL,
            num_heads=NUM_HEADS,
            out_dim=FUSION_OUT_DIM,
            dropout=FUSION_DROPOUT,
            modalities=modalities,
        )
        if drug_mode == "fingerprint":
            self.drug_encoder = FingerprintEncoder(
                fp_size=n_drug_features, hidden_dim=FP_ENCODED_DIM,
                out_dim=FP_ENCODED_DIM, dropout=FUSION_DROPOUT
            )
            drug_dim = FP_ENCODED_DIM
        else:
            self.drug_encoder = None
            drug_dim = n_drug_features

        layers: list[nn.Module] = []
        prev = FUSION_OUT_DIM + drug_dim
        for h in HEAD_HIDDEN_DIMS:
            layers += [nn.Linear(prev, h), nn.BatchNorm1d(h), nn.ReLU(), nn.Dropout(HEAD_DROPOUT)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.head = nn.Sequential(*layers)

    def forward(self, omics: dict[str, torch.Tensor], drug: torch.Tensor) -> torch.Tensor:
        fused = self.fusion(omics)
        drug_repr = self.drug_encoder(drug) if self.drug_encoder is not None else drug
        return self.head(torch.cat([fused, drug_repr], dim=-1)).squeeze(-1)


@torch.no_grad()
def predict(model, omics: dict[str, np.ndarray], drug: np.ndarray, keys: list[str], device) -> np.ndarray:
    model.eval()
    preds = []
    for start in range(0, len(drug), BATCH_SIZE):
        end = start + BATCH_SIZE
        ob = {k: torch.from_numpy(omics[k][start:end]).to(device) for k in keys}
        db = torch.from_numpy(drug[start:end]).to(device)
        preds.append(model(ob, db).cpu().numpy())
    return np.concatenate(preds)


def train(model, loader, omics_val, drug_val, y_val, keys, device) -> tuple[nn.Module, int]:
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=LR_PATIENCE
    )
    best_state = copy.deepcopy(model.state_dict())
    best_rmse, best_epoch, stale = float("inf"), 0, 0

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        for *omics_batches, drug_batch, yb in loader:
            od = {k: t.to(device) for k, t in zip(keys, omics_batches)}
            drug_batch, yb = drug_batch.to(device), yb.to(device)
            optimizer.zero_grad()
            loss = criterion(model(od, drug_batch), yb)
            loss.backward()
            optimizer.step()

        val_rmse = float(np.sqrt(np.mean((y_val - predict(model, omics_val, drug_val, keys, device)) ** 2)))
        scheduler.step(val_rmse)

        if val_rmse < best_rmse - 1e-4:
            best_rmse, best_epoch, stale = val_rmse, epoch, 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            stale += 1

        if epoch % 10 == 0:
            print(f"  [epoch {epoch:>3}] val_rmse={val_rmse:.4f}  best={best_rmse:.4f} @ {best_epoch}")

        if stale >= PATIENCE:
            print(f"  [early stop] epoch {epoch}, best epoch {best_epoch}, val_rmse={best_rmse:.4f}")
            break

    model.load_state_dict(best_state)
    return model, best_epoch


def run_one(modalities: list[str], arm: str, drug_mode: str, restricted: bool,
            threshold: float, device) -> dict:
    label = "+".join(modalities)
    print("\n" + "#" * 82)
    print(f"# CrossAttn | omics={label} | drug={arm}")
    print("#" * 82)

    torch.manual_seed(TORCH_SEED)
    t0 = time.perf_counter()

    omics, cell_ids = load_omics_subset(modalities)
    y_df = load_targets(cell_ids)
    gathered, drug_block, y, groups, n_drug_features, _ = build_pair_tensors(
        omics, cell_ids, y_df, drug_mode, restricted
    )
    keys = list(gathered.keys())
    train_idx, val_idx, test_idx = grouped_split(groups)
    t_prep = time.perf_counter() - t0

    omics_tr = {k: v[train_idx] for k, v in gathered.items()}
    omics_va = {k: v[val_idx] for k, v in gathered.items()}
    omics_te = {k: v[test_idx] for k, v in gathered.items()}

    tensors = [torch.from_numpy(omics_tr[k]) for k in keys] + [
        torch.from_numpy(drug_block[train_idx]),
        torch.from_numpy(y[train_idx]),
    ]
    loader = DataLoader(TensorDataset(*tensors), batch_size=BATCH_SIZE, shuffle=True, drop_last=True)

    model = CrossAttentionRegressor(modalities, n_drug_features, drug_mode).to(device)
    n_params = sum(p.numel() for p in model.parameters())

    t1 = time.perf_counter()
    model, best_epoch = train(
        model, loader, omics_va, drug_block[val_idx], y[val_idx], keys, device
    )
    t_fit = time.perf_counter() - t1

    val = evaluate(y[val_idx], predict(model, omics_va, drug_block[val_idx], keys, device), threshold)
    test = evaluate(y[test_idx], predict(model, omics_te, drug_block[test_idx], keys, device), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])

    print_metric_block(f"CrossAttn | {label} | {arm}", val, test, floor)
    print(f"n_pairs={len(y)}  n_pairs_attn={len(model.fusion.pairs)}  params={n_params:,}  "
          f"best_epoch={best_epoch}  prep={t_prep:.1f}s  fit={t_fit:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {
        "model": "CrossAttention",
        "omics": label,
        "n_modalities": len(modalities),
        "drug_rep": arm,
        "n_pairs": len(y),
        "n_attention_pairs": len(model.fusion.pairs),
        "mean_only_rmse": floor,
        **{f"val_{k}": v for k, v in val.items()},
        **{f"test_{k}": v for k, v in test.items()},
        "params": n_params,
        "best_epoch": best_epoch,
        "fit_seconds": t_fit,
    }


def main() -> None:
    t_start = time.perf_counter()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")
    threshold = compute_shared_threshold()

    results = []
    for modalities in multi_modality_subsets():
        for arm, drug_mode, restricted in DRUG_ARMS:
            results.append(run_one(modalities, arm, drug_mode, restricted, threshold, device))

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"CROSS-ATTENTION MATRIX COMPLETE — {len(df)} runs in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    summary = df[["omics", "drug_rep", "n_pairs", "test_rmse", "test_pcc", "test_r2", "test_auc"]]
    print(summary.sort_values("test_rmse").to_string(index=False))
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
