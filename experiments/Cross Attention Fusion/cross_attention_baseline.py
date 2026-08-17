"""
cross_attention_baseline.py — Cross-attention fusion + MLP head baseline for
continuous IC50 prediction.

Deep-learning bridge between the flat-fusion baselines (RF/MLP on
concatenated omics) and the eventual GAT model: GE / Mut_CNV / Proteomics
PCA embeddings are fused per cell line with `MultiOmicsCrossAttentionFusion`
(learned cross-attention, not concatenation), concatenated with a one-hot
drug vector, and passed through a small MLP head. The fusion module and the
head are trained jointly end-to-end on the regression loss.

Run from the repository root:
    python "experiments/Cross Attention Fusion/cross_attention_baseline.py"
"""

from __future__ import annotations

import copy
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import f1_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit

# Repo root (three levels up from this file: .../<repo>/experiments/Cross Attention Fusion/this.py)
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.models.cross_attention_fusion import MultiOmicsCrossAttentionFusion  # noqa: E402
from src.data.omics_preprocessing import GE_KEY, MUT_CNV_KEY, PROTEOMICS_KEY  # noqa: E402

# --- config ---------------------------------------------------------------
DATA_DIR = Path("data/processed")
ALIGNED_DIR = Path("data/processed/aligned")

# File -> modality-key mapping per src/data/00_run_preprocessing.py:31-36
OMICS_FILES = {
    GE_KEY: DATA_DIR / "transcriptomics_pca.csv",
    MUT_CNV_KEY: DATA_DIR / "genomics_pca.csv",
    PROTEOMICS_KEY: DATA_DIR / "proteomics_pca.csv",
}
TARGET_CSV = ALIGNED_DIR / "gdsc2_response_master.csv"

COL_CELL_LINE = "sanger_model_id"
COL_DRUG = "drug_id"
COL_TARGET = "ln_ic50"

TEST_FRAC = 0.15
VAL_FRAC = 0.15
RANDOM_STATE = 42
DTYPE = np.float32

# --- model / training hyperparameters --------------------------------------
TORCH_SEED = 42
D_MODEL = 128       # per-modality PCA dim (matches OmicsPreprocessingPipeline d_target)
NUM_HEADS = 4
FUSION_OUT_DIM = 256
FUSION_DROPOUT = 0.2
HEAD_HIDDEN_DIMS = [256, 128]
HEAD_DROPOUT = 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-5
BATCH_SIZE = 256
MAX_EPOCHS = 200
PATIENCE = 15
LR_PATIENCE = 5
# --------------------------------------------------------------------------


def peak_rss_gb() -> float:
    """Peak resident set size. Returns 0.0 on Windows to prevent crashes."""
    if platform.system() == "Windows":
        return 0.0

    import resource
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw / 1024**3 if platform.system() == "Darwin" else raw / 1024**2


def load_omics_modalities() -> tuple[dict[str, np.ndarray], pd.Index]:
    """Load the three per-modality PCA CSVs, verify shared cell-line index."""
    arrays: dict[str, np.ndarray] = {}
    index: pd.Index | None = None
    for key, path in OMICS_FILES.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing omics file at {path}.")
        block = pd.read_csv(path, index_col=0)
        if index is None:
            index = block.index
        elif not block.index.equals(index):
            raise ValueError(
                f"Row misalignment in '{key}' ({path.name}): expected {len(index)} "
                f"cell lines matching the first modality, got {len(block.index)}."
            )
        arrays[key] = block.to_numpy(dtype=DTYPE, copy=False)
        print(f"[features] loaded {key:<9} {path.name:<24} {arrays[key].shape}")
    return arrays, index


def load_targets(valid_ids: pd.Index) -> pd.DataFrame:
    if not TARGET_CSV.exists():
        raise FileNotFoundError(
            f"Missing target file at {TARGET_CSV}. "
            "Ensure the DE pipeline (ingest_and_align.py) has been run."
        )

    y = pd.read_csv(TARGET_CSV, usecols=[COL_CELL_LINE, COL_DRUG, COL_TARGET])
    n_raw = len(y)

    y[COL_CELL_LINE] = y[COL_CELL_LINE].astype(str)
    y[COL_DRUG] = y[COL_DRUG].astype(str)
    y[COL_TARGET] = pd.to_numeric(y[COL_TARGET], errors="coerce")

    y = y[np.isfinite(y[COL_TARGET])]
    n_finite = len(y)

    y = y[y[COL_CELL_LINE].isin(set(valid_ids))]
    n_matched = len(y)

    dup = y.duplicated([COL_CELL_LINE, COL_DRUG]).sum()
    if dup:
        print(f"[targets] {dup} duplicate (cell_line, drug) rows -> averaging")
        y = y.groupby([COL_CELL_LINE, COL_DRUG], as_index=False)[COL_TARGET].mean()

    print(f"[targets] {n_raw} rows -> {n_finite} finite -> {n_matched} with omics "
          f"-> {len(y)} unique pairs")
    if y.empty:
        raise ValueError("No target rows survived filtering; check ID formatting.")
    return y.reset_index(drop=True)


def build_pair_tensors(
    omics: dict[str, np.ndarray], cell_ids: pd.Index, y: pd.DataFrame
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, np.ndarray, int]:
    """Gather per-modality cell-line rows per pair, plus a one-hot drug block."""
    row_of = pd.Series(np.arange(len(cell_ids)), index=cell_ids)
    rows = row_of.loc[y[COL_CELL_LINE]].to_numpy()

    drug_codes, drug_levels = pd.factorize(y[COL_DRUG], sort=True)
    n_pairs, n_drugs = len(y), len(drug_levels)

    gathered = {key: arr[rows] for key, arr in omics.items()}

    drug_onehot = np.zeros((n_pairs, n_drugs), dtype=DTYPE)
    drug_onehot[np.arange(n_pairs), drug_codes] = 1.0

    target = y[COL_TARGET].to_numpy(dtype=DTYPE)
    groups = y[COL_CELL_LINE].to_numpy()

    omics_bytes = sum(a.nbytes for a in gathered.values())
    print(f"[design]  {n_pairs} pairs  |  omics: {' + '.join(f'{k}({v.shape[1]})' for k, v in gathered.items())}"
          f"  |  drug one-hot: {n_drugs}  =  {(omics_bytes + drug_onehot.nbytes) / 1024**2:.1f} MB")
    return gathered, drug_onehot, target, groups, n_drugs


def grouped_split(groups: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """70/15/15 by cell line."""
    holdout = VAL_FRAC + TEST_FRAC
    gss1 = GroupShuffleSplit(n_splits=1, test_size=holdout, random_state=RANDOM_STATE)
    train_idx, rest_idx = next(gss1.split(np.zeros(len(groups)), groups=groups))

    gss2 = GroupShuffleSplit(
        n_splits=1, test_size=TEST_FRAC / holdout, random_state=RANDOM_STATE
    )
    rel_val, rel_test = next(
        gss2.split(np.zeros(len(rest_idx)), groups=groups[rest_idx])
    )
    val_idx, test_idx = rest_idx[rel_val], rest_idx[rel_test]

    for name, idx in [("train", train_idx), ("val", val_idx), ("test", test_idx)]:
        print(f"[split]   {name:<5} {len(idx):>7} pairs  "
              f"{len(np.unique(groups[idx])):>5} cell lines  "
              f"({len(idx) / len(groups):.1%})")

    overlap = set(groups[train_idx]) & (set(groups[val_idx]) | set(groups[test_idx]))
    assert not overlap, f"cell line leaked across splits: {sorted(overlap)[:5]}"
    return train_idx, val_idx, test_idx


def evaluate(
    y_true: np.ndarray, y_pred: np.ndarray, threshold: float
) -> tuple[float, float, float, float]:
    """RMSE/PCC on the continuous target, plus AUC/F1 on a threshold-binarized target."""
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    pcc = float(np.corrcoef(y_true, y_pred)[0, 1]) if np.std(y_pred) > 0 else float("nan")

    y_true_bin = (y_true >= threshold).astype(int)
    y_pred_bin = (y_pred >= threshold).astype(int)
    auc = (
        float(roc_auc_score(y_true_bin, y_pred))
        if len(np.unique(y_true_bin)) > 1
        else float("nan")
    )
    f1 = float(f1_score(y_true_bin, y_pred_bin, zero_division=0))
    return rmse, pcc, auc, f1


def compute_shared_threshold() -> float:
    """Median ln_ic50 over all target rows, used as the AUC/F1 binarization cutoff."""
    y_all = pd.read_csv(TARGET_CSV, usecols=[COL_TARGET])
    y_all[COL_TARGET] = pd.to_numeric(y_all[COL_TARGET], errors="coerce")
    threshold = float(y_all[COL_TARGET].median(skipna=True))
    print(f"[threshold] median ln_ic50 across all targets = {threshold:.4f}")
    return threshold


class CrossAttentionRegressor(nn.Module):
    """MultiOmicsCrossAttentionFusion -> concat drug one-hot -> MLP head -> scalar."""

    def __init__(self, n_drugs: int, hidden_dims: list[int], dropout: float):
        super().__init__()
        self.fusion = MultiOmicsCrossAttentionFusion(
            d_model=D_MODEL, num_heads=NUM_HEADS, out_dim=FUSION_OUT_DIM, dropout=FUSION_DROPOUT
        )

        layers: list[nn.Module] = []
        prev = FUSION_OUT_DIM + n_drugs
        for h in hidden_dims:
            layers += [nn.Linear(prev, h), nn.BatchNorm1d(h), nn.ReLU(), nn.Dropout(dropout)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.head = nn.Sequential(*layers)

    def forward(self, omics: dict[str, torch.Tensor], drug_onehot: torch.Tensor) -> torch.Tensor:
        fused = self.fusion(omics)
        x = torch.cat([fused, drug_onehot], dim=-1)
        return self.head(x).squeeze(-1)


def make_loader(
    omics: dict[str, np.ndarray], drug_onehot: np.ndarray, y: np.ndarray,
    keys: list[str], batch_size: int, shuffle: bool,
) -> DataLoader:
    tensors = [torch.from_numpy(omics[k]) for k in keys] + [torch.from_numpy(drug_onehot), torch.from_numpy(y)]
    ds = TensorDataset(*tensors)
    # drop_last avoids a size-1 trailing batch, which BatchNorm1d rejects in train mode
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle, drop_last=shuffle)


@torch.no_grad()
def predict(
    model: nn.Module, omics: dict[str, np.ndarray], drug_onehot: np.ndarray,
    keys: list[str], device: torch.device, batch_size: int,
) -> np.ndarray:
    model.eval()
    preds = []
    n = len(drug_onehot)
    for start in range(0, n, batch_size):
        end = start + batch_size
        omics_batch = {k: torch.from_numpy(omics[k][start:end]).to(device) for k in keys}
        drug_batch = torch.from_numpy(drug_onehot[start:end]).to(device)
        preds.append(model(omics_batch, drug_batch).cpu().numpy())
    return np.concatenate(preds)


def train(
    model: nn.Module,
    train_loader: DataLoader,
    omics_val: dict[str, np.ndarray],
    drug_val: np.ndarray,
    y_val: np.ndarray,
    keys: list[str],
    device: torch.device,
) -> tuple[nn.Module, int]:
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=LR_PATIENCE
    )

    best_state = copy.deepcopy(model.state_dict())
    best_val_rmse = float("inf")
    best_epoch = 0
    epochs_since_improve = 0

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        train_loss = 0.0
        n_seen = 0
        for *omics_batches, drug_batch, yb in train_loader:
            omics_dict = {k: t.to(device) for k, t in zip(keys, omics_batches)}
            drug_batch, yb = drug_batch.to(device), yb.to(device)
            optimizer.zero_grad()
            pred = model(omics_dict, drug_batch)
            loss = criterion(pred, yb)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * len(yb)
            n_seen += len(yb)
        train_loss /= n_seen

        val_pred = predict(model, omics_val, drug_val, keys, device, BATCH_SIZE)
        val_rmse, val_pcc = float(np.sqrt(np.mean((y_val - val_pred) ** 2))), float(
            np.corrcoef(y_val, val_pred)[0, 1] if np.std(val_pred) > 0 else float("nan")
        )
        scheduler.step(val_rmse)

        improved = val_rmse < best_val_rmse - 1e-4
        if improved:
            best_val_rmse = val_rmse
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            epochs_since_improve = 0
        else:
            epochs_since_improve += 1

        if epoch == 1 or epoch % 5 == 0 or improved:
            marker = " *" if improved else ""
            print(f"[epoch {epoch:>3}] train_loss={train_loss:.4f}  "
                  f"val_rmse={val_rmse:.4f}  val_pcc={val_pcc:.4f}{marker}")

        if epochs_since_improve >= PATIENCE:
            print(f"[early stop] no val improvement for {PATIENCE} epochs "
                  f"(best epoch {best_epoch}, val_rmse={best_val_rmse:.4f})")
            break

    model.load_state_dict(best_state)
    return model, best_epoch


def main() -> None:
    t0 = time.perf_counter()

    torch.manual_seed(TORCH_SEED)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    threshold = compute_shared_threshold()

    omics, cell_ids = load_omics_modalities()
    keys = list(omics.keys())
    y_df = load_targets(cell_ids)
    gathered, drug_onehot, y, groups, n_drugs = build_pair_tensors(omics, cell_ids, y_df)
    del omics

    train_idx, val_idx, test_idx = grouped_split(groups)
    t_prep = time.perf_counter() - t0

    omics_train = {k: v[train_idx] for k, v in gathered.items()}
    omics_val = {k: v[val_idx] for k, v in gathered.items()}
    omics_test = {k: v[test_idx] for k, v in gathered.items()}

    train_loader = make_loader(omics_train, drug_onehot[train_idx], y[train_idx], keys, BATCH_SIZE, shuffle=True)

    model = CrossAttentionRegressor(n_drugs=n_drugs, hidden_dims=HEAD_HIDDEN_DIMS, dropout=HEAD_DROPOUT).to(device)
    n_params = sum(p.numel() for p in model.parameters())

    t1 = time.perf_counter()
    model, best_epoch = train(model, train_loader, omics_val, drug_onehot[val_idx], y[val_idx], keys, device)
    t_fit = time.perf_counter() - t1

    t2 = time.perf_counter()
    val_pred = predict(model, omics_val, drug_onehot[val_idx], keys, device, BATCH_SIZE)
    test_pred = predict(model, omics_test, drug_onehot[test_idx], keys, device, BATCH_SIZE)
    val_rmse, val_pcc, val_auc, val_f1 = evaluate(y[val_idx], val_pred, threshold)
    test_rmse, test_pcc, test_auc, test_f1 = evaluate(y[test_idx], test_pred, threshold)
    t_pred = time.perf_counter() - t2

    base_rmse = float(np.sqrt(np.mean((y[test_idx] - y[train_idx].mean()) ** 2)))

    print("\n" + "=" * 68)
    print("CROSS-ATTENTION FUSION + MLP HEAD — IC50 REGRESSION")
    print("=" * 68)
    print(f"{'':<12}{'RMSE':>10}{'PCC':>10}{'AUC':>10}{'F1':>10}")
    print(f"{'Validation':<12}{val_rmse:>10.4f}{val_pcc:>10.4f}{val_auc:>10.4f}{val_f1:>10.4f}")
    print(f"{'Test':<12}{test_rmse:>10.4f}{test_pcc:>10.4f}{test_auc:>10.4f}{test_f1:>10.4f}")
    print(f"{'Mean-only':<12}{base_rmse:>10.4f}{'--':>10}{'--':>10}{'--':>10}")
    print("-" * 68)
    print(f"fusion: d_model={D_MODEL} num_heads={NUM_HEADS} out_dim={FUSION_OUT_DIM}  "
          f"head_hidden={HEAD_HIDDEN_DIMS}  lr={LR}  batch_size={BATCH_SIZE}")
    print(f"params={n_params:,}  best_epoch={best_epoch}  device={device}")
    print(f"Prep {t_prep:.1f}s | Fit {t_fit:.1f}s | Predict {t_pred:.1f}s | "
          f"Total {time.perf_counter() - t0:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")
    print("=" * 68)


if __name__ == "__main__":
    main()
