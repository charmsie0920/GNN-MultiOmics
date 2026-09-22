"""molecular_graph_matrix.py — Cross-attention fusion x atom-level molecular graph drugs.

Phase 1 of the benchmark-alignment work (docs/results.md). The matrix so far
represents a drug either as a one-hot identity or as a 2048-bit Morgan
fingerprint; this adds the representation MoGraphDRP uses (their section
2.2.2) -- the molecule as a graph of atoms and bonds, encoded by a 3-layer GCN.

The controlled comparison is against **E10** (CrossAttention, GE+Proteomics,
fingerprint, RMSE 1.3205): identical fusion, identical omics, identical head,
identical split, with only the drug encoder swapped. Every architectural
constant below is therefore pinned to the value
`experiments/06_full_matrix/cross_attention_matrix.py` uses, and the
population is asserted identical (111,799 pairs over 498 drugs) before
anything trains. One-hot arms are deliberately absent: drug identity cannot
generalize to unseen compounds (R^2 0.024 in
docs/10_leave_drugs_out_results.md), so it is not a representation this
project builds on.

Run from the repository root:
    python "experiments/11_molecular_graph/molecular_graph_matrix.py"
    python "experiments/11_molecular_graph/molecular_graph_matrix.py" --omics GE Proteomics
"""

from __future__ import annotations

import argparse
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
from src.data.drug_graphs import (  # noqa: E402
    BatchedMolGraphs,
    assert_fingerprint_population_parity,
    build_drug_graphs,
    build_graph_pair_tensors,
)
from src.data.experiment_utils import (  # noqa: E402
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
from src.models.cross_attention_fusion import MultiOmicsCrossAttentionFusion  # noqa: E402
from src.models.drug_gcn import MolecularGraphEncoder  # noqa: E402

ALL_MODALITIES = [GE_KEY, MUT_CNV_KEY, PROTEOMICS_KEY]

# Pinned to cross_attention_matrix.py so the fingerprint-vs-graph delta is
# attributable to the drug encoder alone. Do not tune these independently
# without re-running the fingerprint arm under the same values.
TORCH_SEED = 42
D_MODEL = 128
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

RESULTS_CSV = Path("experiments/11_molecular_graph/molecular_graph_matrix_results.csv")


def multi_modality_subsets() -> list[list[str]]:
    """The 4 subsets with >=2 modalities (cross-attention needs a pair to attend across)."""
    subsets: list[list[str]] = []
    for size in (2, 3):
        for combo in combinations(ALL_MODALITIES, size):
            subsets.append(list(combo))
    return subsets


class CrossAttentionGraphRegressor(nn.Module):
    """Fusion over an omics subset -> concat GCN drug embedding -> MLP head -> scalar.

    Structurally identical to `cross_attention_matrix.CrossAttentionRegressor`
    apart from the drug branch: instead of a dense feature row per pair, each
    pair carries an int64 *code* indexing the molecular graph embedding table.
    The code never reaches the head -- it is a lookup key, not a feature, so
    the drug is represented purely by its structure.

    The whole ~498-molecule table is encoded once per forward pass and indexed,
    rather than re-encoding a molecule per pair; see `src/models/drug_gcn.py`.
    """

    def __init__(self, modalities: list[str], batched: BatchedMolGraphs):
        super().__init__()
        self.fusion = MultiOmicsCrossAttentionFusion(
            d_model=D_MODEL,
            num_heads=NUM_HEADS,
            out_dim=FUSION_OUT_DIM,
            dropout=FUSION_DROPOUT,
            modalities=modalities,
        )
        self.drug_encoder = MolecularGraphEncoder(batched)

        layers: list[nn.Module] = []
        prev = FUSION_OUT_DIM + self.drug_encoder.out_dim
        for h in HEAD_HIDDEN_DIMS:
            layers += [nn.Linear(prev, h), nn.BatchNorm1d(h), nn.ReLU(), nn.Dropout(HEAD_DROPOUT)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.head = nn.Sequential(*layers)

    def forward(
        self,
        omics: dict[str, torch.Tensor],
        drug_codes: torch.Tensor,
        drug_table: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """`drug_table` lets inference encode the molecules once for all batches."""
        fused = self.fusion(omics)
        if drug_table is None:
            drug_table = self.drug_encoder()
        return self.head(torch.cat([fused, drug_table[drug_codes]], dim=-1)).squeeze(-1)


@torch.no_grad()
def predict(model, omics: dict[str, np.ndarray], codes: np.ndarray, keys: list[str], device) -> np.ndarray:
    model.eval()
    drug_table = model.drug_encoder()  # constant across batches at inference
    preds = []
    for start in range(0, len(codes), BATCH_SIZE):
        end = start + BATCH_SIZE
        ob = {k: torch.from_numpy(omics[k][start:end]).to(device) for k in keys}
        cb = torch.from_numpy(codes[start:end]).to(device)
        preds.append(model(ob, cb, drug_table).cpu().numpy())
    return np.concatenate(preds)


def train(model, loader, omics_val, codes_val, y_val, keys, device) -> tuple[nn.Module, int]:
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=LR_PATIENCE
    )
    best_state = copy.deepcopy(model.state_dict())
    best_rmse, best_epoch, stale = float("inf"), 0, 0

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        for *omics_batches, code_batch, yb in loader:
            od = {k: t.to(device) for k, t in zip(keys, omics_batches)}
            code_batch, yb = code_batch.to(device), yb.to(device)
            optimizer.zero_grad()
            loss = criterion(model(od, code_batch), yb)
            loss.backward()
            optimizer.step()

        val_rmse = float(np.sqrt(np.mean((y_val - predict(model, omics_val, codes_val, keys, device)) ** 2)))
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


def run_one(modalities: list[str], graphs: dict, threshold: float, device) -> dict:
    label = "+".join(modalities)
    print("\n" + "#" * 82)
    print(f"# CrossAttn | omics={label} | drug=molecular_graph")
    print("#" * 82)

    torch.manual_seed(TORCH_SEED)
    t0 = time.perf_counter()

    omics, cell_ids = load_omics_subset(modalities)
    y_df = load_targets(cell_ids)
    gathered, codes, y, groups, batched, _ = build_graph_pair_tensors(
        omics, cell_ids, y_df, graphs
    )
    keys = list(gathered.keys())
    train_idx, val_idx, test_idx = grouped_split(groups)
    t_prep = time.perf_counter() - t0

    omics_tr = {k: v[train_idx] for k, v in gathered.items()}
    omics_va = {k: v[val_idx] for k, v in gathered.items()}
    omics_te = {k: v[test_idx] for k, v in gathered.items()}

    tensors = [torch.from_numpy(omics_tr[k]) for k in keys] + [
        torch.from_numpy(codes[train_idx]),
        torch.from_numpy(y[train_idx]),
    ]
    loader = DataLoader(TensorDataset(*tensors), batch_size=BATCH_SIZE, shuffle=True, drop_last=True)

    model = CrossAttentionGraphRegressor(modalities, batched).to(device)
    n_params = sum(p.numel() for p in model.parameters())

    t1 = time.perf_counter()
    model, best_epoch = train(model, loader, omics_va, codes[val_idx], y[val_idx], keys, device)
    t_fit = time.perf_counter() - t1

    val = evaluate(y[val_idx], predict(model, omics_va, codes[val_idx], keys, device), threshold)
    test = evaluate(y[test_idx], predict(model, omics_te, codes[test_idx], keys, device), threshold)
    floor = mean_only_floor(y[train_idx], y[test_idx])

    print_metric_block(f"CrossAttn | {label} | molecular_graph", val, test, floor)
    print(f"n_pairs={len(y)}  n_drugs={batched.n_graphs}  n_pairs_attn={len(model.fusion.pairs)}  "
          f"params={n_params:,}  best_epoch={best_epoch}  prep={t_prep:.1f}s  fit={t_fit:.1f}s")
    print(f"Peak RSS: {peak_rss_gb():.2f} GB")

    return {
        "model": "CrossAttention",
        "omics": label,
        "n_modalities": len(modalities),
        "drug_rep": "molecular_graph",
        "n_pairs": len(y),
        "n_attention_pairs": len(model.fusion.pairs),
        "mean_only_rmse": floor,
        **{f"val_{k}": v for k, v in val.items()},
        **{f"test_{k}": v for k, v in test.items()},
        "params": n_params,
        "best_epoch": best_epoch,
        "fit_seconds": t_fit,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--omics",
        nargs="+",
        choices=ALL_MODALITIES,
        default=None,
        help=(
            "Run a single omics subset (>=2 modalities), e.g. --omics GE Proteomics "
            "to reproduce the E10 configuration. Defaults to all 4 multi-modality subsets."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.omics is not None and len(args.omics) < 2:
        raise SystemExit("--omics needs >=2 modalities: cross-attention attends between modalities.")
    subsets = [list(dict.fromkeys(args.omics))] if args.omics else multi_modality_subsets()

    t_start = time.perf_counter()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[device] {device}")

    # Built once and shared across runs: the molecules do not depend on the
    # omics subset, and parsing 621 SMILES per subset would be wasted work.
    graphs = build_drug_graphs()
    assert_fingerprint_population_parity(graphs)
    threshold = compute_shared_threshold()

    results = [run_one(modalities, graphs, threshold, device) for modalities in subsets]

    df = pd.DataFrame(results)
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_CSV, index=False)

    print("\n" + "=" * 100)
    print(f"MOLECULAR GRAPH MATRIX COMPLETE — {len(df)} runs in {time.perf_counter() - t_start:.1f}s")
    print("=" * 100)
    summary = df[["omics", "drug_rep", "n_pairs", "test_rmse", "test_pcc", "test_r2", "test_auc"]]
    print(summary.sort_values("test_rmse").to_string(index=False))
    print(f"\nSaved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
