"""Backend wrapper for the heterogeneous GNN (src/models/train/hetero_gnn.py).

Drives `HeteroIC50GNN` completely unmodified — imports the class, loads the
prebuilt heterogeneous graph, and does all orchestration here in the backend.
None of the science scripts (`src/models/train/hetero_gnn.py`,
`src/models/train/05_train.py`, `src/models/train/dataset.py`) are touched;
this module replicates what they do, in memory, so a UI run never writes to
the repo.

Unlike the cross-attention backend, which must retrain from scratch on every
run, this one loads `models/checkpoints/best_hetero_gnn.pt` when it is
present and predicts immediately -- a run drops from minutes to seconds. If
the checkpoint is missing, it falls back to training the model the same way
`05_train.py` does, streaming real epoch metrics to the UI.
"""

from __future__ import annotations

import gc
import sys
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr
from sklearn.model_selection import GroupShuffleSplit
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

from backend.model_backends._common import enable_mc_dropout, rank_predictions
from backend.model_backends.base import ModelBackend
from backend.model_backends.registry import register

_REPO_ROOT = Path(__file__).resolve().parents[2]

# `src.models.train.hetero_gnn` is imported by package path, so the repo root
# has to be importable regardless of where uvicorn was launched from.
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

_GRAPH_PATH = _REPO_ROOT / "src" / "graph" / "hetero_graph.pt"
_CHECKPOINT_PATH = _REPO_ROOT / "models" / "checkpoints" / "best_hetero_gnn.pt"
# The same file the upload endpoint writes to (backend/routers/dataset.py),
# so a run always scores the dataset the user just uploaded.
_TARGET_CSV = _REPO_ROOT / "data" / "processed" / "aligned" / "gdsc2_response_master.csv"

_COL_CELL_LINE = "sanger_model_id"
_COL_DRUG = "drug_id"
_COL_TARGET = "ln_ic50"

# MC-dropout forward passes per (cell-line, drug) pair used to derive the
# confidence score -- see _predict_drug_panel below.
_MC_SAMPLES = 30

# Training hyperparameters, mirroring src/models/train/05_train.py so the
# fallback path reproduces the checkpoint's recipe rather than inventing a
# new one.
_HIDDEN_DIM = 128
_LR = 1e-3
_WEIGHT_DECAY = 1e-5
_MAX_EPOCHS = 200
_PATIENCE = 15
_TRAIN_BATCH_SIZE = 1024
_EVAL_BATCH_SIZE = 2048
_SEED = 42


def _load_graph(device: torch.device) -> dict:
    """Load the prebuilt HeteroData graph and shape it for `HeteroIC50GNN`."""
    if not _GRAPH_PATH.exists():
        raise FileNotFoundError(
            f"Heterogeneous graph not found at {_GRAPH_PATH}. It is committed to the repo -- "
            "check out the file or rebuild it with src/data/03_graph_construction.py."
        )

    graph = torch.load(_GRAPH_PATH, weights_only=False)

    x_dict = {node_type: graph[node_type].x.to(device, dtype=torch.float32) for node_type in graph.node_types}
    edge_index_dict = {
        edge_type: graph[edge_type].edge_index.to(device, dtype=torch.long) for edge_type in graph.edge_types
    }

    # The graph stores only the forward direction of these two relations, but
    # HeteroIC50GNN declares convolutions for their reverses too (unlike the
    # tracked HeteroGNN, it does not build them itself). Without these, drug
    # and cell_line nodes receive no messages back from the protein layer and
    # every prediction collapses toward a constant. Same construction as
    # src/models/train/05_train.py.
    src, dst = edge_index_dict[("drug", "targets", "protein")]
    edge_index_dict[("protein", "rev_targets", "drug")] = torch.stack([dst, src])

    src, dst = edge_index_dict[("cell_line", "has_mutation", "protein")]
    edge_index_dict[("protein", "rev_has_mutation", "cell_line")] = torch.stack([dst, src])

    return {
        "x_dict": x_dict,
        "edge_index_dict": edge_index_dict,
        # Positional lists -- index i in these lists *is* node index i in the
        # graph. This is the entire id-mapping mechanism; there is no separate
        # encoder artifact that could drift out of sync.
        "cell_to_idx": {str(cell_id): i for i, cell_id in enumerate(graph["cell_line"].node_ids)},
        "drug_to_idx": {str(drug_id): i for i, drug_id in enumerate(graph["drug"].node_ids)},
        "num_proteins": int(x_dict["protein"].shape[0]),
    }


def _load_metadata() -> tuple[dict[str, dict], dict[str, dict]]:
    """Drug and cell-line descriptive metadata, read from the uploaded CSV.

    Keyed by the same raw id strings used for the graph lookup, so names and
    predictions can never be misaligned by a separate mapping step.
    """
    meta = pd.read_csv(
        _TARGET_CSV,
        usecols=[
            _COL_DRUG,
            _COL_CELL_LINE,
            "drug_name",
            "putative_target",
            "pathway_name",
            "cell_line_name",
            "tissue",
            "cancer_type",
        ],
        dtype=str,
        keep_default_na=False,
    )
    drug_meta = meta[[_COL_DRUG, "drug_name", "putative_target", "pathway_name"]].drop_duplicates(_COL_DRUG)
    drug_info = {
        row[_COL_DRUG]: {
            "drug_name": row["drug_name"],
            "putative_target": row["putative_target"],
            "pathway_name": row["pathway_name"],
        }
        for row in drug_meta.to_dict("records")
    }
    cell_meta = meta[[_COL_CELL_LINE, "cell_line_name", "tissue", "cancer_type"]].drop_duplicates(_COL_CELL_LINE)
    cell_info = {
        row[_COL_CELL_LINE]: {
            "cell_line_name": row["cell_line_name"],
            "tissue": row["tissue"],
            "cancer_type": row["cancer_type"],
        }
        for row in cell_meta.to_dict("records")
    }
    return drug_info, cell_info


def _build_pairs(graph_maps: dict, log: Callable[[str], None]) -> dict:
    """Map the uploaded response CSV onto graph node indices and split it.

    Replicates src/models/train/dataset.py in memory (that script writes
    data/raw/aligned_ic50_pairs.csv; a served run must not mutate the repo).
    The split is grouped by cell line so no cell line appears in more than one
    split -- a plain random split would leak the same cell line's other
    drug-response rows into both train and validation.
    """
    frame = pd.read_csv(_TARGET_CSV, usecols=[_COL_CELL_LINE, _COL_DRUG, _COL_TARGET])
    total_rows = len(frame)

    cell_idx = frame[_COL_CELL_LINE].astype(str).map(graph_maps["cell_to_idx"])
    drug_idx = frame[_COL_DRUG].astype(str).map(graph_maps["drug_to_idx"])
    frame = frame.assign(cell_idx=cell_idx, drug_idx=drug_idx).dropna(subset=["cell_idx", "drug_idx"])
    frame["cell_idx"] = frame["cell_idx"].astype(int)
    frame["drug_idx"] = frame["drug_idx"].astype(int)

    if frame.empty:
        raise ValueError(
            "No rows in the uploaded dataset could be matched to the graph. Check that "
            f"{_COL_CELL_LINE}/{_COL_DRUG} values match the graph's cell-line and drug nodes."
        )

    log(f"[data] matched {len(frame)} / {total_rows} response pairs to graph nodes")

    groups = frame[_COL_CELL_LINE].to_numpy()
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.3, random_state=_SEED)
    train_idx, rest_idx = next(splitter.split(frame, groups=groups))
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.5, random_state=_SEED)
    rel_val, rel_test = next(splitter.split(frame.iloc[rest_idx], groups=groups[rest_idx]))
    val_idx, test_idx = rest_idx[rel_val], rest_idx[rel_test]

    log(f"[data] split -> train {len(train_idx)} | val {len(val_idx)} | test {len(test_idx)} (grouped by cell line)")

    return {
        "frame": frame,
        "train": frame.iloc[train_idx],
        "val": frame.iloc[val_idx],
        "test": frame.iloc[test_idx],
        "drug_ids_present": frame[_COL_DRUG].astype(str).unique(),
    }


def _make_loader(frame: pd.DataFrame, batch_size: int, shuffle: bool = False, drop_last: bool = False) -> DataLoader:
    dataset = TensorDataset(
        torch.tensor(frame["cell_idx"].to_numpy(), dtype=torch.long),
        torch.tensor(frame["drug_idx"].to_numpy(), dtype=torch.long),
        torch.tensor(frame[_COL_TARGET].to_numpy(), dtype=torch.float32),
    )
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, drop_last=drop_last)


@torch.no_grad()
def _evaluate(
    model: nn.Module, loader: DataLoader, graph_maps: dict, device: torch.device
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """Deterministic RMSE, Pearson r, and the raw (predicted, actual) arrays over `loader`.

    The raw arrays are needed by `run()` to build the validation-pair scatter
    that stands in for a training curve when the backend loaded a checkpoint
    instead of training (see `_sample_val_scatter`).
    """
    model.eval()
    preds, targets = [], []
    for cell_idx, drug_idx, batch_target in loader:
        out = model(graph_maps["x_dict"], graph_maps["edge_index_dict"], cell_idx.to(device), drug_idx.to(device))
        preds.append(out.cpu().numpy())
        targets.append(batch_target.numpy())

    pred = np.concatenate(preds)
    true = np.concatenate(targets)
    rmse = float(np.sqrt(np.mean((pred - true) ** 2)))
    # pearsonr needs variance in both inputs; a degenerate split would raise.
    pcc = float(pearsonr(true, pred)[0]) if len(pred) > 1 and np.std(pred) > 0 else 0.0
    return rmse, pcc, pred, true


# Cap on how many validation pairs get sent to the UI for the scatter
# fallback -- the validation split can be tens of thousands of rows, far more
# than a small custom-painted widget needs to show the calibration pattern.
_MAX_SCATTER_POINTS = 400


def _sample_val_scatter(pred: np.ndarray, true: np.ndarray, seed: int = _SEED) -> list[dict]:
    """Down-sample (predicted, actual) ln_ic50 pairs into history-shaped dicts.

    Used only when no real per-epoch history exists (checkpoint mode), so the
    results page's training-curve panel still has real, run-specific data to
    plot instead of sitting empty.
    """
    n = len(pred)
    if n > _MAX_SCATTER_POINTS:
        idx = np.random.default_rng(seed).choice(n, size=_MAX_SCATTER_POINTS, replace=False)
        pred, true = pred[idx], true[idx]
    return [
        {"actual_ln_ic50": float(t), "predicted_ln_ic50": float(p)} for t, p in zip(true.tolist(), pred.tolist())
    ]


def _train_model(
    model: nn.Module,
    pairs: dict,
    graph_maps: dict,
    device: torch.device,
    log: Callable[[str], None],
    on_progress: Callable[[int, int], None],
    history: list[dict],
) -> nn.Module:
    """Train `model` following src/models/train/05_train.py's recipe.

    Because this loop lives here rather than in a wrapped script, epoch
    metrics are reported to the UI directly instead of being scraped back out
    of printed text the way the cross-attention backend has to do.
    """
    # drop_last so the head's BatchNorm never sees a 1-sample final batch.
    train_loader = _make_loader(pairs["train"], _TRAIN_BATCH_SIZE, shuffle=True, drop_last=True)
    val_loader = _make_loader(pairs["val"], _EVAL_BATCH_SIZE)

    optimizer = torch.optim.AdamW(model.parameters(), lr=_LR, weight_decay=_WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, factor=0.5, patience=5)
    criterion = nn.MSELoss()

    best_rmse = float("inf")
    best_state: dict | None = None
    epochs_without_improvement = 0

    for epoch in range(1, _MAX_EPOCHS + 1):
        model.train()
        running_loss = 0.0
        for cell_idx, drug_idx, batch_target in train_loader:
            batch_target = batch_target.to(device)
            optimizer.zero_grad()
            out = model(graph_maps["x_dict"], graph_maps["edge_index_dict"], cell_idx.to(device), drug_idx.to(device))
            loss = criterion(out, batch_target)
            loss.backward()
            optimizer.step()
            running_loss += loss.item() * len(batch_target)

        train_loss = running_loss / max(1, len(train_loader.dataset))
        val_rmse, val_pcc, _val_preds, _val_true = _evaluate(model, val_loader, graph_maps, device)
        scheduler.step(val_rmse)

        history.append(
            {"epoch": epoch, "train_loss": float(train_loss), "val_rmse": float(val_rmse), "val_pcc": float(val_pcc)}
        )
        log(f"[epoch {epoch:3d}] train_loss={train_loss:.4f} val_rmse={val_rmse:.4f} val_pcc={val_pcc:.4f}")
        on_progress(epoch, _MAX_EPOCHS)

        if val_rmse < best_rmse:
            best_rmse = val_rmse
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1
            if epochs_without_improvement >= _PATIENCE:
                log(f"[train] early stopping at epoch {epoch} (best val_rmse={best_rmse:.4f})")
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    return model


def _load_or_train(
    pairs: dict,
    graph_maps: dict,
    device: torch.device,
    log: Callable[[str], None],
    on_progress: Callable[[int, int], None],
    history: list[dict],
) -> nn.Module:
    from src.models.train.hetero_gnn import HeteroIC50GNN

    torch.manual_seed(_SEED)
    model = HeteroIC50GNN(num_proteins=graph_maps["num_proteins"], hidden_dim=_HIDDEN_DIM).to(device)

    if _CHECKPOINT_PATH.exists():
        # The convolutions are non-lazy, so the state_dict loads directly with
        # no dummy forward pass needed to materialize shapes first.
        state = torch.load(_CHECKPOINT_PATH, map_location=device)
        model.load_state_dict(state)
        log(f"[model] loaded pretrained checkpoint: {_CHECKPOINT_PATH.name} ({sum(p.numel() for p in model.parameters()):,} params)")
        log("[model] skipping training -- predicting directly from the trained weights")
        on_progress(_MAX_EPOCHS, _MAX_EPOCHS)
        return model

    log(f"[model] no checkpoint at {_CHECKPOINT_PATH.name} -- training from scratch")
    on_progress(0, _MAX_EPOCHS)
    return _train_model(model, pairs, graph_maps, device, log, on_progress, history)


@torch.no_grad()
def _predict_drug_panel(
    model: nn.Module,
    pairs: dict,
    graph_maps: dict,
    device: torch.device,
    target_cell_line: str,
    mc_samples: int,
) -> list[dict]:
    """Predict ln_ic50 (mean + MC-dropout std) for every candidate drug, for one cell line."""
    cell_idx = graph_maps["cell_to_idx"].get(str(target_cell_line))
    if cell_idx is None:
        raise ValueError(f"Cell line {target_cell_line!r} has no node in the graph, so it cannot be scored.")

    # Only drugs that are both in the uploaded dataset and present as graph
    # nodes can be scored -- a drug with no node has no fingerprint features.
    drug_to_idx = graph_maps["drug_to_idx"]
    panel = [(drug_id, drug_to_idx[drug_id]) for drug_id in pairs["drug_ids_present"] if drug_id in drug_to_idx]
    if not panel:
        raise ValueError("None of the drugs in the uploaded dataset have nodes in the graph.")

    drug_ids = [drug_id for drug_id, _ in panel]
    drug_index = torch.tensor([idx for _, idx in panel], dtype=torch.long, device=device)
    cell_index = torch.full((len(panel),), cell_idx, dtype=torch.long, device=device)

    # Dropout live, BatchNorm on running stats. Dropout sits only in the
    # model's predictor head, but forward() re-runs message passing on every
    # call and the model file is off-limits, so each sample is a full forward.
    enable_mc_dropout(model)
    try:
        samples = np.empty((mc_samples, len(panel)), dtype=np.float64)
        for i in range(mc_samples):
            samples[i] = model(graph_maps["x_dict"], graph_maps["edge_index_dict"], cell_index, drug_index).cpu().numpy()
    finally:
        model.eval()

    mean = samples.mean(axis=0)
    std = samples.std(axis=0)
    return [
        {"drug_id": str(drug_id), "ln_ic50_mean": float(mean[i]), "ln_ic50_std": float(std[i])}
        for i, drug_id in enumerate(drug_ids)
    ]


class HeteroGNNBackend(ModelBackend):
    """Heterogeneous GNN over the cell-line / drug / protein graph.

    Cell lines carry 384-dim concatenated omics PCA features, drugs carry
    2048-bit Morgan fingerprints, and proteins are linked by STRING
    interactions, drug-target edges, and cell-line mutation edges. Two
    HeteroConv layers propagate across those relations before an MLP head
    regresses ln(IC50) for a (cell line, drug) pair.
    """

    name = "hetero_gnn"

    @property
    def expected_duration_seconds(self) -> dict[str, float]:
        """Runtime differs by two orders of magnitude between the two paths,
        so the estimate has to be resolved at run time rather than fixed as a
        class attribute -- a 20-minute estimate on a 15-second run would leave
        the UI's progress bar pinned near zero for the whole run."""
        if _CHECKPOINT_PATH.exists():
            return {"cuda": 20.0, "cpu": 45.0}
        return {"cuda": 300.0, "cpu": 1200.0}

    def run(
        self,
        log: Callable[[str], None],
        target_cell_line: str,
        on_results: Callable[[list[dict]], None],
        on_progress: Callable[[int, int], None],
        on_training_history: Callable[[list[dict]], None],
    ) -> None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        log(f"[setup] device: {device}")
        history: list[dict] = []
        graph_maps: dict | None = None
        model: nn.Module | None = None

        try:
            log(f"[graph] loading heterogeneous graph from {_GRAPH_PATH.name}")
            graph_maps = _load_graph(device)
            log(
                f"[graph] {len(graph_maps['cell_to_idx'])} cell lines, "
                f"{len(graph_maps['drug_to_idx'])} drugs, {graph_maps['num_proteins']} proteins"
            )

            pairs = _build_pairs(graph_maps, log)
            drug_info, cell_info = _load_metadata()

            model = _load_or_train(pairs, graph_maps, device, log, on_progress, history)

            # Held-out validation error, both as a quality readout and as the
            # scale the confidence score is calibrated against (see
            # `rank_predictions`). One cheap forward pass in checkpoint mode.
            val_rmse, val_pcc, val_preds, val_true = _evaluate(
                model, _make_loader(pairs["val"], _EVAL_BATCH_SIZE), graph_maps, device
            )
            log(f"[eval] validation RMSE (ln IC50) = {val_rmse:.4f}, PCC = {val_pcc:.4f} -- used as the confidence scale")

            if not history:
                # No epochs ran (checkpoint mode), so there is no curve to
                # show on the results page's training-curve panel. Give it
                # real data anyway: a sample of predicted-vs-actual ln_ic50
                # pairs from the same validation pass above.
                history = _sample_val_scatter(val_preds, val_true)
                log(f"[eval] no training curve to show (loaded from checkpoint) -- sending {len(history)} validation pairs instead")
            on_training_history(history)

            info = cell_info.get(target_cell_line, {})
            log(
                f"[cell-line] {target_cell_line} -> {info.get('cell_line_name', 'unknown')} "
                f"({info.get('tissue', 'unknown')}, {info.get('cancer_type', 'unknown')})"
            )
            log(f"Running inference for target cell line: {target_cell_line}")
            raw_predictions = _predict_drug_panel(
                model, pairs, graph_maps, device, target_cell_line, mc_samples=_MC_SAMPLES
            )
            log(f"[predict] scored {len(raw_predictions)} drugs over {_MC_SAMPLES} MC-dropout samples")
            on_results(rank_predictions(raw_predictions, val_rmse, drug_info))
        finally:
            # The backend process stays alive across runs, and the graph
            # tensors alone are ~21 MB of features plus half a million edges,
            # so drop everything and force a real reclaim before returning
            # rather than letting the next run's allocations do it.
            model = None
            graph_maps = None
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            gc.collect()


register(HeteroGNNBackend())
