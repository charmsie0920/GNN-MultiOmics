"""Backend wrapper for the final model (`final_xattn` in src/final_model/run_ablation.py).

`final_xattn` is `base+proteomics+std_targets+mol_graph+pair_features+cross_attention`:
the frozen MoGraphDRP-aligned base with proteomics as a third omics branch,
the omics fused by cross-attention instead of concatenation, a molecular
graph GCN beside the Morgan fingerprint, five hand-built target-vs-mutation
pair features into the predictor (`PairGraphDRP`), and per-drug standardised
ln(IC50) targets. It replaced `final` (the same model with concatenation) as
the served model on 2026-10-27.

Nothing is trained here. The five seeds were trained by the ablation ladder
with `--save-weights` (src/final_model/results/ablation_results.csv) and
their best-epoch weights are committed, so a run loads them and predicts.
The five are used together as an ensemble: the prediction is their mean, and
the spread between them is what the confidence score measures.

The model's inputs are the fixed reference data the checkpoints were trained
on, never the upload: the per-modality omics PCA files, the drug SMILES, the
heterogeneous graph's target and mutation edges, and the master response file
(which fixes the drug codes, the grouped split and the per-drug target
statistics). The upload decides which cell line and which drugs are scored,
and supplies their display names.

None of the science code is modified. The loading below rebuilds exactly the
rows `run_ablation.load_pairs` builds, without materialising the per-pair
fingerprint matrix, and each checkpoint's validation RMSE is checked against
the one recorded when it was trained.
"""

from __future__ import annotations

import sys
import threading
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from scipy.stats import pearsonr
from torch import nn

from backend.dataset_schema import UPLOADED_CSV_PATH
from backend.model_backends._common import attach_eval_metrics, rank_predictions
from backend.model_backends.base import ModelBackend
from backend.model_backends.registry import register

_REPO_ROOT = Path(__file__).resolve().parents[2]

# `src.*` is imported by package path, so the repo root has to be importable
# regardless of where uvicorn was launched from.
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

_CONFIG = "final_xattn"
# The ladder records the runs, and names their checkpoints, by the spelled-out config.
_RUN_NAME = "base+proteomics+std_targets+mol_graph+pair_features+cross_attention"
_SEEDS = (42, 43, 44, 45, 46)
_RESULTS_DIR = _REPO_ROOT / "src" / "final_model" / "results"
_WEIGHTS_CSV = _RESULTS_DIR / "ablation_results.csv"
_CHECKPOINT_DIR = _RESULTS_DIR / "ablation_results_checkpoints"
_GRAPH_PATH = _REPO_ROOT / "src" / "graph" / "hetero_graph.pt"
# The same file the upload endpoint writes to (backend/routers/dataset.py),
# so a run always scores the dataset the user just uploaded.
_TARGET_CSV = UPLOADED_CSV_PATH

_COL_CELL_LINE = "sanger_model_id"
_COL_DRUG = "drug_id"

# BatchNorm runs on its running statistics in eval mode, so the batch size
# changes speed only, never the output.
_EVAL_BATCH_SIZE = 2048
# A checkpoint whose validation RMSE differs from the recorded one by more
# than this was loaded against different data than it was trained on.
_RECORDED_RMSE_TOLERANCE = 1e-3
# Cap on how many validation pairs get sent to the UI for the scatter that
# stands in for a training curve -- far more than the widget needs otherwise.
_MAX_SCATTER_POINTS = 400
_SCATTER_SEED = 42


def _checkpoint_path(seed: int) -> Path:
    return _CHECKPOINT_DIR / f"{_RUN_NAME}__grouped__seed{seed}.pt"


def _missing_checkpoints() -> list[Path]:
    return [path for path in map(_checkpoint_path, _SEEDS) if not path.exists()]


def _require_repo_root() -> None:
    """The data loaders in `src/` read paths relative to the repo root."""
    if Path.cwd().resolve() != _REPO_ROOT:
        raise RuntimeError(
            f"The backend must be started from the repository root ({_REPO_ROOT}), "
            f"not {Path.cwd()}: run `python -m backend.main` from there."
        )


@dataclass
class _Session:
    """Everything a prediction needs that does not depend on the upload.

    Built once and kept for the life of the process: the reference data and
    the five models together are well under 200 MB, and rebuilding them would
    make every run, and every interpretation click, wait on the same work.
    """

    device: torch.device
    models: list[nn.Module]
    keys: list[str]
    omics: dict[str, np.ndarray]          # modality -> (n_cell_lines, 128)
    cell_row: dict[str, int]              # cell-line id -> row of `omics`
    fingerprints: np.ndarray              # (n_drugs, 2048), row = drug code
    drug_code: dict[str, int]             # drug id -> code the models were trained with
    train_cells: set[str]
    scaler: object                        # PerDrugTargetScaler fitted on the training fold
    adjacency: csr_matrix                 # STRING PPI over protein nodes
    targets_of: dict[str, np.ndarray]     # drug id -> target protein node indices
    mutated_of: dict[str, np.ndarray]     # cell-line id -> mutated protein node indices
    protein_ids: list[str]
    symbol_map: dict[str, str]
    val_rmse: float
    val_pcc: float
    val_scatter: list[dict]
    # (cell line, drug) -> pair features from `build_pair_features`, for the
    # pairs runs have scored.
    run_pairs: dict[tuple[str, str], np.ndarray] = field(default_factory=dict)


_session: _Session | None = None
# Endpoints are sync, so FastAPI serves them from a threadpool; a run and an
# interpretation click must not each build their own session.
_session_lock = threading.Lock()


@torch.no_grad()
def _predict_z(session: _Session, cell_rows: np.ndarray, drug_codes: np.ndarray, pair: np.ndarray) -> np.ndarray:
    """Per-drug z-scores from every model -> (n_models, n_rows)."""
    out = np.empty((len(session.models), len(drug_codes)), dtype=np.float64)
    for start in range(0, len(drug_codes), _EVAL_BATCH_SIZE):
        rows = slice(start, start + _EVAL_BATCH_SIZE)
        omics = {k: torch.from_numpy(session.omics[k][cell_rows[rows]]).to(session.device) for k in session.keys}
        fingerprints = torch.from_numpy(session.fingerprints[drug_codes[rows]]).to(session.device)
        codes = torch.from_numpy(drug_codes[rows]).to(session.device)
        features = torch.from_numpy(pair[rows]).to(session.device)
        for i, model in enumerate(session.models):
            out[i, rows] = model(omics, fingerprints, codes, features).cpu().numpy()
    return out


def _to_ln(session: _Session, z: np.ndarray, drug_codes: np.ndarray) -> np.ndarray:
    """(n_models, n_rows) z-scores -> ln(IC50), with the training fold's statistics."""
    return np.stack([session.scaler.inverse_transform(row, drug_codes) for row in z]).astype(np.float64)


def _rmse_pcc(pred: np.ndarray, true: np.ndarray) -> tuple[float, float]:
    rmse = float(np.sqrt(np.mean((pred - true) ** 2)))
    # pearsonr needs variance in both inputs; a degenerate split would raise.
    pcc = float(pearsonr(true, pred)[0]) if len(pred) > 1 and np.std(pred) > 0 else 0.0
    return rmse, pcc


def _sample_val_scatter(pred: np.ndarray, true: np.ndarray) -> list[dict]:
    """Down-sample (predicted, actual) ln(IC50) pairs into history-shaped dicts.

    Nothing is trained, so there is no per-epoch curve; the results page's
    training-curve panel shows real validation pairs instead.
    """
    if len(pred) > _MAX_SCATTER_POINTS:
        idx = np.random.default_rng(_SCATTER_SEED).choice(len(pred), size=_MAX_SCATTER_POINTS, replace=False)
        pred, true = pred[idx], true[idx]
    return [{"actual_ln_ic50": float(t), "predicted_ln_ic50": float(p)} for t, p in zip(true.tolist(), pred.tolist())]


def _build_session(log: Callable[[str], None], on_progress: Callable[[int, int], None]) -> _Session:
    from src.data.drug_graphs import build_drug_graphs, collate_drug_graphs
    from src.data.experiment_utils import (
        COL_TARGET,
        build_morgan_fingerprints,
        cell_row_indices,
        grouped_split,
        load_omics_subset,
        load_targets,
    )
    from src.data.pair_features import PPI, _load_pair_graph, build_pair_features
    from src.data.target_scaling import PerDrugTargetScaler
    from src.final_model.run_ablation import DROPOUT, resolve_config
    from src.interpretation.genes import load_symbol_map
    from src.models.pair_graph_drp import PairGraphDRP

    _require_repo_root()
    missing = _missing_checkpoints()
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} of the final model's checkpoints are missing, e.g. {missing[0].name}. "
            "They are committed under src/final_model/results/ablation_results_checkpoints/ and "
            "were saved by `python src/final_model/run_ablation.py --configs final_xattn "
            "--protocols grouped --seeds 42 43 44 45 46 --save-weights`."
        )

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    cfg = resolve_config(_CONFIG)
    keys = list(cfg["modalities"])
    total_steps = 4

    # -- reference data: the rows `run_ablation.load_pairs` builds ----------
    omics, cell_ids = load_omics_subset(keys)
    log(f"[omics] {len(cell_ids)} cell lines | " + " + ".join(f"{k} ({omics[k].shape[1]})" for k in keys))
    fingerprint_of = build_morgan_fingerprints()
    graphs = build_drug_graphs()
    y_used = load_targets(cell_ids)
    # The filter `build_drug_features(..., "fingerprint")` applies, without
    # building its (n_pairs, 2048) fingerprint matrix.
    y_used = y_used[y_used[_COL_DRUG].isin(fingerprint_of.keys())].reset_index(drop=True)
    codes, drug_levels = pd.factorize(y_used[_COL_DRUG], sort=True)
    codes = codes.astype(np.int64)
    batched = collate_drug_graphs(graphs, list(drug_levels))
    y = y_used[COL_TARGET].to_numpy(dtype=np.float32)
    groups = y_used[_COL_CELL_LINE].to_numpy()
    rows = cell_row_indices(cell_ids, y_used)
    on_progress(1, total_steps)

    pair_graph = _load_pair_graph(y_used, _GRAPH_PATH)
    graph = pair_graph.data
    n_proteins = int(graph["protein"].x.shape[0])
    ppi = graph[PPI].edge_index.numpy()
    adjacency = csr_matrix((np.ones(ppi.shape[1], dtype=np.int8), (ppi[0], ppi[1])), shape=(n_proteins, n_proteins))
    log(
        f"[graph] {len(cell_ids)} cell lines, {len(drug_levels)} drugs, {n_proteins} proteins "
        f"| {ppi.shape[1]} STRING edges for the pair features"
    )
    pair = build_pair_features(y_used, _GRAPH_PATH).matrix

    # -- the split and target statistics every checkpoint was trained under --
    train_idx, val_idx, _test_idx = grouped_split(groups)
    scaler = PerDrugTargetScaler().fit(y[train_idx], codes[train_idx])
    log(f"[data] {len(y)} reference pairs | grouped split: train {len(train_idx)}, val {len(val_idx)}")
    on_progress(2, total_steps)

    # -- the five seeds ------------------------------------------------------
    models: list[nn.Module] = []
    for seed in _SEEDS:
        model = PairGraphDRP(
            modalities=keys, drug_mode=cfg["drug_mode"], fusion=cfg["fusion"], head=cfg["head"],
            batched=batched, omics_in_dim=omics[keys[0]].shape[1], dropout=DROPOUT, pair_dim=pair.shape[1],
        )
        model.load_state_dict(torch.load(_checkpoint_path(seed), map_location="cpu"))
        models.append(model.to(device).eval())
    n_params = sum(p.numel() for p in models[0].parameters())
    log(f"[model] loaded {len(models)} checkpoints of '{_CONFIG}' (seeds {_SEEDS[0]}-{_SEEDS[-1]}, {n_params:,} params each)")
    on_progress(3, total_steps)

    session = _Session(
        device=device, models=models, keys=keys,
        omics=omics, cell_row={str(c): i for i, c in enumerate(cell_ids)},
        fingerprints=np.stack([fingerprint_of[d] for d in drug_levels]).astype(np.float32),
        drug_code={str(d): i for i, d in enumerate(drug_levels)},
        train_cells=set(groups[train_idx].astype(str)), scaler=scaler, adjacency=adjacency,
        targets_of=pair_graph.targets_of, mutated_of=pair_graph.mutated_of,
        protein_ids=[str(p) for p in graph["protein"].node_ids], symbol_map=load_symbol_map(),
        val_rmse=float("nan"), val_pcc=float("nan"), val_scatter=[],
    )

    # -- validation: each seed against its recorded figure, then the ensemble -
    val_pred = _to_ln(session, _predict_z(session, rows[val_idx], codes[val_idx], pair[val_idx]), codes[val_idx])
    runs = pd.read_csv(_WEIGHTS_CSV)
    recorded = runs[(runs["config"] == _RUN_NAME) & (runs["protocol"] == "grouped")].set_index("seed")["val_rmse"]
    for seed, pred in zip(_SEEDS, val_pred):
        rmse, _ = _rmse_pcc(pred, y[val_idx])
        expected = recorded.get(seed, float("nan"))
        if not abs(rmse - expected) <= _RECORDED_RMSE_TOLERANCE:
            log(f"[eval] WARNING seed {seed}: validation RMSE {rmse:.4f}, recorded {expected:.4f} -- "
                "the reference data differs from what this checkpoint was trained on")
    ensemble = val_pred.mean(axis=0)
    session.val_rmse, session.val_pcc = _rmse_pcc(ensemble, y[val_idx])
    session.val_scatter = _sample_val_scatter(ensemble, y[val_idx].astype(np.float64))
    per_seed = ", ".join(f"{_rmse_pcc(p, y[val_idx])[0]:.4f}" for p in val_pred)
    log(f"[eval] validation RMSE per seed: {per_seed} (match the recorded runs)")
    on_progress(4, total_steps)
    return session


def _get_session(log: Callable[[str], None], on_progress: Callable[[int, int], None]) -> _Session:
    global _session
    with _session_lock:
        if _session is None:
            _session = _build_session(log, on_progress)
        else:
            log("[model] reusing the loaded final model and reference data")
        return _session


def _pair_vector(hop_rows: dict[int, np.ndarray], targets: np.ndarray, mutated: np.ndarray) -> np.ndarray:
    """The five pair features for one pair, given each target's row of hop counts.

    The same formula as `src/data/pair_features.build_pair_features`, for one
    pair at a time, so a single gene can be taken out of either set. Its
    output for the unchanged sets is checked against that function's.
    """
    has_target, has_mutation = len(targets) > 0, len(mutated) > 0
    hops, near = np.inf, 0
    if has_target and has_mutation:
        to_target = np.min([hop_rows[int(t)] for t in targets], axis=0)[mutated]
        hops, near = to_target.min(), int((to_target <= 1).sum())
    return np.array(
        [has_target, has_mutation, hops == 0, 1.0 / (1.0 + hops), np.log1p(near)], dtype=np.float32
    )


class FinalModelBackend(ModelBackend):
    """The project's final model, as five trained seeds used as an ensemble."""

    name = "final_model"
    # Genes enter this model by name only through the pair features (the
    # omics are PCA components), and those can be attributed exactly by
    # taking one gene out at a time. See `explain`.
    supports_interpretation = True

    @property
    def expected_duration_seconds(self) -> dict[str, float]:
        """The first run builds the session; later runs only predict."""
        if _session is not None:
            return {"cuda": 5.0, "cpu": 8.0}
        return {"cuda": 45.0, "cpu": 120.0}

    def explain(self, run_id: str, target_cell_line: str, drug_id: str, top_k: int = 50) -> list[dict]:
        """Each gene's effect on one prediction, by removing it from the pair.

        The candidates are the drug's target proteins and the cell line's
        driver-mutated proteins: the only genes the model sees by name. For
        each one the five pair features are recomputed as if the gene were in
        neither set, and the gene's score is the ensemble's prediction with
        it minus the prediction without it, in ln(IC50). The omics inputs are
        unchanged, so a mutation still reaches the model through Mut_CNV;
        what is measured is the pair-feature route only.
        """
        from src.interpretation.genes import annotate_genes

        session = _get_session(lambda _line: None, lambda _current, _total: None)
        cell, drug = str(target_cell_line), str(drug_id)
        if cell not in session.cell_row:
            raise ValueError(f"Cell line {cell!r} has no omics profile, so it cannot be explained.")
        if drug not in session.drug_code:
            raise ValueError(f"Drug {drug!r} is not one the final model was trained on.")

        empty = np.empty(0, dtype=np.int64)
        targets = session.targets_of.get(drug, empty)
        mutated = session.mutated_of.get(cell, empty)
        if len(targets) == 0:
            raise ValueError(
                "This drug has no known target protein, so no gene enters its prediction by name. "
                "The prediction rests on the omics profile and the drug's structure."
            )

        hops = shortest_path(session.adjacency, unweighted=True, indices=targets)
        hop_rows = {int(t): hops[i] for i, t in enumerate(targets)}

        actual = _pair_vector(hop_rows, targets, mutated)
        reference = session.run_pairs.get((cell, drug))
        if reference is not None and not np.allclose(actual, reference):
            raise RuntimeError("Recomputed pair features disagree with src/data/pair_features.py.")

        candidates = np.union1d(targets, mutated)
        variants = [actual] + [
            _pair_vector(hop_rows, targets[targets != gene], mutated[mutated != gene]) for gene in candidates
        ]
        n = len(variants)
        code = session.drug_code[drug]
        drug_codes = np.full(n, code, dtype=np.int64)
        z = _predict_z(session, np.full(n, session.cell_row[cell]), drug_codes, np.stack(variants))
        pred = _to_ln(session, z, drug_codes).mean(axis=0)
        scores = pred[0] - pred[1:]

        moved = np.abs(scores) > 1e-6
        if not moved.any():
            raise ValueError(
                "None of this drug's targets or this cell line's mutated genes changes the prediction: "
                "no mutation lies close enough to a target in the interaction network."
            )
        node_ids = session.protein_ids
        return annotate_genes(
            scores[moved].tolist(),
            [node_ids[i] for i in candidates[moved]],
            session.symbol_map,
            top_k,
            driver_proteins={node_ids[i] for i in mutated},
            target_proteins={node_ids[i] for i in targets},
        )

    def run(
        self,
        log: Callable[[str], None],
        target_cell_line: str,
        on_results: Callable[[list[dict]], None],
        on_progress: Callable[[int, int], None],
        on_training_history: Callable[[list[dict]], None],
    ) -> None:
        from src.data.pair_features import build_pair_features

        log(f"[setup] device: {'cuda' if torch.cuda.is_available() else 'cpu'}")
        session = _get_session(log, on_progress)
        log(f"[eval] ensemble validation RMSE (ln IC50) = {session.val_rmse:.4f}, PCC = {session.val_pcc:.4f} "
            "-- used as the confidence scale")
        on_training_history(attach_eval_metrics(list(session.val_scatter), session.val_rmse, session.val_pcc))

        cell = str(target_cell_line)
        if cell not in session.cell_row:
            raise ValueError(f"Cell line {cell!r} has no omics profile, so it cannot be scored.")

        uploaded = pd.read_csv(_TARGET_CSV, dtype=str, keep_default_na=False)
        drug_info, cell_info = _metadata(uploaded)
        panel = sorted({d for d in uploaded[_COL_DRUG].unique() if d in session.drug_code}, key=session.drug_code.get)
        if not panel:
            raise ValueError("None of the drugs in the uploaded dataset are ones the final model was trained on.")
        n_uploaded = uploaded[_COL_DRUG].nunique()
        log(f"[data] {len(panel)} of {n_uploaded} uploaded drugs are in the model's {len(session.drug_code)}")

        info = cell_info.get(cell, {})
        log(
            f"[cell-line] {cell} -> {info.get('cell_line_name', 'unknown')} "
            f"({info.get('tissue', 'unknown')}, {info.get('cancer_type', 'unknown')})"
        )
        if cell in session.train_cells:
            log("[cell-line] note: this cell line was in the training fold, so measured responses "
                "for it were seen in training")
        else:
            log("[cell-line] held out: the model never saw this cell line's responses")

        log(f"Running inference for target cell line: {cell}")
        frame = pd.DataFrame({_COL_CELL_LINE: cell, _COL_DRUG: panel})
        pair = build_pair_features(frame, _GRAPH_PATH).matrix
        # What `explain` checks its own recomputation against.
        session.run_pairs.update({(cell, d): pair[i] for i, d in enumerate(panel)})
        drug_codes = np.array([session.drug_code[d] for d in panel], dtype=np.int64)
        z = _predict_z(session, np.full(len(panel), session.cell_row[cell]), drug_codes, pair)
        pred = _to_ln(session, z, drug_codes)
        mean, std = pred.mean(axis=0), pred.std(axis=0)
        raw = [
            {"drug_id": d, "ln_ic50_mean": float(mean[i]), "ln_ic50_std": float(std[i])}
            for i, d in enumerate(panel)
        ]
        log(f"[predict] scored {len(raw)} drugs with the {len(session.models)}-seed ensemble")
        on_results(rank_predictions(raw, session.val_rmse, drug_info))


def _metadata(uploaded: pd.DataFrame) -> tuple[dict[str, dict], dict[str, dict]]:
    """Drug and cell-line display names from the upload, keyed by raw id.

    Optional columns: an upload with only the three required ones still runs,
    and its drugs are shown by id.
    """
    def table(key: str, fields: list[str]) -> dict[str, dict]:
        present = [f for f in fields if f in uploaded.columns]
        rows = uploaded[[key, *present]].drop_duplicates(key)
        return {row[key]: {f: row.get(f, "") for f in fields} for row in rows.to_dict("records")}

    return (
        table(_COL_DRUG, ["drug_name", "putative_target", "pathway_name"]),
        table(_COL_CELL_LINE, ["cell_line_name", "tissue", "cancer_type"]),
    )


register(FinalModelBackend())
