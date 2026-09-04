"""Shared loaders/split/eval helpers for the omics x drug-rep x architecture experiment matrix.

`rf_baseline.py`, `mlp_baseline.py`, and `cross_attention_baseline.py` each
duplicated the same `load_targets`/`grouped_split`/`evaluate`/`peak_rss_gb`
boilerplate near-verbatim. This module centralizes that logic so the ~40-run
experiment matrix (docs/plan/experiment_matrix_plan.md) doesn't copy-paste it
a fourth, fifth, and sixth time, and adds the two axes none of the existing
scripts support yet: loading an arbitrary *subset* of the three omics
modalities, and a Morgan-fingerprint drug representation (in addition to the
existing one-hot).

Every matrix cell (RF/MLP/cross-attention x omics-subset x drug-rep) is
expected to call these functions so every run sits on an identical
`GroupShuffleSplit`-by-cell-line 70/15/15 partition with `random_state=42`,
which is what makes the final comparison table valid.
"""

from __future__ import annotations

import csv
import platform
from pathlib import Path
from typing import Dict, List, Literal, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import f1_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit
from scipy.stats import spearmanr

# --- shared paths -----------------------------------------------------------
DATA_DIR = Path("data/processed")
ALIGNED_DIR = DATA_DIR / "aligned"

GE_KEY = "GE"
MUT_CNV_KEY = "Mut_CNV"
PROTEOMICS_KEY = "Proteomics"

# File -> modality-key mapping per src/data/00_run_preprocessing.py:31-36
OMICS_FILES: Dict[str, Path] = {
    GE_KEY: DATA_DIR / "transcriptomics_pca.csv",
    MUT_CNV_KEY: DATA_DIR / "genomics_pca.csv",
    PROTEOMICS_KEY: DATA_DIR / "proteomics_pca.csv",
}

DRUG_SMILES_PATH = Path("data/raw/pubchem/gdsc_drug_smiles.csv")
TARGET_CSV = ALIGNED_DIR / "gdsc2_response_master.csv"

COL_CELL_LINE = "sanger_model_id"
COL_DRUG = "drug_id"
COL_TARGET = "ln_ic50"

TEST_FRAC = 0.15
VAL_FRAC = 0.15
RANDOM_STATE = 42
DTYPE = np.float32

MORGAN_RADIUS = 2
MORGAN_FP_SIZE = 2048

DrugRepMode = Literal["onehot", "fingerprint"]

# The three arms every matrix script sweeps. `onehot_restricted` is the
# population control: one-hot drug identity, but evaluated on exactly the
# pairs that fingerprint mode can cover, so drug-representation comparisons
# are measured on identical rows instead of across two different populations.
DRUG_ARMS: List[Tuple[str, DrugRepMode, bool]] = [
    ("onehot", "onehot", False),
    ("onehot_restricted", "onehot", True),
    ("fingerprint", "fingerprint", False),
]


# --- system / reporting ------------------------------------------------------
def peak_rss_gb() -> float:
    """Peak resident set size. Returns 0.0 on Windows to prevent crashes."""
    if platform.system() == "Windows":
        return 0.0

    import resource

    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw / 1024**3 if platform.system() == "Darwin" else raw / 1024**2


# --- omics loading (subset-capable) ------------------------------------------
def load_omics_subset(modalities: Sequence[str]) -> Tuple[Dict[str, np.ndarray], pd.Index]:
    """Load the requested subset of the 3 per-modality PCA CSVs.

    `modalities` must be a non-empty subset of {GE_KEY, MUT_CNV_KEY,
    PROTEOMICS_KEY}. All requested CSVs share the identical 532-cell-line
    index (asserted, not assumed), so any subset yields the same row
    population regardless of how many modalities are requested.
    """
    if not modalities:
        raise ValueError("modalities must be a non-empty subset of OMICS_FILES")

    arrays: Dict[str, np.ndarray] = {}
    index: pd.Index | None = None
    for key in modalities:
        path = OMICS_FILES[key]
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
    assert index is not None
    return arrays, index


def concat_omics(arrays: Dict[str, np.ndarray], modalities: Sequence[str]) -> np.ndarray:
    """Concatenate the requested modalities column-wise, in a fixed order."""
    return np.concatenate([arrays[key] for key in modalities], axis=1)


# --- targets ------------------------------------------------------------------
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

    print(
        f"[targets] {n_raw} rows -> {n_finite} finite -> {n_matched} with omics "
        f"-> {len(y)} unique pairs"
    )
    if y.empty:
        raise ValueError("No target rows survived filtering; check ID formatting.")
    return y.reset_index(drop=True)


def compute_shared_threshold() -> float:
    """Median ln_ic50 over all target rows, used as the AUC/F1 binarization cutoff."""
    y_all = pd.read_csv(TARGET_CSV, usecols=[COL_TARGET])
    y_all[COL_TARGET] = pd.to_numeric(y_all[COL_TARGET], errors="coerce")
    threshold = float(y_all[COL_TARGET].median(skipna=True))
    print(f"[threshold] median ln_ic50 across all targets = {threshold:.4f}")
    return threshold


# --- drug representations ------------------------------------------------------
def build_morgan_fingerprints() -> Dict[str, np.ndarray]:
    """Parse cached SMILES into 2048-bit Morgan fingerprints, keyed by GDSC drug_id.

    Drugs with missing/unparseable SMILES are simply absent from the returned
    dict (not zero-filled) -- callers must drop pairs referencing an absent
    drug_id rather than silently zero-featurizing them.
    """
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    if not DRUG_SMILES_PATH.exists():
        raise FileNotFoundError(
            f"Missing {DRUG_SMILES_PATH}. Run src/data/fetch_drug_smiles.py first."
        )

    with DRUG_SMILES_PATH.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=MORGAN_RADIUS, fpSize=MORGAN_FP_SIZE)

    fingerprints: Dict[str, np.ndarray] = {}
    unparseable = 0
    for row in rows:
        smiles = row["canonical_smiles"]
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            unparseable += 1
            continue
        fp = mfpgen.GetFingerprint(mol)
        arr = np.zeros((MORGAN_FP_SIZE,), dtype=DTYPE)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fingerprints[row["drug_id"]] = arr

    print(
        f"[fingerprints] {len(fingerprints)} / {len(rows)} drug_ids resolved to a "
        f"Morgan fingerprint ({unparseable} missing/unparseable SMILES dropped)"
    )
    return fingerprints


def build_drug_features(
    y: pd.DataFrame, mode: DrugRepMode, restrict_to_fingerprintable: bool = False
) -> Tuple[pd.DataFrame, np.ndarray, int, List[str]]:
    """Build per-row drug features for `y`, returning the (possibly filtered) y.

    - `mode="onehot"`: one-hot drug identity. With
      `restrict_to_fingerprintable=False` every pair is kept (the historical
      baseline population, ~134,764 pairs). With it `True`, pairs are first
      filtered to the drugs fingerprint mode can cover, giving the population
      control arm -- the *same* rows fingerprint mode sees, so the
      drug-representation comparison isn't confounded by a differing row set.
    - `mode="fingerprint"`: rows whose `drug_id` has no resolved Morgan
      fingerprint (unresolved PubChem name or unparseable SMILES) are always
      dropped, since there is no feature vector to give them.

    Returns (filtered_y, drug_feature_matrix, feature_dim, drug_levels).
    """
    if mode == "onehot":
        if restrict_to_fingerprintable:
            fingerprints = build_morgan_fingerprints()
            n_before = len(y)
            y = y[y[COL_DRUG].isin(fingerprints.keys())].reset_index(drop=True)
            print(
                f"[drug_features] onehot (restricted): {n_before - len(y)} pairs dropped "
                f"to match the fingerprint-covered population -> {len(y)} pairs remain"
            )
        drug_codes, drug_levels = pd.factorize(y[COL_DRUG], sort=True)
        n_pairs, n_drugs = len(y), len(drug_levels)
        drug_onehot = np.zeros((n_pairs, n_drugs), dtype=DTYPE)
        drug_onehot[np.arange(n_pairs), drug_codes] = 1.0
        return y, drug_onehot, n_drugs, list(drug_levels)

    if mode == "fingerprint":
        fingerprints = build_morgan_fingerprints()
        mask = y[COL_DRUG].isin(fingerprints.keys())
        n_before = len(y)
        y_filtered = y[mask].reset_index(drop=True)
        n_after = len(y_filtered)
        print(
            f"[drug_features] fingerprint mode: {n_before - n_after} pairs dropped "
            f"(drug has no resolved fingerprint) -> {n_after} pairs remain"
        )
        drug_matrix = np.stack(
            [fingerprints[d] for d in y_filtered[COL_DRUG]], axis=0
        ).astype(DTYPE)
        drug_levels = sorted(set(y_filtered[COL_DRUG]))
        return y_filtered, drug_matrix, MORGAN_FP_SIZE, drug_levels

    raise ValueError(f"Unknown drug representation mode: {mode!r}")


# --- pair-matrix construction (flat models: RF / MLP) --------------------------
def build_pair_matrix(
    X_cell: np.ndarray,
    cell_ids: pd.Index,
    y: pd.DataFrame,
    drug_mode: DrugRepMode,
    restrict_to_fingerprintable: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], pd.DataFrame]:
    """Gather cell-line features per pair and append the requested drug block.

    Returns (X, target, groups, feature_names, y_used) -- `y_used` is the row
    set that actually ended up in `X` (possibly filtered, see
    `build_drug_features`), since downstream group assertions need to
    reference the same rows.
    """
    y_used, drug_block, n_drug_features, drug_levels = build_drug_features(
        y, drug_mode, restrict_to_fingerprintable
    )

    row_of = pd.Series(np.arange(len(cell_ids)), index=cell_ids)
    rows = row_of.loc[y_used[COL_CELL_LINE]].to_numpy()

    n_pairs, n_omics = len(y_used), X_cell.shape[1]
    X = np.empty((n_pairs, n_omics + n_drug_features), dtype=DTYPE)
    X[:, :n_omics] = X_cell[rows]
    X[:, n_omics:] = drug_block

    target = y_used[COL_TARGET].to_numpy(dtype=DTYPE)
    groups = y_used[COL_CELL_LINE].to_numpy()
    names = [f"omics_{i}" for i in range(n_omics)] + [
        f"drug={d}" if drug_mode == "onehot" else f"fp_{i}"
        for i, d in enumerate(drug_levels if drug_mode == "onehot" else range(n_drug_features))
    ]

    print(
        f"[design]  {n_pairs} pairs x {X.shape[1]} features "
        f"({n_omics} omics + {n_drug_features} drug/{drug_mode}) = {X.nbytes / 1024**2:.1f} MB"
    )
    return X, target, groups, names, y_used


# --- pair tensors (dict-of-modality models: cross-attention) -------------------
def build_pair_tensors(
    omics: Dict[str, np.ndarray],
    cell_ids: pd.Index,
    y: pd.DataFrame,
    drug_mode: DrugRepMode,
    restrict_to_fingerprintable: bool = False,
) -> Tuple[Dict[str, np.ndarray], np.ndarray, np.ndarray, np.ndarray, int, pd.DataFrame]:
    """Gather per-modality cell-line rows per pair, plus the requested drug block.

    Returns (gathered_omics, drug_block, target, groups, n_drug_features, y_used).
    """
    y_used, drug_block, n_drug_features, _ = build_drug_features(
        y, drug_mode, restrict_to_fingerprintable
    )

    row_of = pd.Series(np.arange(len(cell_ids)), index=cell_ids)
    rows = row_of.loc[y_used[COL_CELL_LINE]].to_numpy()

    gathered = {key: arr[rows] for key, arr in omics.items()}
    target = y_used[COL_TARGET].to_numpy(dtype=DTYPE)
    groups = y_used[COL_CELL_LINE].to_numpy()

    omics_bytes = sum(a.nbytes for a in gathered.values())
    print(
        f"[design]  {len(y_used)} pairs  |  omics: "
        f"{' + '.join(f'{k}({v.shape[1]})' for k, v in gathered.items())}"
        f"  |  drug/{drug_mode}: {n_drug_features}  =  "
        f"{(omics_bytes + drug_block.nbytes) / 1024**2:.1f} MB"
    )
    return gathered, drug_block, target, groups, n_drug_features, y_used


# --- split ----------------------------------------------------------------------
def random_pair_split(
    n_rows: int, train_frac: float = 0.8, val_frac: float = 0.1
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Random split over (cell line, drug) PAIRS -- the MoGraphDRP protocol (their §2.5).

    Deliberately NOT leakage-free: a cell line's ~279 measurements are scattered
    across all three folds, so at test time the model has already seen ~223
    other drug responses for that exact cell line. That makes the task
    *imputation* (fill in a partially-observed row) rather than *generalization*
    (predict for an unseen cell line).

    This exists solely to produce a like-for-like number against the paper's
    published RMSE, and to test whether their XGBoost refinement gain depends on
    this leakage. `grouped_split` remains the protocol for every real result --
    see docs/09_split_protocol_comparison.md.
    """
    rng = np.random.default_rng(RANDOM_STATE)
    perm = rng.permutation(n_rows)
    n_train = int(round(train_frac * n_rows))
    n_val = int(round(val_frac * n_rows))
    train_idx = perm[:n_train]
    val_idx = perm[n_train : n_train + n_val]
    test_idx = perm[n_train + n_val :]

    for name, idx in [("train", train_idx), ("val", val_idx), ("test", test_idx)]:
        print(f"[split*]  {name:<5} {len(idx):>7} pairs  ({len(idx) / n_rows:.1%})  [RANDOM, leaky]")
    return train_idx, val_idx, test_idx


def leave_drugs_out_split(drug_groups: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """70/15/15 grouped by DRUG -- every test drug is unseen during training.

    The protocol that tests generalization to new *compounds* rather than new
    cell lines. It is the only setting in which the Morgan-fingerprint drug
    representation can demonstrate its purpose: a one-hot vector has no column
    for a drug outside the training vocabulary, so the model structurally
    cannot represent it, whereas a fingerprint is computed from structure and
    transfers to any molecule with a SMILES string.

    Mechanically identical to `grouped_split`, just keyed on `drug_id` instead
    of `sanger_model_id`. Note that a cell line WILL appear on both sides here
    (that is the point -- we are holding out drugs, not cell lines).
    """
    holdout = VAL_FRAC + TEST_FRAC
    gss1 = GroupShuffleSplit(n_splits=1, test_size=holdout, random_state=RANDOM_STATE)
    train_idx, rest_idx = next(gss1.split(np.zeros(len(drug_groups)), groups=drug_groups))

    gss2 = GroupShuffleSplit(n_splits=1, test_size=TEST_FRAC / holdout, random_state=RANDOM_STATE)
    rel_val, rel_test = next(gss2.split(np.zeros(len(rest_idx)), groups=drug_groups[rest_idx]))
    val_idx, test_idx = rest_idx[rel_val], rest_idx[rel_test]

    for name, idx in [("train", train_idx), ("val", val_idx), ("test", test_idx)]:
        print(
            f"[split-D] {name:<5} {len(idx):>7} pairs  "
            f"{len(np.unique(drug_groups[idx])):>4} drugs  ({len(idx) / len(drug_groups):.1%})"
        )

    overlap = set(drug_groups[train_idx]) & (set(drug_groups[val_idx]) | set(drug_groups[test_idx]))
    assert not overlap, f"drug leaked across splits: {sorted(overlap)[:5]}"
    return train_idx, val_idx, test_idx


def leakage_report(groups: np.ndarray, train_idx: np.ndarray, test_idx: np.ndarray) -> Dict[str, float]:
    """Quantify how much cell-line information a split leaks from train into test."""
    train_cells = set(groups[train_idx])
    test_cells = set(groups[test_idx])
    shared = train_cells & test_cells
    leaked_rows = int(np.isin(groups[test_idx], list(shared)).sum()) if shared else 0

    train_counts = pd.Series(groups[train_idx]).value_counts()
    seen_per_test_row = (
        float(np.mean([train_counts.get(c, 0) for c in groups[test_idx]])) if len(test_idx) else 0.0
    )
    return {
        "test_cell_lines": len(test_cells),
        "also_in_train": len(shared),
        "test_rows_with_seen_cell_line_pct": 100.0 * leaked_rows / max(len(test_idx), 1),
        "mean_train_rows_per_test_cell_line": seen_per_test_row,
    }


def grouped_split(groups: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """70/15/15 by cell line, `random_state=RANDOM_STATE` pinned for every caller."""
    holdout = VAL_FRAC + TEST_FRAC
    gss1 = GroupShuffleSplit(n_splits=1, test_size=holdout, random_state=RANDOM_STATE)
    train_idx, rest_idx = next(gss1.split(np.zeros(len(groups)), groups=groups))

    gss2 = GroupShuffleSplit(n_splits=1, test_size=TEST_FRAC / holdout, random_state=RANDOM_STATE)
    rel_val, rel_test = next(gss2.split(np.zeros(len(rest_idx)), groups=groups[rest_idx]))
    val_idx, test_idx = rest_idx[rel_val], rest_idx[rel_test]

    for name, idx in [("train", train_idx), ("val", val_idx), ("test", test_idx)]:
        print(
            f"[split]   {name:<5} {len(idx):>7} pairs  "
            f"{len(np.unique(groups[idx])):>5} cell lines  ({len(idx) / len(groups):.1%})"
        )

    overlap = set(groups[train_idx]) & (set(groups[val_idx]) | set(groups[test_idx]))
    assert not overlap, f"cell line leaked across splits: {sorted(overlap)[:5]}"
    return train_idx, val_idx, test_idx


# --- evaluation -------------------------------------------------------------------
def evaluate(y_true: np.ndarray, y_pred: np.ndarray, threshold: float) -> Dict[str, float]:
    """Full metric set: RMSE, MAE, R^2, PCC, SCC (regression) + AUC, F1 (binarized)."""
    err = y_true - y_pred
    rmse = float(np.sqrt(np.mean(err**2)))
    mae = float(np.mean(np.abs(err)))
    ss_res = float(np.sum(err**2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    pcc = float(np.corrcoef(y_true, y_pred)[0, 1]) if np.std(y_pred) > 0 else float("nan")
    scc = float(spearmanr(y_true, y_pred).statistic) if np.std(y_pred) > 0 else float("nan")

    y_true_bin = (y_true >= threshold).astype(int)
    y_pred_bin = (y_pred >= threshold).astype(int)
    auc = (
        float(roc_auc_score(y_true_bin, y_pred))
        if len(np.unique(y_true_bin)) > 1
        else float("nan")
    )
    f1 = float(f1_score(y_true_bin, y_pred_bin, zero_division=0))

    return {"rmse": rmse, "mae": mae, "r2": r2, "pcc": pcc, "scc": scc, "auc": auc, "f1": f1}


def mean_only_floor(y_train: np.ndarray, y_test: np.ndarray) -> float:
    """RMSE of always predicting the training-fold mean -- the floor every model must beat."""
    return float(np.sqrt(np.mean((y_test - y_train.mean()) ** 2)))


def print_metric_block(title: str, val: Dict[str, float], test: Dict[str, float], floor: float) -> None:
    """Standard formatted console report, matching the existing baseline scripts' style."""
    cols = ["rmse", "mae", "r2", "pcc", "scc", "auc", "f1"]
    header = "".join(f"{c.upper():>10}" for c in cols)
    print("\n" + "=" * (12 + 10 * len(cols)))
    print(title)
    print("=" * (12 + 10 * len(cols)))
    print(f"{'':<12}{header}")
    print(f"{'Validation':<12}" + "".join(f"{val[c]:>10.4f}" for c in cols))
    print(f"{'Test':<12}" + "".join(f"{test[c]:>10.4f}" for c in cols))
    print(f"{'Mean-only':<12}{floor:>10.4f}" + "".join(f"{'--':>10}" for _ in cols[1:]))
    print("-" * (12 + 10 * len(cols)))
