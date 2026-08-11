"""
Early fusion (concatenation) of PCA-compressed multi-omics matrices.

Assumes each CSV is z-scored + PCA-compressed, with cell line IDs in column 0
and rows already aligned across all three modalities.

Memory strategy:
  1. Read only the header of each file to get block widths (cheap).
  2. Pre-allocate the final float32 array once.
  3. Stream each modality in, copy its block, free it immediately.
Peak RSS ~= final array + one modality block (no 2x hstack copy).
"""

from __future__ import annotations

import gc
from pathlib import Path

import numpy as np
import pandas as pd

# --- config ---------------------------------------------------------------
DATA_DIR = Path("data/processed")

MODALITIES = {
    "genomics": DATA_DIR / "genomics_pca.csv",
    "transcriptomics": DATA_DIR / "transcriptomics_pca.csv",
    "proteomics": DATA_DIR / "proteomics_pca.csv",
}

DTYPE = np.float32  # half the RAM of float64; PCA scores don't need f64 precision
# --------------------------------------------------------------------------


def peek_shape(path: Path) -> tuple[int, list[str]]:
    """Return (n_feature_cols, column_names) by reading the header only."""
    header = pd.read_csv(path, index_col=0, nrows=0)
    return header.shape[1], list(header.columns)


def load_block(path: Path, n_cols: int) -> pd.DataFrame:
    """Load one modality directly as float32."""
    return pd.read_csv(
        path,
        index_col=0,
        dtype={i: DTYPE for i in range(1, n_cols + 1)},
        engine="c",
    )


def fuse(modalities: dict[str, Path]) -> tuple[np.ndarray, pd.Index, list[str]]:
    widths, col_names = {}, {}
    for name, path in modalities.items():
        widths[name], col_names[name] = peek_shape(path)

    total_cols = sum(widths.values())

    index: pd.Index | None = None
    fused: np.ndarray | None = None
    feature_names: list[str] = []
    offset = 0

    for name, path in modalities.items():
        block = load_block(path, widths[name])

        if index is None:
            index = block.index
            fused = np.empty((len(index), total_cols), dtype=DTYPE)
            est_gb = fused.nbytes / 1024**3
            print(f"[alloc] fused array {fused.shape} -> {est_gb:.3f} GB ({DTYPE.__name__})")
        else:
            # rows assumed pre-aligned; fail loudly if that assumption breaks
            if len(block.index) != len(index) or not block.index.equals(index):
                raise ValueError(
                    f"Row misalignment in '{name}': expected {len(index)} cell lines "
                    f"matching the first modality, got {len(block.index)}."
                )

        fused[:, offset : offset + widths[name]] = block.to_numpy(dtype=DTYPE, copy=False)
        feature_names.extend(f"{name}_{c}" for c in col_names[name])
        offset += widths[name]

        print(f"[load]  {name:<16} {block.shape[0]} x {block.shape[1]}  "
              f"-> cols [{offset - widths[name]}:{offset}]")

        del block
        gc.collect()

    return fused, index, feature_names


def main() -> None:
    X, cell_lines, feature_names = fuse(MODALITIES)

    n_bytes = X.nbytes
    print("\n=== Early fusion complete ===")
    print(f"Final shape          : {X.shape[0]} cell lines x {X.shape[1]} features")
    print(f"dtype                : {X.dtype}")
    print(f"Feature-matrix memory: {n_bytes / 1024**2:.2f} MB  ({n_bytes / 1024**3:.4f} GB)")
    print(f"Per-sample footprint : {n_bytes / X.shape[0] / 1024:.2f} KB")
    print(f"Cell lines           : {cell_lines[0]} ... {cell_lines[-1]}")
    print(f"First/last feature   : {feature_names[0]} / {feature_names[-1]}")

    budget_gb = 4.0
    print(f"Within {budget_gb:.0f} GB budget: {n_bytes / 1024**3 < budget_gb}")

    # Optional: cache as .npy so downstream runs skip CSV parsing entirely.
    # np.save(DATA_DIR / "fused_early.npy", X)
    # pd.Series(feature_names).to_csv(DATA_DIR / "fused_feature_names.csv", index=False, header=False)
    # cell_lines.to_series().to_csv(DATA_DIR / "fused_cell_lines.csv", index=False, header=False)


if __name__ == "__main__":
    main()