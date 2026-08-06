"""Turn ingest_and_align.py's long-format aligned CSVs into wide omics matrices.

`ingest_and_align.py` produces one row per (cell line, gene) for RNA-seq/CNV/
mutations, joined to a shared Sanger `standard_model_id`. `OmicsPreprocessingPipeline`
(src/data/omics_preprocessing.py) instead expects one row per cell line and one
column per feature. This module bridges the two:

- GE / CNV are pivoted directly (cell line x gene, values averaged over any
  duplicate rows for the same gene).
- Mutations are pivoted into a binary presence/absence matrix per gene, since a
  cell line can carry several distinct variants in the same gene and
  `OmicsPreprocessingPipeline` needs one numeric value per (cell line, gene).
- Mutation and CNV wide matrices are combined column-wise into a single
  Mut_CNV matrix, matching the proposal report's "mutations and CNVs combined
  into a feature matrix" preprocessing step.
- Proteomics is loaded directly from the raw CCLE/Sanger proteomics archive,
  which already ships as a wide matrix (rows=cell line, columns=protein), so
  no pivot is needed there, just a bespoke header parse.

Output matrices are written to `data/processed/wide/` and are ready to be
passed straight into `OmicsPreprocessingPipeline.fit_transform()`.
"""

from __future__ import annotations

import argparse
import zipfile
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

ALIGNED_ROOT = Path("data/processed/aligned")
WIDE_OUTPUT_ROOT = Path("data/processed/wide")
PROTEOMICS_ZIP_PATH = Path("data/raw/auxiliary/Proteomics_20250211.zip")
PROTEOMICS_MEMBER = "Protein_matrix_averaged_20250211.tsv"

INDEX_COL = "standard_model_id"


def pivot_long_to_wide(path: Path, columns_col: str, value_col: str) -> pd.DataFrame:
    """Pivot a long aligned CSV into a (cell line x feature) matrix.

    Duplicate (cell line, feature) rows are averaged rather than erroring, since
    the source data carries an explicit `duplicate` flag for exactly this case.
    """
    df = pd.read_csv(
        path,
        usecols=[INDEX_COL, columns_col, value_col],
        dtype={INDEX_COL: "category", columns_col: "category"},
    )
    df = df[(df[INDEX_COL].astype(str) != "") & (df[columns_col].astype(str) != "")]
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
    wide = df.pivot_table(index=INDEX_COL, columns=columns_col, values=value_col, aggfunc="mean")
    return wide.astype(np.float32)


def pivot_mutation_presence(path: Path, columns_col: str = "gene_symbol") -> pd.DataFrame:
    """Pivot mutation call rows into a binary (cell line x gene) presence matrix."""
    df = pd.read_csv(path, usecols=[INDEX_COL, columns_col], dtype={INDEX_COL: "category", columns_col: "category"})
    df = df[(df[INDEX_COL].astype(str) != "") & (df[columns_col].astype(str) != "")]
    df["present"] = np.float32(1.0)
    wide = df.pivot_table(index=INDEX_COL, columns=columns_col, values="present", aggfunc="max", fill_value=0.0)
    return wide.astype(np.float32)


def combine_mutation_and_cnv(mutation_wide: pd.DataFrame, cnv_wide: pd.DataFrame) -> pd.DataFrame:
    """Column-wise union of the binary mutation matrix and the continuous CNV matrix."""
    combined = mutation_wide.add_prefix("MUT_").join(cnv_wide.add_prefix("CNV_"), how="outer")
    return combined.fillna(0.0).astype(np.float32)


def load_proteomics_wide(zip_path: Path = PROTEOMICS_ZIP_PATH, member: str = PROTEOMICS_MEMBER) -> pd.DataFrame:
    """Load the CCLE/Sanger proteomics matrix, which ships already wide.

    The TSV has a 3-row header block: row 0 holds uniprot_id column labels, row 1
    holds gene symbol column labels, row 2 labels the row-key columns
    ("model_name", "model_id"). Data rows start at row 3 as
    [model_name, model_id, <protein intensity values...>].
    """
    with zipfile.ZipFile(zip_path) as archive, archive.open(member) as handle:
        raw = pd.read_csv(handle, sep="\t", header=None, low_memory=False)

    uniprot_ids = raw.iloc[0, 2:].tolist()
    data = raw.iloc[3:].copy()
    data.columns = ["model_name", "model_id"] + uniprot_ids
    data = data.rename(columns={"model_id": INDEX_COL}).set_index(INDEX_COL)
    wide = data[uniprot_ids].apply(pd.to_numeric, errors="coerce")
    return wide.astype(np.float32)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Build and write the GE / Mut_CNV / Proteomics wide matrices."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligned-root", type=Path, default=ALIGNED_ROOT, help="Output directory from ingest_and_align.py.")
    parser.add_argument("--proteomics-zip", type=Path, default=PROTEOMICS_ZIP_PATH, help="Path to the raw proteomics archive.")
    parser.add_argument("--output-dir", type=Path, default=WIDE_OUTPUT_ROOT, help="Directory for the wide matrix outputs.")
    args = parser.parse_args(argv)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    ge_wide = pivot_long_to_wide(args.aligned_root / "rnaseq_aligned.csv", "gene_symbol", "rsem_tpm")
    cnv_wide = pivot_long_to_wide(args.aligned_root / "cnv_aligned.csv", "symbol", "total_copy_number")
    mutation_wide = pivot_mutation_presence(args.aligned_root / "mutations_aligned.csv")
    mut_cnv_wide = combine_mutation_and_cnv(mutation_wide, cnv_wide)
    proteomics_wide = load_proteomics_wide(args.proteomics_zip)

    ge_path = args.output_dir / "GE_wide.csv"
    mut_cnv_path = args.output_dir / "Mut_CNV_wide.csv"
    proteomics_path = args.output_dir / "Proteomics_wide.csv"

    ge_wide.to_csv(ge_path)
    mut_cnv_wide.to_csv(mut_cnv_path)
    proteomics_wide.to_csv(proteomics_path)

    print(f"GE_wide: {ge_wide.shape} -> {ge_path}")
    print(f"Mut_CNV_wide: {mut_cnv_wide.shape} (mutation cols: {mutation_wide.shape[1]}, cnv cols: {cnv_wide.shape[1]}) -> {mut_cnv_path}")
    print(f"Proteomics_wide: {proteomics_wide.shape} -> {proteomics_path}")

    common_ids = set(ge_wide.index) & set(mut_cnv_wide.index) & set(proteomics_wide.index)
    print(f"Cell lines common to all 3 modalities: {len(common_ids)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
