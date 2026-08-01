"""Preprocess raw multi-omics matrices into a uniform PCA-reduced representation.

Each modality (transcriptomics, genomics, proteomics) arrives as a wide
DataFrame indexed by cell line ID with one column per feature. This module
aligns modalities to a shared set of cell lines, applies modality-specific
normalization, and compresses every modality to the same target dimension
`d_target` via PCA so they can be fed into the cross-attention fusion module
in `src/models/cross_attention_fusion.py`.

Fitted scalers/imputers/PCA models are saved to disk so the exact training-time
transform can be replayed at inference without re-fitting on new data.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import joblib
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler

GE_KEY = "GE"
MUT_CNV_KEY = "Mut_CNV"
PROTEOMICS_KEY = "Proteomics"
MODALITIES = (GE_KEY, MUT_CNV_KEY, PROTEOMICS_KEY)


class OmicsPreprocessingPipeline:
    """Aligns, normalizes, and PCA-compresses GE / Mut_CNV / Proteomics matrices."""

    def __init__(self, d_target: int = 128, random_state: int = 42, model_dir: str | Path = "models/omics_pca"):
        self.d_target = d_target
        self.random_state = random_state
        self.model_dir = Path(model_dir)
        self.pipelines_: Dict[str, Pipeline] = {}

    def _build_pipeline(self, modality: str, n_samples: int, n_features: int) -> Pipeline:
        """Construct the modality-specific normalization + PCA pipeline.

        n_components is capped at min(d_target, n_samples, n_features) since PCA
        cannot produce more components than the smaller of the sample or feature count.
        """
        n_components = min(self.d_target, n_samples, n_features)
        if n_components < self.d_target:
            print(
                f"[OmicsPreprocessingPipeline] WARNING: {modality} capped PCA components "
                f"to {n_components} (requested d_target={self.d_target}, "
                f"n_samples={n_samples}, n_features={n_features})."
            )

        steps: List[tuple] = []
        if modality == PROTEOMICS_KEY:
            steps.append(("impute_median", SimpleImputer(strategy="median")))
        if modality == GE_KEY:
            # Variance stabilization for raw-count-scale transcriptomics. Skip this
            # step if the source column is already log-transformed upstream (e.g.
            # CCLE's "...TPMLogp1" file) to avoid double log1p-ing the data.
            steps.append(("log1p", FunctionTransformer(np.log1p, validate=True)))
        steps.append(("zscore", StandardScaler()))
        steps.append(("pca", PCA(n_components=n_components, random_state=self.random_state)))
        return Pipeline(steps)

    def align_by_cell_line(self, omics: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        """Intersect and reorder modality DataFrames onto a shared, sorted cell line ID index."""
        common_ids = None
        for df in omics.values():
            ids = set(df.index)
            common_ids = ids if common_ids is None else (common_ids & ids)
        common_ids = sorted(common_ids or set())
        if not common_ids:
            raise ValueError("No overlapping cell line IDs found across the provided omics modalities.")
        return {name: df.loc[common_ids] for name, df in omics.items()}

    def fit_transform(self, omics: Dict[str, pd.DataFrame]) -> Dict[str, np.ndarray]:
        """Align modalities, fit a pipeline per modality, and return PCA-compressed arrays."""
        aligned = self.align_by_cell_line(omics)
        transformed: Dict[str, np.ndarray] = {}
        for modality, df in aligned.items():
            pipeline = self._build_pipeline(modality, n_samples=df.shape[0], n_features=df.shape[1])
            transformed[modality] = pipeline.fit_transform(df.values).astype(np.float32)
            self.pipelines_[modality] = pipeline
        return transformed

    def transform(self, omics: Dict[str, pd.DataFrame]) -> Dict[str, np.ndarray]:
        """Apply already-fitted pipelines (inference phase) to new, aligned omics data."""
        if not self.pipelines_:
            raise RuntimeError("No fitted pipelines found. Call fit_transform() or load() first.")
        aligned = self.align_by_cell_line(omics)
        return {
            modality: self.pipelines_[modality].transform(df.values).astype(np.float32)
            for modality, df in aligned.items()
        }

    def save(self) -> None:
        """Persist fitted pipelines to `self.model_dir` for reuse at inference time."""
        self.model_dir.mkdir(parents=True, exist_ok=True)
        for modality, pipeline in self.pipelines_.items():
            joblib.dump(pipeline, self.model_dir / f"{modality}_pipeline.joblib")

    def load(self) -> None:
        """Load previously fitted pipelines from `self.model_dir`."""
        self.pipelines_ = {
            modality: joblib.load(self.model_dir / f"{modality}_pipeline.joblib")
            for modality in MODALITIES
            if (self.model_dir / f"{modality}_pipeline.joblib").exists()
        }
