"""Per-drug standardization of the ln(IC50) regression target.

Every model so far regresses raw ln(IC50) directly, so most of the loss is
spent learning *which compound this is* rather than *how this cell line
responds to it*: potency varies far more between drugs than between cell
lines within a drug. Standardizing the target per drug removes that offset,
leaving the model to predict the part that actually depends on the omics
profile.

The transform is:

    z = (y - mean_d) / std_d      for drug d

with the inverse applied to predictions before scoring, so RMSE/R2/PCC stay
in ln(IC50) units and remain directly comparable to every row in
docs/results.md.

**No leakage by construction.** `fit` accepts only the training rows -- there
is no argument through which test labels could reach it. Under the
cell-line-grouped split every drug appears in training, so this is legitimate:
the statistics describe compounds the model has already seen, not the held-out
cell lines it is being asked to predict.

**Not usable for leave-drugs-out.** There, held-out compounds have no training
rows and therefore no statistics, so the fallback would apply to every test
pair and the transform would degenerate to a global z-score. Keep this as a
separate arm of the cell-line-grouped protocol rather than a global default.
"""

from __future__ import annotations

from typing import Dict, Hashable, Tuple

import numpy as np

# Below this many training rows a drug's own mean/std are too unstable to
# trust, so the global training statistics are used instead.
MIN_TRAIN_SAMPLES = 5
# Absolute floor on a drug's std, so an exactly-constant drug cannot divide by
# zero.
MIN_STD = 1e-3
# Relative floor, as a fraction of the global training std. This matters more
# than the absolute one: a drug whose training responses happen to be nearly
# identical would otherwise divide by a tiny number, amplifying any deviation
# in its held-out rows into a huge z and letting that one compound dominate
# the loss. Flooring at a fraction of the global spread caps that blow-up.
MIN_STD_RATIO = 0.1


class PerDrugTargetScaler:
    """Standardize ln(IC50) within each drug, using training rows only.

    Typical use, after the split::

        scaler = PerDrugTargetScaler().fit(y[train_idx], drugs[train_idx])
        z = scaler.transform(y, drugs)          # train the model on z
        preds = scaler.inverse_transform(z_hat, drugs)   # score in ln(IC50)
    """

    def __init__(self, min_samples: int = MIN_TRAIN_SAMPLES, min_std: float = MIN_STD,
                 min_std_ratio: float = MIN_STD_RATIO):
        self.min_samples = min_samples
        self.min_std = min_std
        self.min_std_ratio = min_std_ratio
        self.n_floored_groups_: int = 0
        self.means_: Dict[Hashable, float] = {}
        self.stds_: Dict[Hashable, float] = {}
        self.global_mean_: float = 0.0
        self.global_std_: float = 1.0
        self.n_fallback_groups_: int = 0

    def fit(self, y_train: np.ndarray, groups_train: np.ndarray) -> "PerDrugTargetScaler":
        """Compute per-drug statistics from **training rows only**."""
        if len(y_train) != len(groups_train):
            raise ValueError(
                f"y_train ({len(y_train)}) and groups_train ({len(groups_train)}) "
                f"must be the same length"
            )
        y_train = np.asarray(y_train, dtype=np.float64)

        self.global_mean_ = float(y_train.mean())
        self.global_std_ = max(float(y_train.std()), self.min_std)
        std_floor = max(self.min_std, self.min_std_ratio * self.global_std_)

        order = np.argsort(groups_train, kind="stable")
        sorted_groups = np.asarray(groups_train)[order]
        sorted_y = y_train[order]
        boundaries = np.flatnonzero(np.r_[True, sorted_groups[1:] != sorted_groups[:-1]])
        splits = np.split(sorted_y, boundaries[1:])

        self.means_, self.stds_ = {}, {}
        self.n_fallback_groups_ = 0
        self.n_floored_groups_ = 0
        for group, values in zip(sorted_groups[boundaries], splits):
            if len(values) < self.min_samples:
                self.n_fallback_groups_ += 1
                continue
            raw_std = float(values.std())
            if raw_std < std_floor:
                self.n_floored_groups_ += 1
            self.means_[group] = float(values.mean())
            self.stds_[group] = max(raw_std, std_floor)

        print(
            f"[target_scaling] fitted on {len(y_train)} training rows: "
            f"{len(self.means_)} drugs with own statistics, "
            f"{self.n_fallback_groups_} below min_samples={self.min_samples} "
            f"-> global (mean={self.global_mean_:.4f}, std={self.global_std_:.4f}); "
            f"{self.n_floored_groups_} had std floored at {std_floor:.4f}"
        )
        return self

    def params_for(self, groups: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Per-row (mean, std), falling back to global stats for unknown drugs."""
        if not self.means_ and self.global_std_ == 1.0 and self.global_mean_ == 0.0:
            raise RuntimeError("PerDrugTargetScaler.fit must be called before use")
        means = np.array(
            [self.means_.get(g, self.global_mean_) for g in groups], dtype=np.float64
        )
        stds = np.array(
            [self.stds_.get(g, self.global_std_) for g in groups], dtype=np.float64
        )
        return means, stds

    def transform(self, y: np.ndarray, groups: np.ndarray) -> np.ndarray:
        """ln(IC50) -> per-drug z-score."""
        means, stds = self.params_for(groups)
        return ((np.asarray(y, dtype=np.float64) - means) / stds).astype(np.float32)

    def inverse_transform(self, z: np.ndarray, groups: np.ndarray) -> np.ndarray:
        """Per-drug z-score -> ln(IC50), for scoring in the original units."""
        means, stds = self.params_for(groups)
        return (np.asarray(z, dtype=np.float64) * stds + means).astype(np.float32)


if __name__ == "__main__":
    rng = np.random.default_rng(42)

    # 40 drugs with deliberately different potency offsets and spreads, plus
    # one rare drug that must fall back and one with no variance at all.
    n_drugs, per_drug = 40, 60
    offsets = rng.normal(3.0, 2.5, size=n_drugs)
    spreads = rng.uniform(0.3, 2.0, size=n_drugs)
    groups = np.repeat(np.arange(n_drugs), per_drug)
    y = rng.normal(offsets[groups], spreads[groups]).astype(np.float32)

    groups = np.r_[groups, [98, 98, 99, 99, 99, 99, 99, 99, 99, 99]]
    y = np.r_[y, rng.normal(1.0, 1.0, 2).astype(np.float32), np.full(8, 7.5, np.float32)]

    idx = rng.permutation(len(y))
    train_idx, test_idx = idx[: int(0.7 * len(idx))], idx[int(0.7 * len(idx)) :]

    scaler = PerDrugTargetScaler().fit(y[train_idx], groups[train_idx])

    # Round trip must be lossless (to float32 precision).
    z = scaler.transform(y, groups)
    back = scaler.inverse_transform(z, groups)
    print(f"\nround-trip max abs error: {np.abs(back - y).max():.2e}")
    assert np.allclose(back, y, atol=1e-4)

    # Standardized training targets should be ~N(0, 1) per drug.
    z_train = scaler.transform(y[train_idx], groups[train_idx])
    common = [g for g in np.unique(groups[train_idx]) if g in scaler.means_]
    per_drug_means = [z_train[groups[train_idx] == g].mean() for g in common]
    varying = [g for g in common if y[train_idx][groups[train_idx] == g].std() > 1e-6]
    per_drug_stds = [z_train[groups[train_idx] == g].std() for g in varying]
    print(f"per-drug standardized mean: {np.mean(per_drug_means):+.4f} "
          f"(max |mean| {np.max(np.abs(per_drug_means)):.4f})")
    print(f"per-drug standardized std:  {np.mean(per_drug_stds):.4f}")
    assert np.max(np.abs(per_drug_means)) < 0.2
    assert np.allclose(per_drug_stds, 1.0, atol=0.2)

    # Statistics must come from training rows only.
    for g in list(scaler.means_)[:5]:
        mask = (groups[train_idx] == g)
        assert np.isclose(scaler.means_[g], y[train_idx][mask].mean(), atol=1e-4)
    print("per-drug statistics match the training rows exactly")

    # The rare drug (2 rows) and an entirely unseen drug both fall back.
    assert 98 not in scaler.means_, "drug with 2 rows should fall back"
    m, s = scaler.params_for(np.array([98, 12345]))
    assert np.allclose(m, scaler.global_mean_) and np.allclose(s, scaler.global_std_)
    print("rare and unseen drugs fall back to global statistics")

    # Zero-variance drug must not produce inf/nan.
    z99 = scaler.transform(y[groups == 99], groups[groups == 99])
    assert np.isfinite(z99).all(), "zero-variance drug produced non-finite targets"
    print("zero-variance drug handled without inf/nan")

    print("\nAll per-drug target scaling checks passed.")
