# Seed Variance — Which Results Survive Repetition

**Script:** [`experiments/13_seed_variance/seed_variance.py`](../experiments/13_seed_variance/seed_variance.py)
**Raw results:** `experiments/13_seed_variance/seed_variance_results.csv`, `seed_variance_summary.csv`
**Companions:** [`results.md`](./results.md), [`11_molecular_graph_results.md`](./11_molecular_graph_results.md), [`12_full_architecture_results.md`](./12_full_architecture_results.md)
**Date:** 2026-09-22 · **Runtime:** 15 runs

## What this covers

Every row in [`results.md`](./results.md) is a single run at `TORCH_SEED = 42`.
That is adequate for ranking configurations separated by a lot, but the
differences this project turned on are small — 0.018 between the molecular
graph and fingerprint arms, 0.038 between the full architecture and the same
model without the PPI graph.

This re-runs the three configurations those claims depend on across five seeds
(42–46) and reports mean ± std.

**The split is deliberately not reseeded.** `grouped_split` keeps
`RANDOM_STATE = 42`, so every run sees the identical partition and the only
things varying are weight initialization, batch shuffling and CUDA
nondeterminism. Because the split is shared, comparisons are **paired**: at a
fixed seed, two runs differ only in architecture, so the spread of the per-seed
*difference* is a tighter test than whether two independent means separate.

The original runners are reused via `load_module` with `TORCH_SEED` overridden,
so a swept run is the same code path as the recorded matrix run.

## Results

Test RMSE per seed:

| Seed | E04 (MLP, Proteomics, fp) | E10 (CrossAttn, GE+Prot, fp) | MG (CrossAttn, GE+Prot, mol. graph) |
|---|---|---|---|
| 42 | 1.3155 | 1.2997 | 1.3044 |
| 43 | 1.3257 | 1.3087 | 1.3049 |
| 44 | 1.3630 | 1.3203 | 1.3262 |
| 45 | 1.2878 | 1.2999 | 1.3208 |
| 46 | 1.3446 | 1.2861 | 1.3408 |

Summary:

| Config | Mean RMSE | Std | Min | Max | Recorded (single seed) | Params |
|---|---|---|---|---|---|---|
| **E10** fingerprint | **1.3029** | **0.0126** | 1.2861 | 1.3203 | 1.3205 | 742,017 |
| MG molecular graph | 1.3194 | 0.0153 | 1.3044 | 1.3408 | 1.3021 | 501,377 |
| E04 MLP proteomics | 1.3273 | 0.0286 | 1.2878 | 1.3630 | 1.2843 | 1,280,769 |

Paired per-seed differences (negative = first config better):

| Comparison | Mean Δ | Std Δ | Seeds won | Separable |
|---|---|---|---|---|
| MG − E10 | +0.0165 | 0.0231 | 1 / 5 | No |
| MG − E04 | −0.0079 | 0.0260 | 4 / 5 | No |
| E10 − E04 | −0.0244 | 0.0272 | 4 / 5 | No |

`Separable` flags whether the mean difference exceeds its own spread. It is a
readability aid, **not** a significance test at n=5.

## 1. Two previously reported findings do not survive

**Molecular graphs do not beat fingerprints.** The single-seed comparison
(MG 1.3021 vs E10 1.3205) suggested a 1.4% improvement. Across five seeds the
molecular graph is **0.0165 worse** on average and loses 4 of 5 seeds. The
original gap was noise, and it pointed in the wrong direction.

**Proteomics-alone is not the best configuration.** E04's recorded 1.2843 was
the favourable end of a wide distribution: its mean is **1.3273, the worst of
the three**, with more than double E10's spread. Its `best_epoch` across seeds
(6, 10, 3, 9, 6) shows it converges and begins overfitting within ten epochs,
consistent with being the least stable of the three.

## 2. What is actually true

No pair of configurations is separable at n=5 — every `std Δ` exceeds its
`mean Δ`. The defensible statement is that all three perform **equivalently**
on the cell-line-grouped protocol, at roughly **RMSE 1.30–1.33**.

They differ in stability and cost rather than accuracy. E10 has both the best
mean and the tightest spread; the molecular graph matches it at 32% fewer
parameters; E04 is the largest model, the least stable, and the fastest to
overfit.

## 3. The practical noise threshold is ~0.03, not ~0.01

Three sources of variance, in increasing size:

| Source | Magnitude | Evidence |
|---|---|---|
| CUDA nondeterminism, same seed | ~0.008 | Two runs of one config at seed 42: 1.3021, 1.3099 |
| Seed-to-seed, same session | ~0.013–0.029 | The std column above |
| Cross-session (different hardware/build) | ~0.021 | E10 recorded 1.3205 vs 1.2997 at seed 42 here |

E04 spans 1.2878–1.3630 across five seeds — a range of **0.075**.

**Any difference below ~0.03 RMSE between single runs should be treated as
undetectable.** This applies retroactively to every row in
[`results.md`](./results.md): the ranking of closely-spaced rows is not
reliable.

## 4. Status of the project's claims under this threshold

| Claim | Gap | Status |
|---|---|---|
| GAT fails on this data | ~1.0 | **Solid** — 3 independent runs plus the benchmark's Table 3 |
| `n_tokens` has no effect | 0.009 | **Confirmed null** — spread below the same-seed noise floor |
| PPI graph costs RMSE | 0.038 | **Borderline** — ~2.5× seed spread, not repeated |
| Molecular graph beats fingerprint | 0.018 | **Refuted** — reverses across seeds |
| Proteomics-alone is best | 0.018 | **Refuted** — worst mean of the three |

## Caveats

- **The split is never varied.** This measures training variance only. Which
  80 cell lines land in the test fold is a single draw at `random_state=42`,
  and split variance is plausibly larger than anything measured here. A
  grouped k-fold cross-validation would quantify it and remains unrun.
- n=5 is small. The paired differences have wide intervals and none reach
  separability; a larger sweep could resolve E10 vs E04, whose paired
  difference (−0.0244 ± 0.0272) is the closest to separating.
- Only three configurations were swept. The full architecture
  ([`12_full_architecture_results.md`](./12_full_architecture_results.md)) and
  the GAT variants were not, so their single-seed numbers carry the same
  uncertainty the three configurations above were found to have.
