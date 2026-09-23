"""build_results_table.py — consolidate every matrix run into docs/results.md.

Phase 5 of the experiment matrix (docs/plan/experiment_matrix_plan.md). Reads
the per-family result CSVs and emits one master markdown table, so the summary
can never drift from the numbers the experiments actually produced. The emitted
table is also the artifact to paste into Notion (it imports markdown tables
directly).

Run from the repository root:
    python experiments/build_results_table.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# The original matrix. Ranked on its own so the E01-E59 identifiers these
# docs already cross-reference stay stable as later experiments are added.
LEGACY_SOURCES = [
    ("experiments/06_full_matrix/rf_matrix_results.csv", "06_rf_ablation_results.md"),
    ("experiments/06_full_matrix/mlp_matrix_results.csv", "06_mlp_ablation_results.md"),
    ("experiments/06_full_matrix/cross_attention_matrix_results.csv", "06_cross_attention_ablation_results.md"),
    ("experiments/07_gnn_ablation/gnn_results.csv", "07_gnn_ablation_results.md"),
    ("experiments/07_gnn_ablation/hetero_ic50_gnn_results.csv", "hetero_gnn_test_bugfixes.md"),
]
# Appended after the legacy block, so these take identifiers from E60 onward.
LATER_SOURCES = [
    ("experiments/11_molecular_graph/molecular_graph_matrix_results.csv", "11_molecular_graph_results.md"),
    ("experiments/12_full_architecture/full_architecture_results.csv", "12_full_architecture_results.md"),
]
SEED_VARIANCE_CSV = Path("experiments/13_seed_variance/seed_variance_summary.csv")

# Drug-identity arms are retained as a reference section, never ranked as
# headline results: one-hot cannot generalize to unseen compounds (R^2 0.024 in
# 10_leave_drugs_out_results.md), so it is not a representation this project
# builds on. See the "Drug-identity reference arms" section below.
REFERENCE_DRUG_REPS = {"onehot", "onehot_restricted"}

# Below this, a difference between single runs is not distinguishable from
# run-to-run variance -- see 13_seed_variance_results.md.
NOISE_THRESHOLD = 0.03
# Predicting each drug's training mean, with no omics input at all, on the
# 111,799-pair population (computed in 11_molecular_graph_results.md).
PER_DRUG_MEAN_RMSE = 1.4889
ENSEMBLE_CSV = Path("experiments/08_ensemble_refinement/ensemble_results.csv")
OUTPUT = Path("docs/results.md")

METRICS = ["test_rmse", "test_mae", "test_r2", "test_pcc", "test_scc", "test_auc", "test_f1"]


def _load(sources) -> pd.DataFrame:
    frames = []
    for path, doc in sources:
        if not Path(path).exists():
            print(f"[skip] {path} not found")
            continue
        df = pd.read_csv(path)
        df["doc"] = doc
        if "mutation_edges" in df.columns:
            df["notes"] = df["mutation_edges"].map(
                {True: "+mutation edges", False: "no mutation edges"}
            )
        else:
            df["notes"] = ""
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True).sort_values("test_rmse").reset_index(drop=True)


def load_all() -> pd.DataFrame:
    """Legacy runs keep E01-E59; later experiments are appended from E60."""
    legacy = _load(LEGACY_SOURCES)
    later = _load(LATER_SOURCES)
    legacy["exp_id"] = [f"E{i:02d}" for i in range(1, len(legacy) + 1)]
    if not later.empty:
        later["exp_id"] = [f"E{i:02d}" for i in range(len(legacy) + 1, len(legacy) + len(later) + 1)]
    return pd.concat([legacy, later], ignore_index=True)


def fmt_row(r: pd.Series) -> str:
    cells = [
        r["exp_id"],
        str(r["model"]),
        str(r["omics"]),
        str(r["drug_rep"]),
        str(r["notes"]) or "—",
        f"{int(r['n_pairs']):,}",
        f"{r['test_rmse']:.4f}",
        f"{r['test_mae']:.4f}",
        f"{r['test_r2']:.4f}",
        f"{r['test_pcc']:.4f}",
        f"{r['test_scc']:.4f}",
        f"{r['test_auc']:.4f}",
        f"{r['test_f1']:.4f}",
        f"{int(r['params']):,}" if pd.notna(r.get("params")) else "—",
        f"{r['fit_seconds']:.1f}" if pd.notna(r.get("fit_seconds")) else "—",
        f"[{r['doc'].replace('.md', '')}](./{r['doc']})",
    ]
    return "| " + " | ".join(cells) + " |"


def main() -> None:
    df = load_all()
    is_reference = df["drug_rep"].isin(REFERENCE_DRUG_REPS)
    primary = df[~is_reference].sort_values("test_rmse").reset_index(drop=True)
    reference = df[is_reference].sort_values("test_rmse").reset_index(drop=True)

    header = (
        "| # | Model | Omics | Drug rep | Graph | Pairs | RMSE | MAE | R² | PCC | SCC | AUC | F1 | Params | Fit (s) | Details |\n"
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"
    )
    rows = "\n".join(fmt_row(r) for _, r in primary.iterrows())
    reference_rows = "\n".join(fmt_row(r) for _, r in reference.iterrows())

    seeds = pd.read_csv(SEED_VARIANCE_CSV) if SEED_VARIANCE_CSV.exists() else pd.DataFrame()
    seed_rows = "\n".join(
        # label contains " | " separators, which would break the table columns
        f"| {str(r['label']).replace(' | ', ', ')} | {int(r['n_seeds'])} | {r['recorded_rmse']:.4f} | "
        f"{r['rmse_mean']:.4f} | {r['rmse_std']:.4f} | {r['rmse_min']:.4f} | {r['rmse_max']:.4f} |"
        for _, r in seeds.iterrows()
    )

    best = primary.iloc[0]
    near_best = primary[primary["test_rmse"] < best["test_rmse"] + NOISE_THRESHOLD]
    per_family = primary.loc[primary.groupby("model")["test_rmse"].idxmin()].sort_values("test_rmse")
    family_rows = "\n".join(
        f"| {r['model']} | {r['omics']} | {r['drug_rep']} | {r['test_rmse']:.4f} | "
        f"{r['test_pcc']:.4f} | {r['test_r2']:.4f} |"
        for _, r in per_family.iterrows()
    )

    ens = pd.read_csv(ENSEMBLE_CSV) if ENSEMBLE_CSV.exists() else pd.DataFrame()
    ens_rows = "\n".join(
        # config contains " | " separators, which would break the table columns
        f"| {r['base_model']} | {str(r['config']).replace(' | ', ', ')} | {r['base_test_rmse']:.4f} | "
        f"{r['refined_test_rmse']:.4f} | {r['rmse_delta']:+.4f} | {r['rmse_pct_change']:+.2f}% |"
        for _, r in ens.iterrows()
    )

    content = f"""# Master Results — Full Experiment Matrix

Auto-generated by [`experiments/build_results_table.py`](../experiments/build_results_table.py)
from the per-family result CSVs. **Do not hand-edit** — re-run the script instead.

**{len(df)} runs** ({len(primary)} ranked below, {len(reference)} drug-identity
reference arms in a separate section), all on an identical
`GroupShuffleSplit`-by-cell-line 70/15/15 partition (`random_state=42`),
532 tri-omics-complete cell lines. Lower RMSE/MAE is better; higher
R²/PCC/SCC/AUC/F1 is better.

## Read this before comparing any two rows

Every row is a **single run**, and repeated runs of the same configuration vary
by **±0.013–0.029 RMSE** (a 0.075 range across five seeds; see
[13_seed_variance_results](./13_seed_variance_results.md)). **Differences below
~{NOISE_THRESHOLD:.2f} RMSE are not distinguishable from run-to-run variance.**
{len(near_best)} of the {len(primary)} ranked runs sit within that margin of the
top row, so the ordering among them is not a ranking.

Two findings originally drawn from single runs did not survive repetition —
molecular graphs beating fingerprints, and proteomics-alone being best. Both
reversed:

| Configuration | n seeds | Recorded | Mean | Std | Min | Max |
|---|---|---|---|---|---|---|
{seed_rows}

## Baselines

| Baseline (no model) | Test RMSE |
|---|---|
| Global training mean, 134,764-pair population | 2.7097 |
| Global training mean, 111,799-pair population | 2.7690 |
| **Per-drug training mean** (no omics at all), 111,799 pairs | **{PER_DRUG_MEAN_RMSE:.4f}** |

The per-drug mean is the baseline that matters: **drug identity alone accounts
for 71% of the reducible error.** Measured against it, the best model here
explains roughly a quarter of the remaining, omics-dependent variance. Measured
against the global mean instead, the same model shows R² ≈ 0.78 — which mostly
reflects knowing which compound was screened, not the multi-omics profile.

## Best per model family

| Model | Omics | Drug rep | RMSE | PCC | R² |
|---|---|---|---|---|---|
{family_rows}

**Lowest single-run RMSE: {best['exp_id']} — {best['model']}, {best['omics']},
{best['drug_rep']} — {best['test_rmse']:.4f}.** This is *not* "the best model":
re-run across five seeds, this configuration averages **1.3273 ± 0.0286**, the
worst of the three configurations tested, and its recorded value is the
favourable end of its own distribution. See
[13_seed_variance_results](./13_seed_variance_results.md) before quoting any
single number from this table.

## Ensemble refinement (XGBoost, applied post-hoc)

| Base model | Config | RMSE without | RMSE with | Δ | Change |
|---|---|---|---|---|---|
{ens_rows}

Refinement made both models worse — see
[ensemble_refinement_results](./08_ensemble_refinement_results.md) for why
(the paper's +19.7% gain likely depends on their random split allowing
cell-line leakage).

## All runs, ranked by test RMSE

{header}
{rows}

## Drug-identity reference arms (not ranked)

These runs represent a drug by its **identity**, not its structure. They are
kept for reference and for the population control they provide, but they are
deliberately excluded from the ranking above: a one-hot drug vector carries no
information about a compound the model has not seen, so these numbers cannot
support the use case the project targets. On unseen compounds one-hot collapses
to R² 0.024 while fingerprints hold at 0.422
([leave_drugs_out_results](./10_leave_drugs_out_results.md)).

They also score well here for a reason that flatters them: drug identity alone
accounts for 71% of the reducible error under this split, and a one-hot vector
hands that to the model directly.

{header}
{reference_rows}

- **`onehot`** — 295-dim drug identity, all 134,764 pairs. Matches the
  historical baselines.
- **`onehot_restricted`** — same features, restricted to the 111,799 pairs
  whose drug has a resolvable SMILES. The **population control**: comparing
  this against `fingerprint` isolates the representation effect, since only
  498 of 621 GDSC drug IDs resolved to a structure.

Comparing `onehot` directly against `fingerprint` conflates representation with
population and overstates the fingerprint effect; use `onehot_restricted` as the
baseline for that comparison.

## Caveats that apply to every row

- All cells use the 532 cell lines with complete tri-omics coverage (~55% of
  GDSC2's 969), including single-modality cells — this isolates the modality
  effect but means single-omics rows are not "what you'd get using every cell
  line that modality covers."
- The split is grouped by **cell line**, not by drug, so these numbers measure
  generalization to unseen cell lines. They say nothing about generalization to
  unseen *drugs* — the case where Morgan fingerprints matter most. That is
  measured separately in
  [leave_drugs_out_results](./10_leave_drugs_out_results.md): on compounds never
  screened during training, one-hot collapses to R² 0.024 while fingerprints
  hold at R² 0.422. **Do not read the fingerprint rows below as evidence
  against the fingerprint representation** — they measure the axis on which it
  is not expected to win.
- Numbers are **not comparable to MoGraphDRP's published RMSE 0.6622**, which
  uses a random 80/10/10 split where the same cell line can appear in both
  training and test. Measured under *their* protocol our best model reaches
  RMSE 0.8971 — so of the 0.582 apparent gap, 0.347 (60%) is protocol and
  0.235 (40%) is genuine architectural difference. See
  [split_protocol_comparison](./09_split_protocol_comparison.md).
"""

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(content, encoding="utf-8")
    print(f"Wrote {len(df)} runs -> {OUTPUT}")
    print(f"Best: {best['exp_id']} {best['model']} / {best['omics']} / {best['drug_rep']} "
          f"RMSE={best['test_rmse']:.4f}")


if __name__ == "__main__":
    main()
