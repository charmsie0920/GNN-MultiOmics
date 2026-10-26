"""Figures for the final written report (FIT3162 Report 2).

Reads only committed result files and the git history; writes PNGs next to
this script. Re-run from the repository root:

    .venv/bin/python reports/final_report/figures/make_report_figures.py

Figures:
    fig_ablation_ladder.png    test RMSE change vs `base`, per rung, paired by seed
    fig_protocols.png          `base` against the final model under three splits
    fig_direct_hit.png         signed error on pairs whose drug target is mutated
    fig_commit_timeline.png    commits per week against the planned milestones
"""

from __future__ import annotations

import subprocess
from datetime import date, timedelta
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "src" / "final_model" / "results"
OUT = Path(__file__).resolve().parent

# Reference palette (dataviz skill), validated light-mode: blue / orange.
BLUE = "#2a78d6"
ORANGE = "#eb6834"
GREY = "#8a8985"
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"

NOISE_BAND = 0.03
PER_DRUG_FLOOR = 1.4889
FINAL = "base+proteomics+std_targets+mol_graph+pair_features"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.edgecolor": GRID,
        "axes.labelcolor": TEXT_2,
        "axes.titlecolor": TEXT,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "xtick.color": TEXT_2,
        "ytick.color": TEXT_2,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def _save(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUT / name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / name}")


def paired_delta(rows: pd.DataFrame, config: str, baseline: str = "base") -> pd.Series:
    """Per-seed test RMSE of `config` minus `baseline`, on the seeds both have."""
    cfg = rows[rows.config == config].set_index("seed").test_rmse
    base = rows[rows.config == baseline].set_index("seed").test_rmse
    shared = cfg.index.intersection(base.index)
    return cfg.loc[shared] - base.loc[shared]


def fig_ablation_ladder() -> None:
    ladder = pd.read_csv(RESULTS / "ablation_results.csv")
    gate = pd.read_csv(RESULTS / "pair_gate_results.csv")
    rows = pd.concat([ladder, gate[gate.config != "base"]], ignore_index=True)
    rows = rows[rows.protocol == "grouped"]

    # (config, label, verdict); verdict picks the colour and the printed tag.
    rungs = [
        ("base+proteomics", "+ proteomics", "kept"),
        ("base+std_targets", "+ per-drug standardisation", "kept"),
        ("base+mol_graph", "+ molecular-graph GCN", "none"),
        ("base+pair_features", "+ PPI pair features", "none"),
        ("base+cross_attention", "+ cross-attention fusion", "rejected"),
        ("base+pair_attention", "+ pair attention", "rejected"),
        ("base+bilinear", "+ bilinear head", "rejected"),
        ("aligned", "MoGraphDRP-aligned", "none"),
        ("full", "all five Phase 2 parts", "none"),
        ("base+proteomics+std_targets", "proteomics + standardisation", "kept"),
        (FINAL, "final model", "final"),
    ]
    colours = {"kept": BLUE, "final": BLUE, "none": GREY, "rejected": ORANGE}
    tags = {
        "kept": "kept",
        "final": "final",
        "none": "no effect beyond noise",
        "rejected": "rejected",
    }

    clip = 0.16
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.axvspan(-NOISE_BAND, NOISE_BAND, color=GRID, zorder=0)
    ax.axvline(0, color=TEXT_2, linewidth=0.8, zorder=1)
    ax.text(0, len(rungs) - 0.45, "noise band ±0.03", ha="center", va="top", color=TEXT_2, fontsize=7)

    for y, (config, label, verdict) in enumerate(rungs):
        delta = paired_delta(rows, config)
        mean, std, n = delta.mean(), delta.std(), len(delta)
        better = int((delta < 0).sum())
        shown = min(mean, clip)
        ax.barh(y, shown, height=0.62, color=colours[verdict], zorder=2)
        if mean <= clip:
            ax.errorbar(shown, y, xerr=std, color=TEXT, elinewidth=0.9, capsize=2.5, zorder=3)
        note = f"{mean:+.3f}  ({better}/{n} seeds better) · {tags[verdict]}"
        if mean > clip:
            note = f"{mean:+.3f} (axis clipped)  ({better}/{n}) · {tags[verdict]}"
        x_text = max(shown, 0) + (std if mean <= clip else 0) + 0.006
        ax.text(x_text, y, note, va="center", fontsize=7.5, color=TEXT)

    ax.set_yticks(range(len(rungs)), [r[1] for r in rungs])
    ax.invert_yaxis()
    ax.set_xlim(-0.09, 0.42)
    ax.set_xlabel("Change in test RMSE vs base (ln IC50); negative is better")
    ax.grid(axis="x", color=GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    base = rows[rows.config == "base"].test_rmse
    ax.set_title(
        f"Additive ablation, grouped split: base test RMSE {base.mean():.4f} ± {base.std():.4f}, "
        f"per-drug floor {PER_DRUG_FLOOR}",
        fontsize=9,
    )
    _save(fig, "fig_ablation_ladder.png")


def fig_protocols() -> None:
    grouped = pd.read_csv(RESULTS / "ablation_results.csv")
    random_split = pd.read_csv(RESULTS / "random_split_results.csv")
    ldo = pd.read_csv(RESULTS / "ldo_results.csv")

    panels = [
        ("Unseen cell lines\n(grouped split, headline)", grouped, "base", FINAL, PER_DRUG_FLOOR, "per-drug mean"),
        ("Random pairs\n(sanity check)", random_split, "base", "final", None, None),
        ("Unseen drugs\n(leave-drugs-out)", ldo, "base", "base+proteomics+mol_graph+pair_features", 2.5500, "mean only"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9))
    for ax, (title, frame, base_cfg, final_cfg, floor, floor_label) in zip(axes, panels):
        stats = []
        for cfg in (base_cfg, final_cfg):
            values = frame[frame.config == cfg].test_rmse
            stats.append((values.mean(), values.std(), len(values)))
        colours = [GREY, BLUE]
        names = ["base", "final" if final_cfg != panels[2][3] else "final, no std."]
        for x, ((mean, std, n), colour, name) in enumerate(zip(stats, colours, names)):
            ax.bar(x, mean, width=0.6, color=colour, zorder=2)
            ax.errorbar(x, mean, yerr=std, color=TEXT, elinewidth=0.9, capsize=2.5, zorder=3)
            ax.text(x, 0.08, f"{mean:.3f}\nn={n}", ha="center", va="bottom", fontsize=7.5, color="#ffffff", fontweight="bold")
        if floor:
            ax.axhline(floor, color=ORANGE, linestyle="--", linewidth=1.2, zorder=4)
            ax.text(-0.45, floor + 0.04, f"{floor_label} floor {floor:.4f}", ha="left", va="bottom", fontsize=7, color=TEXT_2)
        ax.set_xticks([0, 1], names)
        ax.set_ylim(0, 2.9)
        ax.set_title(title, fontsize=8.5)
        ax.grid(axis="y", color=GRID, linewidth=0.6, zorder=0)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Test RMSE (ln IC50)")
    fig.tight_layout()
    _save(fig, "fig_protocols.png")


def fig_direct_hit() -> None:
    subsets = pd.read_csv(RESULTS / "ablation_results_subsets_graph_decision.csv")
    subsets = subsets[subsets.subset == "direct_hit"]
    configs = [
        ("base+proteomics+std_targets+mol_graph", "without pair features", GREY),
        (FINAL, "final (with pair features)", BLUE),
    ]
    folds = [("val", "Validation (213 pairs)"), ("test", "Test (211 pairs)")]

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
    for ax, metric, label in zip(axes, ["mean_error", "rmse"], ["Mean signed error (ln IC50)", "RMSE (ln IC50)"]):
        width = 0.36
        for i, (cfg, name, colour) in enumerate(configs):
            for j, (fold, _) in enumerate(folds):
                values = subsets[(subsets.config == cfg) & (subsets.fold == fold)][metric]
                x = j + (i - 0.5) * width
                ax.bar(x, values.mean(), width=width - 0.02, color=colour, label=name if j == 0 else None, zorder=2)
                ax.errorbar(x, values.mean(), yerr=values.std(), color=TEXT, elinewidth=0.9, capsize=2.5, zorder=3)
                offset = 0.03 if values.mean() >= 0 else -0.03
                ax.text(x, values.mean() + values.std() * (1 if offset > 0 else -1) + offset, f"{values.mean():+.2f}" if metric == "mean_error" else f"{values.mean():.2f}",
                        ha="center", va="bottom" if offset > 0 else "top", fontsize=7.5, color=TEXT)
        ax.axhline(0, color=TEXT_2, linewidth=0.8)
        ax.set_xticks(range(len(folds)), [f[1] for f in folds])
        ax.set_ylabel(label)
        ax.grid(axis="y", color=GRID, linewidth=0.6, zorder=0)
        ax.set_axisbelow(True)
    axes[0].set_ylim(-0.25, 0.65)
    axes[1].set_ylim(0, 2.4)
    axes[0].set_title("Over-prediction removed", fontsize=8.5)
    axes[1].set_title("Error on direct-hit pairs", fontsize=8.5)
    axes[1].legend(frameon=False, fontsize=7.5, loc="upper left", ncol=2, bbox_to_anchor=(-0.02, 1.0))
    fig.suptitle("Pairs whose drug target is mutated in the cell line (1.26% of pairs), 5 seeds", x=0.01, ha="left",
                 fontsize=9, fontweight="bold", color=TEXT)
    fig.tight_layout()
    _save(fig, "fig_direct_hit.png")


def fig_commit_timeline() -> None:
    log = subprocess.run(
        ["git", "log", "--all", "--no-merges", "--date=short", "--format=%ad"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout.split()
    days = pd.to_datetime(pd.Series(log))
    start = date(2026, 4, 20)
    end = date(2026, 11, 2)
    weeks = days.dt.to_period("W-SUN").dt.start_time
    counts = weeks.value_counts()
    index = pd.date_range(start, end, freq="W-MON")
    counts = counts.reindex(index, fill_value=0)

    fig, ax = plt.subplots(figsize=(7.4, 2.8))
    ax.bar(counts.index, counts.values, width=5, color=BLUE, zorder=2)
    milestones = [
        (date(2026, 8, 24), "data freeze 24 Aug", True, 0),
        (date(2026, 9, 9), "mid-sem pitch 9 Sep", True, 1),
        (date(2026, 9, 28), "professor review 28 Sep", True, 2),
        (date(2026, 10, 12), "final model chosen 12 Oct", False, 3),
        (date(2026, 11, 2), "report due 2 Nov", True, 4),
    ]
    top = counts.max() + 2 + 3 * len(milestones)
    for when, label, planned, row in milestones:
        ax.axvline(pd.Timestamp(when), color=ORANGE if not planned else TEXT_2,
                   linestyle="-" if not planned else "--", linewidth=1.0, zorder=3)
        ax.text(pd.Timestamp(when) - timedelta(days=1.5), top - 3 * row, label, fontsize=7, color=TEXT_2,
                va="top", ha="right", bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.6))
    ax.annotate("no commits\n4 May – 24 Jul", xy=(pd.Timestamp(2026, 6, 14), 1), ha="center", fontsize=7.5, color=TEXT_2)
    ax.set_ylabel("Commits per week\n(all branches, no merges)")
    ax.set_ylim(0, top + 1)
    ax.grid(axis="y", color=GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title(f"Repository activity, {len(log)} commits; dashed = planned milestone, solid = actual", fontsize=9)
    fig.autofmt_xdate()
    _save(fig, "fig_commit_timeline.png")


def fig_app_architecture() -> None:
    """Box diagram of the delivered desktop application, serving the final model."""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    boxes = {
        "ui": (0.0, 2.3, 3.3, 1.7, "PySide6 desktop client", "UILauncher.py, pages/",
               "Dataset Initialisation · Execution Log\nFinal Results · Model Analytics"),
        "workers": (0.0, 0.2, 3.3, 1.6, "Client layer", "client/api_client.py, workers.py",
                    "QThread workers: upload, start run,\nstatus poller (1 s), results, attribution"),
        "api": (4.6, 2.3, 3.3, 1.7, "FastAPI backend", "backend/main.py, routers/",
                "/dataset/upload · /run/start\n/run/status · /results/*"),
        "runs": (4.6, 0.2, 3.3, 1.6, "Run manager + model backend", "backend/model_backends/",
                 "FinalModelBackend: 5-seed ensemble\nof final, rank_predictions()"),
        "data": (9.2, 2.3, 2.6, 1.7, "Local artefacts", "",
                 "omics PCA files, SMILES\nhetero_graph.pt, 5 checkpoints\nuploads/ (deleted on exit)"),
        "ext": (9.2, 0.2, 2.6, 1.6, "Enrichr (external)", "",
                "pathway enrichment of\nattributed genes, cached"),
    }
    fig, ax = plt.subplots(figsize=(9.6, 3.6))
    for key, (x, y, w, h, title, path, body) in boxes.items():
        colour = BLUE if key in ("ui", "workers") else (ORANGE if key in ("api", "runs") else GREY)
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                    facecolor=SURFACE, edgecolor=colour, linewidth=1.6))
        ax.text(x + 0.12, y + h - 0.15, title, fontsize=7.6, fontweight="bold", color=TEXT, va="top")
        if path:
            ax.text(x + 0.12, y + h - 0.45, path, fontsize=6.6, color=TEXT_2, va="top", family="DejaVu Sans Mono")
        ax.text(x + 0.12, y + h - 0.8, body, fontsize=6.8, color=TEXT_2, va="top", linespacing=1.35)

    def arrow(start, end, label, both=True, dy=0.08):
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="<|-|>" if both else "-|>", mutation_scale=9,
                                     color=TEXT_2, linewidth=1.0))
        if label:
            mx, my = (start[0] + end[0]) / 2, (start[1] + end[1]) / 2
            ax.text(mx, my + dy, label, fontsize=6.6, color=TEXT_2, ha="center", va="bottom")

    arrow((1.65, 2.3), (1.65, 1.8), "")
    arrow((3.3, 1.2), (4.6, 2.9), "", dy=0.0)
    ax.text(3.75, 2.1, "HTTP\nlocalhost", fontsize=6.6, color=TEXT_2, ha="right", va="center")
    arrow((6.25, 2.3), (6.25, 1.8), "")
    arrow((7.9, 1.4), (9.2, 3.0), "loads", both=False, dy=0.0)
    arrow((7.9, 0.8), (9.2, 0.8), "genes")
    ax.set_xlim(-0.1, 11.9)
    ax.set_ylim(0.0, 4.1)
    ax.axis("off")
    _save(fig, "fig_app_architecture.png")


if __name__ == "__main__":
    fig_app_architecture()
    fig_ablation_ladder()
    fig_protocols()
    fig_direct_hit()
    fig_commit_timeline()
