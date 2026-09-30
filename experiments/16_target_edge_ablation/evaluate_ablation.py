"""Compare models trained with and without the probe drugs' target edges.

Models (all scored on the same test split, data/raw/aligned_ic50_pairs.csv):
- production:      models/checkpoints/best_hetero_gnn.pt on the original graph
- zero_shot:       the same weights on the ablated graph (edge vanishes, no chance to adapt)
- control_seed42:  retrained at seed 42 on the original graph
- ablated_seed42:  retrained at seed 42 on the ablated graph

control_seed42 vs ablated_seed42 is the clean comparison: same seed, same
data, differing only in the 7 removed edges.

Attribution: for driver-mutant cell lines on their gene's inhibitors (BRAF,
EGFR), the Check 1 statistic -- Spearman(target protein's gradient x embedding
score, measured ln IC50) -- per model, against a floor of 50 untrained models
on the *ablated* graph (what wiring alone gives once the edge is gone).

Run from the repository root, after build_ablated_graph.py and both train_ablated.py arms:
    python experiments/16_target_edge_ablation/evaluate_ablation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "15_interpretability_validation"))

import random_baseline as rb  # noqa: E402
from build_ablated_graph import PROBE_DRUGS  # noqa: E402
from train_ablated import GRAPHS, HIDDEN_DIM, SPLIT_CSV, load_graph  # noqa: E402

from src.data.experiment_utils import compute_shared_threshold, evaluate  # noqa: E402
from src.models.train.hetero_gnn import HeteroIC50GNN  # noqa: E402

CHECKPOINTS = {
    "production": (REPO_ROOT / "models" / "checkpoints" / "best_hetero_gnn.pt", "original"),
    "zero_shot": (REPO_ROOT / "models" / "checkpoints" / "best_hetero_gnn.pt", "ablated"),
    "control_seed42": (HERE / "checkpoints" / "hetero_gnn_original_seed42.pt", "original"),
    "ablated_seed42": (HERE / "checkpoints" / "hetero_gnn_ablated_seed42.pt", "ablated"),
}
RMSE_CSV = HERE / "rmse_comparison.csv"
ATTRIBUTION_CSV = HERE / "attribution_comparison.csv"
N_CONTROLS = 50
TOP_K = 50


def graph_maps(name: str, device: torch.device) -> dict:
    x_dict, edge_index_dict = load_graph(GRAPHS[name], device)
    graph = torch.load(GRAPHS[name], weights_only=False)
    return {
        "x_dict": x_dict,
        "edge_index_dict": edge_index_dict,
        "cell_to_idx": {str(c): i for i, c in enumerate(graph["cell_line"].node_ids)},
        "drug_to_idx": {str(d): i for i, d in enumerate(graph["drug"].node_ids)},
        "num_proteins": int(x_dict["protein"].shape[0]),
        "protein_node_ids": [str(p) for p in graph["protein"].node_ids],
    }


def load_model(path: Path, maps: dict, device: torch.device) -> torch.nn.Module:
    model = HeteroIC50GNN(num_proteins=maps["num_proteins"], hidden_dim=HIDDEN_DIM).to(device)
    model.load_state_dict(torch.load(path, map_location=device))
    return model.eval()


@torch.no_grad()
def predict(model, maps: dict, frame: pd.DataFrame, device) -> np.ndarray:
    out = []
    for start in range(0, len(frame), 4096):
        chunk = frame.iloc[start : start + 4096]
        out.append(model(maps["x_dict"], maps["edge_index_dict"],
                         torch.tensor(chunk["cell_idx"].to_numpy(), device=device),
                         torch.tensor(chunk["drug_idx"].to_numpy(), device=device)).cpu().numpy())
    return np.concatenate(out)


def case_study_pairs(gene: str, split: pd.DataFrame) -> pd.DataFrame:
    """Every measured pair of a `gene`-driver-mutant line on one of that gene's probe inhibitors."""
    drugs = [d for d, (_, g) in PROBE_DRUGS.items() if g == gene]
    mutants = rb.driver_mutant_cell_lines(gene)
    pairs = split[split["drug_id"].isin(drugs) & split["sanger_model_id"].isin(mutants)]
    return pairs.drop_duplicates(["sanger_model_id", "drug_id"]).reset_index(drop=True)


def attribution_stats(scores: np.ndarray, protein_idx: int, pairs: pd.DataFrame) -> dict:
    target = scores[:, protein_idx]
    ranks = (np.abs(scores) > np.abs(target)[:, None]).sum(axis=1) + 1
    held_out = pairs["split"].isin(["val", "test"]).to_numpy()
    ln_ic50 = pairs["ln_ic50"].to_numpy()

    def rho(mask):
        return float(spearmanr(target[mask], ln_ic50[mask]).statistic) if np.std(target[mask]) > 0 else np.nan

    return {
        "spearman_all": rho(np.ones(len(pairs), dtype=bool)),
        "spearman_held_out": rho(held_out),
        "sensitising_share": float((target < 0).mean()),
        "target_in_top50_share": float((ranks <= TOP_K).mean()),
        "median_target_rank": float(np.median(ranks)),
        "mean_abs_target_score": float(np.abs(target).mean()),
    }


def main() -> int:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[setup] device={device}")
    maps = {name: graph_maps(name, device) for name in GRAPHS}
    split = pd.read_csv(SPLIT_CSV, usecols=["sanger_model_id", "drug_id", "cell_idx", "drug_idx", "ln_ic50", "split"],
                        dtype={"sanger_model_id": str, "drug_id": str})
    models = {name: (load_model(path, maps[g], device), g) for name, (path, g) in CHECKPOINTS.items()}

    # --- RMSE ---------------------------------------------------------------
    threshold = compute_shared_threshold()
    test = split[split["split"] == "test"].reset_index(drop=True)
    probe_gene = test["drug_id"].map({d: g for d, (_, g) in PROBE_DRUGS.items()})
    subsets = {
        "all test rows": np.ones(len(test), dtype=bool),
        "non-probe drugs": probe_gene.isna().to_numpy(),
        "probe drugs (BRAF+EGFR)": probe_gene.notna().to_numpy(),
        "BRAF probes": (probe_gene == "BRAF").to_numpy(),
        "EGFR probes": (probe_gene == "EGFR").to_numpy(),
    }
    rmse_rows = []
    for name, (model, g) in models.items():
        pred = predict(model, maps[g], test, device)
        for label, mask in subsets.items():
            metrics = evaluate(test["ln_ic50"].to_numpy()[mask], pred[mask], threshold)
            rmse_rows.append({"model": name, "graph": g, "subset": label, "n": int(mask.sum()), **metrics})
    rmse = pd.DataFrame(rmse_rows)
    rmse.to_csv(RMSE_CSV, index=False)

    # --- Attribution ----------------------------------------------------------
    attribution_rows = []
    for gene in ("BRAF", "EGFR"):
        pairs = case_study_pairs(gene, split)
        protein_idx = rb.protein_index_for_symbol(maps["original"], gene)
        assert maps["ablated"]["protein_node_ids"][protein_idx] == maps["original"]["protein_node_ids"][protein_idx]
        cells, drugs = pairs["cell_idx"].tolist(), pairs["drug_idx"].tolist()
        print(f"[{gene}] {len(pairs)} pairs, {pairs['sanger_model_id'].nunique()} cell lines "
              f"({pairs.loc[pairs['split'] != 'train', 'sanger_model_id'].nunique()} held out)")

        for name, (model, g) in models.items():
            scores = rb.score_pairs(model, maps[g], cells, drugs).cpu().numpy()
            attribution_rows.append({"gene": gene, "model": name, "graph": g, "n_pairs": len(pairs),
                                     **attribution_stats(scores, protein_idx, pairs)})

        control_stats = []
        for model in rb.build_random_models(maps["ablated"], device, n=N_CONTROLS):
            scores = rb.score_pairs(model, maps["ablated"], cells, drugs).cpu().numpy()
            control_stats.append(attribution_stats(scores, protein_idx, pairs))
        controls = pd.DataFrame(control_stats)
        for stat in ("mean", "std"):
            attribution_rows.append({"gene": gene, "model": f"untrained_x{N_CONTROLS}_{stat}", "graph": "ablated",
                                     "n_pairs": len(pairs), **getattr(controls, stat)(numeric_only=True).to_dict()})
        # Share of controls matching or beating each trained model's correlation.
        for row in attribution_rows:
            if row["gene"] == gene and row["model"] in models:
                row["empirical_p_vs_untrained"] = float((controls["spearman_all"] >= row["spearman_all"]).mean())

    attribution = pd.DataFrame(attribution_rows)
    attribution.to_csv(ATTRIBUTION_CSV, index=False)

    with pd.option_context("display.width", 220, "display.max_columns", None, "display.float_format", "{:.4f}".format):
        print("\n" + "=" * 110 + "\nRMSE on the test split\n" + "=" * 110)
        print(rmse.pivot(index="subset", columns="model", values="rmse")[list(CHECKPOINTS)]
              .loc[list(subsets)].to_string())
        print("\nPCC:")
        print(rmse.pivot(index="subset", columns="model", values="pcc")[list(CHECKPOINTS)]
              .loc[list(subsets)].to_string())
        print("\n" + "=" * 110 + "\nTarget-protein attribution on driver-mutant x own-inhibitor pairs\n" + "=" * 110)
        print(attribution.to_string(index=False))
    print(f"\nSaved {RMSE_CSV.name}, {ATTRIBUTION_CSV.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
