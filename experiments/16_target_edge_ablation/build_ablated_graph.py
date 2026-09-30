"""Remove the drug->target edges of 7 probe drugs and save the result as a separate graph.

Drug nodes reach the protein graph only through `targets`, so a probe drug
with its edges removed is fully cut off from every protein. Everything else
(PPI edges, mutation edges, node features, other drugs' target edges) is left
byte-identical, so any change in the retrained model comes from these edges.

Run from the repository root:
    python experiments/16_target_edge_ablation/build_ablated_graph.py
"""

from __future__ import annotations

from pathlib import Path

import torch

REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_GRAPH = REPO_ROOT / "src" / "graph" / "hetero_graph.pt"
ABLATED_GRAPH = REPO_ROOT / "experiments" / "16_target_edge_ablation" / "hetero_graph_ablated.pt"
TARGET_EDGE = ("drug", "targets", "protein")

PROBE_DRUGS = {
    "1036": ("PLX-4720", "BRAF"),
    "1061": ("SB590885", "BRAF"),
    "1373": ("Dabrafenib", "BRAF"),
    "1010": ("Gefitinib", "EGFR"),
    "1168": ("Erlotinib", "EGFR"),
    "1915": ("AZD3759", "EGFR"),
    "1919": ("Osimertinib", "EGFR"),
}


def main() -> int:
    graph = torch.load(SOURCE_GRAPH, weights_only=False)
    drug_ids = [str(d) for d in graph["drug"].node_ids]
    protein_ids = list(graph["protein"].node_ids)

    missing = [d for d in PROBE_DRUGS if d not in drug_ids]
    if missing:
        raise ValueError(f"Probe drugs absent from the graph: {missing}")
    probe_idx = torch.tensor([drug_ids.index(d) for d in PROBE_DRUGS])

    edge_index = graph[TARGET_EDGE].edge_index
    removed = torch.isin(edge_index[0], probe_idx)
    for src, dst in edge_index[:, removed].t().tolist():
        name, gene = PROBE_DRUGS[drug_ids[src]]
        print(f"[remove] {drug_ids[src]:>5} {name:<12} -> {protein_ids[dst]}  (annotated target {gene})")

    # Every probe must lose at least one edge, or it was never wired and the ablation is vacuous for it.
    cut = {drug_ids[s] for s in edge_index[0, removed].tolist()}
    unwired = set(PROBE_DRUGS) - cut
    if unwired:
        raise ValueError(f"Probe drugs had no target edge to remove: {sorted(unwired)}")

    graph[TARGET_EDGE].edge_index = edge_index[:, ~removed]
    print(f"\n[targets] {edge_index.shape[1]} -> {graph[TARGET_EDGE].edge_index.shape[1]} edges "
          f"({int(removed.sum())} removed across {len(PROBE_DRUGS)} probe drugs)")
    for edge_type in graph.edge_types:
        print(f"  {edge_type}: {graph[edge_type].edge_index.shape[1]}")

    torch.save(graph, ABLATED_GRAPH)
    print(f"\nSaved -> {ABLATED_GRAPH.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
