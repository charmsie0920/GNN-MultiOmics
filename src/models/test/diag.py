import pandas as pd
import torch

graph = torch.load("src/graph/hetero_graph.pt", weights_only=False)
df = pd.read_csv("data/aligned/gdsc2_response_master.csv")

print("--- GRAPH NODE ID SAMPLES ---")
print("Cell Line Graph Sample:", list(graph["cell_line"].node_ids)[:3])
print("Drug Graph Sample:     ", list(graph["drug"].node_ids)[:3])
print("Drug Graph ID Types:   ", type(graph["drug"].node_ids[0]))

print("\n--- CSV COLUMN SAMPLES ---")
print("CSV standard_model_id:", list(df["standard_model_id"])[:3])
print("CSV Sanger ID:        ", list(df["sanger_model_id"])[:3])
print("CSV drug_id:          ", list(df["drug_id"])[:3])
print("CSV drug_id Type:     ", type(df["drug_id"].iloc[0]))

print(graph.edge_types)
print(graph.node_types)