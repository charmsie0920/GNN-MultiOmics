import pandas as pd
import torch
from sklearn.model_selection import train_test_split

# 1. Load Graph and Master Response
graph = torch.load("src/graph/hetero_graph.pt", weights_only=False)
df = pd.read_csv("data/aligned/gdsc2_response_master.csv")

# 2. Build Lookup Tables (all keys as strings)
cell_to_idx = {str(cell_id): idx for idx, cell_id in enumerate(graph["cell_line"].node_ids)}
drug_to_idx = {str(drug_id): idx for idx, drug_id in enumerate(graph["drug"].node_ids)}

# 3. Determine Correct Sanger Column in CSV
sanger_col = "sanger_model_id" if "sanger_model_id" in df.columns else "gdsc_sanger_model_id"

# 4. Map CSV Values to Graph Indices
df["cell_idx"] = df[sanger_col].astype(str).map(cell_to_idx)
df["drug_idx"] = df["drug_id"].astype(str).map(drug_to_idx)

# 5. Drop Unmatched Rows & Cast Indices to Int
df_valid = df.dropna(subset=["cell_idx", "drug_idx"]).copy()
df_valid["cell_idx"] = df_valid["cell_idx"].astype(int)
df_valid["drug_idx"] = df_valid["drug_idx"].astype(int)

# Identify Target Label Column
target_col = "nlme_result" if "nlme_result" in df_valid.columns else df_valid.columns[-1]

# 6. Train / Val / Test Splits (80 / 10 / 10)
train_df, test_df = train_test_split(df_valid, test_size=0.2, random_state=42)
val_df, test_df = train_test_split(test_df, test_size=0.5, random_state=42)

df_valid.loc[train_df.index, "split"] = "train"
df_valid.loc[val_df.index, "split"] = "val"
df_valid.loc[test_df.index, "split"] = "test"

# 7. Save Processed Dataset
df_valid.to_csv("data/raw/aligned_ic50_pairs.csv", index=False)

print(f"Dataset successfully created!")
print(f"Matched pairs: {len(df_valid)} / {len(df)}")
print(f"Splits -> Train: {len(train_df)} | Val: {len(val_df)} | Test: {len(test_df)}")