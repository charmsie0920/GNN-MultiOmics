import pandas as pd
import torch
from sklearn.model_selection import GroupShuffleSplit

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

# 6. Train / Val / Test Splits (70 / 15 / 15), grouped by cell line so no
# cell line appears in more than one split (a plain random split leaks the
# same cell line's other drug-response rows into both train and test).
groups = df_valid[sanger_col].to_numpy()
gss1 = GroupShuffleSplit(n_splits=1, test_size=0.3, random_state=42)
train_idx, rest_idx = next(gss1.split(df_valid, groups=groups))

gss2 = GroupShuffleSplit(n_splits=1, test_size=0.5, random_state=42)
rel_val, rel_test = next(gss2.split(df_valid.iloc[rest_idx], groups=groups[rest_idx]))
val_idx, test_idx = rest_idx[rel_val], rest_idx[rel_test]

train_df, val_df, test_df = df_valid.iloc[train_idx], df_valid.iloc[val_idx], df_valid.iloc[test_idx]
df_valid.loc[train_df.index, "split"] = "train"
df_valid.loc[val_df.index, "split"] = "val"
df_valid.loc[test_df.index, "split"] = "test"

overlap = set(groups[train_idx]) & (set(groups[val_idx]) | set(groups[test_idx]))
assert not overlap, f"cell line leaked across splits: {sorted(overlap)[:5]}"

# 7. Save Processed Dataset
df_valid.to_csv("data/raw/aligned_ic50_pairs.csv", index=False)

print(f"Dataset successfully created!")
print(f"Matched pairs: {len(df_valid)} / {len(df)}")
print(f"Splits -> Train: {len(train_df)} | Val: {len(val_df)} | Test: {len(test_df)}")