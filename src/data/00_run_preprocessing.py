import pandas as pd
from pathlib import Path
from omics_preprocessing import OmicsPreprocessingPipeline, GE_KEY, MUT_CNV_KEY, PROTEOMICS_KEY

# 1. Define paths
WIDE_DIR = Path("data/processed/wide")
FINAL_OUT_DIR = Path("data/processed")
FINAL_OUT_DIR.mkdir(parents=True, exist_ok=True)

def main():
    print("Loading wide matrices...")
    omics_data = {
        GE_KEY: pd.read_csv(WIDE_DIR / "GE_wide.csv", index_col=0),
        MUT_CNV_KEY: pd.read_csv(WIDE_DIR / "Mut_CNV_wide.csv", index_col=0),
        PROTEOMICS_KEY: pd.read_csv(WIDE_DIR / "Proteomics_wide.csv", index_col=0)
    }

    # 2. Initialize Aditya's Pipeline (compressing to 128 dimensions)
    print("Running OmicsPreprocessingPipeline (Imputation, Z-score, PCA)...")
    pipeline = OmicsPreprocessingPipeline(d_target=128)
    
    # 3. Fit and transform the data
    transformed_data = pipeline.fit_transform(omics_data)
    
    # Save the fitted pipelines for inference later
    pipeline.save()

    # 4. Save to the exact filenames Charmaine's early fusion script expects
    print("Saving final PCA-compressed CSVs for baseline fusion...")
    
    # Map Aditya's keys to Charmaine's expected filenames
    file_mapping = {
        GE_KEY: "transcriptomics_pca.csv",
        MUT_CNV_KEY: "genomics_pca.csv",
        PROTEOMICS_KEY: "proteomics_pca.csv"
    }

    for modality_key, file_name in file_mapping.items():
        # Convert numpy array back to DataFrame to save with cell line indices
        df_pca = pd.DataFrame(
            transformed_data[modality_key], 
            index=omics_data[modality_key].index
        )
        out_path = FINAL_OUT_DIR / file_name
        df_pca.to_csv(out_path)
        print(f"Saved: {out_path} with shape {df_pca.shape}")

    print("Data Engineering pipeline complete! Ready for 01_early_fusion.py.")

if __name__ == "__main__":
    main()