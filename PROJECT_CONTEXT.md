# PROJECT CONTEXT: Multi-Omics Graph Neural Network (MCS16) for Drug Sensitivity Prediction

## 1. Project Overview & Objective

- **Goal:** Build a multi-modal Graph Neural Network (GNN) platform that predicts cellular drug response ($IC_{50}$ / $LN(IC_{50})$) by fusing multi-omics cell line profiles with chemical compound topologies.
- **Key Target:** Achieve a lower Root Mean Squared Error (RMSE) than standard flat vector concatenation models while operating within local consumer GPU hardware limits.
- **Core Architecture:**
  1. **Omics Input:** Transcriptomics (RNA-seq TPM), Somatic Mutations (MAF/binary indicators), Copy Number Variations (CNV).
  2. **Chemical Input:** SMILES strings converted into molecular graphs via RDKit.
  3. **Biological Network:** Protein-Protein Interaction (PPI) topology extracted from the STRING database.
  4. **Deep Learning Engine:** PyTorch Geometric (PyG) heterogeneous GNN utilizing multi-head cross-attention feature fusion and Graph Attention Network (GAT) message passing layers.

---

## 2. Team Member Roles & Code Ownership

- **[AA] Bioinformatics Data Engineer (Current User Scope):**
  - **Ownership:** Data ingestion, batch effect correction, imputation, dimension projection (PCA), STRING network extraction, and serializing PyTorch Geometric graph data objects (`.pt`).
- **[WH] Machine Learning Engineer:**
  - **Ownership:** PyTorch Geometric heterogeneous GNN architecture, multi-head cross-attention fusion tensor math, GAT message-passing scripts, $IC_{50}$ regression training.
- **[CC] Technical & Infrastructure Lead:**
  - **Ownership:** End-to-end integration framework, cloud execution clusters, strictly enforcing local hardware VRAM footprints (**4GB to 8GB VRAM limit**).
- **[SM] Full-Stack & Data Visualization Engineer:**
  - **Ownership:** Standalone dashboard GUI, asynchronous REST API integration, interactive node-link sub-graph neighborhood rendering.

---

## 3. Data Engineering Strategy & Dataset Specs (Pending Local Upload)

The Data Engineer is opting for **frozen bulk downloads** instead of live REST APIs to ensure pipeline stability, reproducibility, and prevent rate limits.

### Incoming Files to Expect in Data Directory (`./data/raw/`):

1. **Broad DepMap / CCLE (Multi-Omics Features):**
   - `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (Transcriptomics TPM)
   - `OmicsSomaticMutations.csv` (Genomic mutations)
   - `OmicsCNV.csv` (Copy number variations)
   - `Model.csv` (Metadata mapping `ModelID` e.g., `ACH-000162` to `COSMIC_ID` / cell lines)
2. **Sanger GDSC (Drug Response Labels):**
   - `GDSC1_fitted_dose_response.csv` or `GDSC2_fitted_dose_response.csv` (Mapping `COSMIC_ID` to drug names and $LN(IC_{50})$ values)
3. **STRING Database (Network Topology):**
   - `9606.protein.links.v12.0.txt.gz` (Human PPI interaction scores)

### Data Processing Pipeline Sequence:

1. **Ingestion & ID Matching:** Join GDSC `COSMIC_ID` labels with Broad `ModelID` features via `Model.csv`.
2. **Cleaning & Batch Correction:** KNN Imputation for missing values $\rightarrow$ ComBat / $z$-score normalization to remove laboratory artifacts.
3. **Dimension Projections:** Binarize mutations ($0/1$) $\rightarrow$ Run independent PCA pipelines on continuous omics layers to reduce column shape to a uniform hidden vector dimension ($D$).
4. **Graph Construction:** Filter STRING PPI network to match selected genes $\rightarrow$ Generate coordinate-format sparse index arrays (`edge_index` of shape `[2, Num_Edges]`).
5. **Handoff Export:** Package into `torch_geometric.data.Data` objects and export as binary serialized PyTorch tensors (`.pt`) to `./data/processed/`.

---

## 4. Coding & Environment Constraints

- **Language & Frameworks:** Python 3.10+, PyTorch Geometric (PyG), PyTorch, RDKit, Pandas, Scikit-learn, NeuroCombat / Combat.
- **Hardware Footprint Envelope:** Code routines (especially tensor batching and PCA matrix operations) must be memory-efficient, targeting **4GB to 8GB GPU VRAM** and standard CPU RAM limits.
- **Formatting Standards:**
  - Use modular python scripts (`src/data/`, `src/models/`).
  - Strict tensor shape logging (e.g., print shape checks after PCA projections and before PyG object compilation).
  - Explicit function docstrings adhering to standard scientific python conventions.

---

## 5. Next Action

- **Awaiting Upload:** User will upload the downloaded dataset bulk files into the workspace directory (`./data/raw/`).
- **Immediate Task for Agent:** Stand by to write the initial Python ingestion and cleaning script (`src/data/ingest_and_align.py`) to process `Model.csv`, GDSC dose response tables, and CCLE matrices once files are placed.
