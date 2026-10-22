# GNN-Based Multi-Omics Integration for Predicting Cancer Drug Response

**Final Project Report** · FIT3162 Computer Science Project 2 · Team MCS16
Supervisor: Dr Ong Huey Fang · Submitted 2 November 2026

**Word count:** [fill in from Google Docs, sections 1–6 only; target ≤ 3,500]

## Team members

| Initials | Name | Student ID |
|---|---|---|
| AA | Aditya Ansh | 33520070 |
| CC | Charmaine Chee | 33520399 |
| WH | Wasif Hussain | 33621837 |
| SM | Shaahid Marikar | 34086242 |

## Authorship statement

> **PROPOSED, for the team to confirm or change before submission. A missing statement costs 50%.** Each member drafted the sections listed and edited and finalised other members' contributions to them.

- **AA (Aditya Ansh)** led Section 3.2 (Results), Section 3.4 (Key decisions) and Sections 4.2–4.4 (Data pipeline, Training, Evaluation protocol). AA prepared Figures 2, 3, 4 and 6 from the result files and Table D1 (Appendix D).
- **CC (Charmaine Chee)** led Section 2.2 (Related work), Section 4.1 (Design rationale) and Section 3.5 (Limitations and future work). CC prepared Figure 5 (final architecture) and the reference list, and was the document's proof editor.
- **SM (Shaahid Marikar)** led Section 1 (Introduction), Section 3.1 (Delivered product), Section 3.3 (Requirements) and Sections 4.5–4.6 (Software, tools and application architecture). SM prepared Figure 1, Table A1 (RTM), the screenshots in Appendix B and the sample code in Appendix C.
- **WH (Wasif Hussain)** led Section 2.1 (Context), Section 5 (Critique) and Section 6 (Conclusion). WH prepared Table 4 (proposal against outcome), the front matter and the word count.

Google Doc (edit log): [link]

## Table of contents

[Insert automatically in Google Docs: Insert → Table of contents.]

**Figures:**
1. Application architecture
2. Additive ablation ladder
3. Direct-hit pairs
4. Three evaluation protocols
5. Final model architecture
6. Repository activity against milestones

Appendix figures: B1–B6, application screenshots (Appendix B).

**Tables:**
1. Main results
2. Requirements summary
3. Tools
4. Proposal against outcome

Appendix tables: A1 (RTM), D1 (full ablation).

---

## 1. Introduction

Patients with the same cancer diagnosis often respond differently to the same drug. Matching therapy to a tumour's molecular profile pays off: across 570 phase II trials, biomarker-matched therapy produced a 31% response rate against 10.5% without it (Schwaederle et al., 2015). Testing every drug on every tumour is not possible, so computational drug-response prediction (DRP) estimates a drug's potency, measured as the half-maximal inhibitory concentration (IC50), from a tumour's molecular data and the drug's structure.

This project, MCS16, set out to predict ln(IC50) for cancer cell lines from three omics layers, gene expression (GE), mutation and copy-number variation (Mut+CNV) and proteomics, together with drug structure, and to put the predictions in a local desktop application for cancer pharmacology researchers and bioinformaticians. The proposal planned a graph attention network with cross-attention fusion that would improve on the published MoGraphDRP benchmark (Ahmadi et al., 2026).

We delivered three things:
- a desktop application that ranks the drugs in GDSC2 by predicted IC50 for a chosen cell line;
- a final model that predicts for cell lines never seen in training with a test RMSE of 1.2737, against 1.3252 for the reproduced MoGraphDRP-style base and 1.4889 for a per-drug-mean baseline;
- an additive ablation that measures every component we added, over five random seeds.

The ablation also shows which proposed components did not help, and the report says so.

Section 2 summarises the background and related work. Section 3 presents the product, the results, the requirements and the key decisions. Section 4 describes the method in enough detail to reproduce it. Section 5 critiques the project and its management, and Section 6 concludes.

## 2. Project background

### 2.1 Context

Precision oncology depends on reading a tumour's molecular profile, yet that information is under-used. In a 2026 cohort of 63,294 US patients with advanced cancer, many were never given genomic testing (Hage Chehade et al., 2026), and clinicians face more multi-omics data than one person can integrate (Aldea et al., 2023). Public pharmacogenomic screens make the computational side of the problem tractable. GDSC2 measures IC50 for hundreds of drugs across about a thousand cancer cell lines (Yang et al., 2013), and the Cell Model Passports provide matched omics for those lines (van der Meer et al., 2019). A model trained on these screens can rank candidate drugs for a cell line before costly wet-lab assays. That is the use case our application serves; clinical use on patients is out of scope.

### 2.2 Related work, updated since the proposal

The proposal reviewed flat-vector baselines (random forest, support vector regression), deep multi-omics models with concatenated features (Wang et al., 2021; Li et al., 2026), and graph models (Peng et al., 2022). It chose MoGraphDRP as the benchmark. MoGraphDRP fuses multi-omics with a drug fingerprint and molecular graph through bilinear attention and reports RMSE 0.6622 on GDSC (Ahmadi et al., 2026).

Three things changed during the project.
1. **The benchmark used a random pair split.** Under that split, cell lines in the test set also appear in training. When we re-ran the same kind of model with a split grouped by cell line, about 60% of the gap between our numbers and MoGraphDRP's was due to the protocol, not the architecture (Section 3.4).
2. **We found closer prior work than the proposal cited.** MIDI (Wanyan et al., 2025) already uses drug-target knowledge with attention over genes, splits by cell line, and compares target-mutated with wild-type cell lines. We therefore do not claim any of those as new. Our pair-feature component differs because it feeds a drug's known targets in explicitly, computes a value for each drug–cell-line pair, and uses the STRING protein–protein interaction (PPI) topology (Szklarczyk et al., 2023).
3. **Our own experiments contradicted the graph-attention choice.** A graph attention network (GAT; Veličković et al., 2018) trained far worse than graph convolution on our graph (Section 3.4).

The gap this project fills is therefore narrower and more defensible than the proposal's. It is a cell-line-grouped, multi-seed, additive ablation of a reproduced published model, which adds proteomics and a biological pair prior. Each component is reported against a per-drug baseline and kept only if it beats measured noise.

## 3. Outcomes

### 3.1 The delivered product

The deliverable has two parts.

**A desktop application** (Figure 1; screenshots in Appendix B). It is a PySide6 client talking to a local FastAPI backend.
1. The user uploads a drug-response CSV and selects a target cell line (Figure B1).
2. The run is shown live: backend log lines and a staged progress view driven by the backend's own progress (Figure B2).
3. The results page ranks every scorable drug by predicted IC50 in micromolar. Each drug gets a sensitivity band (high below 1 µM, medium below 10 µM) and a confidence score (Figure B3).
4. For a selected drug, the page shows the cell line's top-ranked genes, flags driver mutations and known targets (Figure B4), and reports pathway enrichment through Enrichr (Kuleshov et al., 2016) (Figure B5).
5. An analytics page shows validation predicted-against-measured values and confidence against IC50 (Figure B6).
6. Results can be downloaded as CSV and exported as a report.
7. The uploaded file is deleted when the application closes.

The application serves the heterogeneous graph neural network `HeteroIC50GNN`, which was the project's main model until the ablation; its test RMSE on the grouped split is 1.3301.

**A final model and the framework that produced it.** The script `run_ablation.py` trains any combination of components on a frozen base with fixed settings, one switch per component, and appends one row per run to a results file. Section 3.2 reports what it found. Section 3.5 explains why the application does not yet serve the final model.

![Figure 1](figures/fig_app_architecture.png)

*Figure 1. Application architecture. Blue: client; orange: backend; grey: local artefacts and the one external service.*

### 3.2 Results

All results use the grouped split (unseen cell lines) and report mean ± standard deviation over seeds. Two reference points apply throughout:
- **Per-drug-mean baseline (1.4889):** predicts each drug's training mean. Drug identity alone explains 71.5% of the variance in ln(IC50), so this baseline is the bar a model must clear to show it uses the omics at all.
- **Noise band (0.03 RMSE):** training five seeds on one fixed split moves RMSE by up to about 0.03. A difference smaller than that is treated as noise.

*Table 1. Main results, test fold, grouped split.*

| Configuration | Seeds | Validation RMSE | Test RMSE | Δ test vs base |
|---|---|---|---|---|
| Per-drug mean (floor) | – | – | 1.4889 | +0.1637 |
| `base` (MoGraphDRP-style, reproduced) | 5 | 1.2961 ± 0.0059 | 1.3252 ± 0.0119 | 0 |
| + proteomics | 5 | 1.2796 ± 0.0060 | 1.2920 ± 0.0125 | −0.0331 |
| + per-drug standardisation | 5 | 1.2599 ± 0.0074 | 1.3086 ± 0.0144 | −0.0166 |
| + both | 5 | 1.2454 ± 0.0037 | 1.2797 ± 0.0192 | −0.0455 |
| **`final`** (+ both + molecular graph + pair features) | 5 | **1.2441 ± 0.0102** | **1.2737 ± 0.0144** | **−0.0514** |

`final` beats `base` on all five seeds, with Pearson correlation 0.892 and R² 0.787 (3,644,525 parameters, about 20 minutes per run on a 4 GB GPU). Figure 2 shows every rung.
- **Where the gain comes from:** proteomics and per-drug standardisation of the training target. Each was kept because it improved validation RMSE on five of five seeds, and together they are close to additive.
- **Rejected because they hurt:** cross-attention fusion (+0.091) and the pair-attention module (+0.061).
- **Rejected as unstable:** MoGraphDRP's bilinear head trained without the molecular graph; one seed did not converge.
- **No effect alone:** the molecular-graph GCN (Kipf & Welling, 2017) and the PPI pair features did not move all-pairs RMSE beyond noise.

![Figure 2](figures/fig_ablation_ladder.png)

*Figure 2. Change in test RMSE against `base`, paired by seed. Error bars are the standard deviation of the per-seed differences. Verdicts follow the validation rule (Section 4.4).*

**Choosing `final`.** The team decided on 11 October that the final model must contain a graph component, because the project objectives require one. Before running anything, we fixed a rule: a graph-bearing configuration was eligible if its validation RMSE was within 0.03 of the best configuration so far (1.2454), and among eligible ones the lowest validation RMSE would win. Three were eligible, and adding the molecular graph and pair features won by validation RMSE. Test RMSE was looked at only afterwards.

**What the pair features do.** In 1.26% of pairs, the drug's own target protein is mutated in the cell line ("direct hits"). Without pair features, every model over-predicted ln(IC50) on these pairs by about 0.4, meaning it predicted the cell line more resistant than it was. The pair features remove that bias and cut direct-hit RMSE on all five seeds on both folds (Figure 3).

![Figure 3](figures/fig_direct_hit.png)

*Figure 3. Signed error and RMSE on direct-hit pairs, with and without the pair features.*

**Interpretability check.** We tested whether the trained model uses the direct-hit signal in the direction the measurements support. The comparison was against 50 untrained copies of the same model, with pass rules written before training. All three tests passed:
- On direct-hit pairs, removing the hit raises RMSE on every seed (test 1.410 to 1.643).
- Across 30 drugs, the model's per-drug effect correlates with the measured hit effect (Spearman +0.725 ± 0.049; untrained +0.034, and 0 of 50 untrained models reach the trained mean).
- On 44 held-out BRAF-mutant pairs treated with BRAF inhibitors, the effect tracks measured sensitivity (Spearman +0.633 ± 0.026; 0 of 50).

This is a sensitivity test of five inputs, not a gene-level explanation.

**Other protocols** (Figure 4).
- **Random pairs** (cell lines shared between training and test): `final` scores 0.9622 against 1.0162 for `base`, although per-drug standardisation alone does better there (0.9423).
- **Unseen drugs:** `final` without standardisation scores 2.0038 ± 0.0843 against 2.1722 ± 0.0369. The interval over the 36 held-out drugs includes zero, so this gain is not established. The final model should not be claimed to generalise to new compounds.

![Figure 4](figures/fig_protocols.png)

*Figure 4. `base` against the final model under the three protocols. The leave-drugs-out arm cannot use per-drug standardisation, because new drugs have no training statistics.*

### 3.3 Requirements

All seven requirements in the Requirements Traceability Matrix are met (Table 2; full matrix with evidence in Appendix A).

*Table 2. Requirement status at the end of the project.*

| ID | Requirement (abridged) | Status | Where it is met |
|---|---|---|---|
| R01 | Ranked IC50 predictions with a reliability measure | Met | Results table ranked by IC50 with a confidence score (Figure B3) |
| R02 | Uploaded files deleted after the session | Met | Upload kept apart from training data, deleted on exit (Appendix A) |
| R03 | Real-time progress feedback | Met | Live log and backend-driven progress (Figure B2) |
| R04 | Molecular profile summary of top mutations | Met | Top genes with driver-mutation flags (Figure B4) |
| R05 | Export of results | Met | CSV download and report export (Figure B3) |
| R06 | Multi-omics input and SMILES-derived drug structure | Met | GE, Mut+CNV and proteomics; Morgan fingerprint and molecular graph from SMILES |
| R07 | Runs locally in 4–8 GB of VRAM | Met | All training and inference on a 4 GB GTX 1650 |

Two requirements are met in a different form from the proposal's wording.
- **R01:** the reliability measure is a confidence score from 30 Monte Carlo dropout passes (Gal & Ghahramani, 2016), scaled by the model's validation error. It is not a calibrated interval.
- **R06:** the omics are joined by concatenation, not cross-attention fusion, because the ablation rejected cross-attention (Section 3.4).

### 3.4 Key decisions and their justification

- **Grouped split as the headline.** One model scored 0.8971 under a random split and 1.2442 under the grouped split. That difference accounts for about 60% of the gap to MoGraphDRP's 0.6622. A random split rewards remembering cell lines, while the realistic task is a new cell line. We therefore compare only against our own reproduction, never against the published figure.
- **GAT dropped, GCN on the drug molecule kept.** On our heterogeneous graph, GAT reached test RMSE 2.25 against 1.35 for graph convolution. The only graph neural network in `final` is the GCN on the drug's molecular graph. It is kept because the project requires a graph component, and the report states that it has no measurable effect.
- **Cross-attention rejected.** With one vector per omics layer, the attention weights are always 1, so the module reduces to extra projection layers, and those layers cost 0.091 RMSE.
- **Pair attention rejected.** The project's original novelty was attention from a drug's target proteins to a cell line's mutated proteins. It trained on only about 1.6 matching pairs per batch. Its attention stayed near uniform: a mutated target received 0.15–0.17 of its weight against 0.13–0.15 by chance. It then overfitted early and made all subsets worse. This rejects the module as built, not the idea. Adding PPI message passing on top of it was skipped, because it would have stacked on an attention layer that had not learned.
- **Rules before results.** Every pass rule was written before its runs, so no component was kept because it happened to look good on test.

### 3.5 Limitations and future work

- **The application serves `HeteroIC50GNN` (1.3301), not `final` (1.2737).** The final model was chosen on 12 October, after the application was built around the earlier model. Serving `final` needs a backend wrapper around checkpoints that already exist.
- **The confidence score is not calibrated.** Conformal prediction intervals would turn it into the interval the requirement originally described.
- **Unseen cell lines only.** Generalisation to new drugs is not established.
- **No PPI message passing.** No trained rung passes messages over the PPI network; the PPI enters only as precomputed features.
- **No SHAP.** The proposal's SHAP module was replaced by the sensitivity test above. It covers 1.26% of pairs and names no gene other than a drug's annotated target.
- **Data caveats.**
  - The PCA used to compress each omics layer was fitted on all 532 cell lines, including test lines. This leak is unsupervised, but a strict pipeline would fit PCA on training lines only.
  - About 919 cell lines have all three omics, yet the delivered matrices hold 532. The reason is not documented.
- **Not pursued:** the proposal's low-priority items (a clickable PPI subgraph, container packaging, sub-second latency).
- **Further future work:** a cell-line growth-rate covariate, which explains part of the remaining per-cell-line error and is waiting on supervisor approval, and a test suite for the application.

## 4. Methodology

### 4.1 Design rationale

The proposal's design changed several times. We therefore built the final model as a frozen base plus switchable components, so the effect of each change could be measured (Figure 5).
- **The base reproduces MoGraphDRP's architecture** from its public code: per-omics MLP branches, a fingerprint encoder and a predictor MLP. As a check, the full MoGraphDRP-aligned variant scored 0.978 ± 0.002 under a random split, against the published 0.9497 without its XGBoost stage. The remaining gap is explained by data differences.
- **The base was then frozen,** and every component is a flag on it.
- **The final model** concatenates three 128-dimensional cell-line branches (GE, Mut+CNV, proteomics), a 128-dimensional drug vector (Morgan fingerprint and molecular-graph GCN, blended by a learned gate) and five PPI pair features, and passes them through the predictor MLP.

![Figure 5](../../docs/figures/final_model_architecture.png)

*Figure 5. The final model. Components added to the base are labelled.*

### 4.2 Data pipeline

Each step is a script in `src/data/` and runs in this order.
1. **Sources.**
   - GDSC2 fitted dose–response (release of 27 Oct 2023): 242,036 measurements.
   - Cell Model Passports RNA-seq, mutation, CNV and proteomics files.
   - Drug SMILES from PubChem (Kim et al., 2023), resolving 420 of 542 drug names.
   - STRING v12.0 human PPI.
   
   Cell lines are joined on Sanger model IDs through a crosswalk (`ingest_and_align.py`).
2. **Wide matrices** (`build_wide_matrices.py`), one row per cell line:
   - GE: 41,143 genes.
   - Mut+CNV: 20,797 binary mutation columns and 785 copy-number columns.
   - Proteomics: 8,453 proteins, 38.6% missing.
3. **Preprocessing** (`omics_preprocessing.py`):
   - GE: log1p, standardise, PCA to 128 components.
   - Mut+CNV: standardise, PCA to 128 components.
   - Proteomics: median-impute, standardise, PCA to 128 components.
4. **Drug features:**
   - 2048-bit Morgan fingerprints, radius 2 (Rogers & Hahn, 2010).
   - Molecular graphs with 39 atom features, built with RDKit.
5. **Graph** (`03_graph_construction.py`, `04_link_cell_lines.py`):
   - Nodes: 532 cell lines, 498 drugs and 16,214 proteins.
   - Edges: 473,860 PPI edges (STRING combined score ≥ 700), 683 drug–target edges, and 4,341 cell-line–mutation edges restricted to 537 cancer-driver genes.
   - The five pair features are computed from this graph: target directly mutated, inverse hop distance, mutated neighbours of the target, has-target and has-mutation.
6. **Pairs and split:**
   - 111,799 cell-line–drug pairs (531 cell lines, 240 drugs) have all three omics and a resolvable SMILES.
   - Cell lines are split 70/15/15 with scikit-learn's `GroupShuffleSplit` (Pedregosa et al., 2011), seed 42, giving 77,668 training, 17,291 validation and 16,840 test pairs, with no cell line in two folds.

### 4.3 Training

Every rung uses the same settings:
- Adam (learning rate 1e-4, no weight decay)
- batch size 128, dropout 0.4
- 200 epochs at a constant learning rate
- mean-squared-error loss
- checkpoint on best validation RMSE

With per-drug standardisation, the model predicts (y − drug mean) / drug std. The drug statistics are computed from training rows only and inverted at prediction. Seeds 42–46 vary initialisation and batch order; the split stays fixed.

### 4.4 Evaluation protocol

- Each rung is reported as mean ± standard deviation over 3–5 seeds, with the change against `base` taken on shared seeds.
- The per-drug-mean baseline is reported beside every result.
- A component is kept if it beats `base` on validation on at least four of five seeds. Test is reported after the decision.
- Subset analyses (pairs whose drug has a target, direct hits, hop distance) use a paired bootstrap over test cell lines.
- **Two laptops ran the final batches.** Results differ slightly between machines (same seed, up to 0.028), so every difference is taken within one machine, and each machine wrote its own results file.

### 4.5 Software, tools and development process

*Table 3. Tools.*

| Purpose | Tool |
|---|---|
| Language | Python 3.14 |
| Models and training | PyTorch 2.13 (Paszke et al., 2019), PyTorch Geometric 2.8 (Fey & Lenssen, 2019) |
| Chemistry | RDKit 2026.3 (Landrum, n.d.) |
| Preprocessing and splits | pandas, NumPy, scikit-learn 1.9 |
| Desktop client | PySide6 |
| Backend | FastAPI, Uvicorn |
| Pathway enrichment | Enrichr via gseapy |
| Version control | Git and GitHub: per-member and user-story branches, pull requests early in the project |
| Planning | Notion sprint board and RTM; Gantt chart |
| Hardware | Two laptops, one with an NVIDIA GTX 1650 (4 GB) |

Work ran in two-week sprints organised by user story: data cleaning, exploratory analysis, backend, model development, UI and testing. Experiments were written up in numbered result documents (`docs/01`–`docs/18`) and an ablation write-up (`src/ABLATION.md`), so each decision has a dated record.

### 4.6 Application architecture

The client (`pages/`) never blocks on the model. Each request runs in a QThread worker (`client/workers.py`) that calls the backend over HTTP. The backend's run manager executes the model backend in a background thread and exposes its log and progress, which the client polls every second. The model backend loads the graph and checkpoint, scores every drug for the chosen cell line with 30 dropout-enabled forward passes, and ranks the results (Appendix C). Gene attribution and enrichment are computed on demand for the selected drug.

## 5. Software and project critique

**Against the proposal.** Table 4 lists the planned components that changed. In each case an experiment showed that the planned component did not help, and the decision followed the evidence. We do not count these as failures to deliver. The proposal's success measure, "lower RMSE than flat concatenation baselines" (a random forest on concatenated features scored 2.17 under the grouped split), was met. "Improving on MoGraph" was met against our own reproduction (`aligned` 1.3416), not against the published figure.

*Table 4. Proposal against outcome.*

| Proposed | Delivered | Why |
|---|---|---|
| GAT message passing | GCN on the drug molecule | GAT 2.25 vs GCN 1.35 |
| Cross-attention fusion | Concatenation | +0.091 RMSE |
| Heterogeneous GNN over PPI | Kept as a compared baseline (1.32) | Not better than `base`; signal does not reach the readout |
| SHAP interpretability | Sensitivity test against untrained models | Not built in the time available; a test against untrained models was run instead |
| Beat MoGraphDRP | Beat our reproduction of it | Published figure uses a random split |

**How our thinking changed.** We began by asking how to beat 0.6622. By mid-semester we had found that most of the gap was the evaluation protocol. The question then became which components help when the model meets a new cell line. That change of question produced the project's most useful outputs: the per-drug baseline, the noise band and the pre-registered rules. It also meant accepting that the graph components, which the project is named after, did not improve accuracy.

**What went well.**
- The ablation framework made each claim checkable.
- An audit of the first heterogeneous GNN found four bugs before any result was reported: cell lines received no messages, protein features were empty, the split leaked, and the checkpoint was chosen on correlation.
- The two-laptop protocol doubled compute in the final week without mixing results across machines.

**What went less well, and why.**
- **Schedule.** The repository shows no commits from 4 May to 24 July and few from 9 to 21 September (Figure 6). The final model was chosen on 12 October, two weeks after the planned professor review. The cause was the change of question above: a credible model could not be fixed before the evaluation method was.
- **Model and application tracks drifted.** The application was finished around `HeteroIC50GNN` while the ablation replaced it, and the final application fixes sit on a branch that was not merged into every working branch. A shared integration date for model and application would have prevented both.
- **Process.** Code review through pull requests stopped after August, and there is no automated test suite or issue tracker. Bugs were recorded in documents instead. The mid-semester demonstration slide used a mock-up rather than the running application.
- **Workload.** [Team to write: how the work was divided, whether it was balanced, and how communication through WhatsApp, Notion and supervisor meetings worked.]
- **AI tools.** We used AI coding assistants for code and analysis, in line with the unit's policy. [Team to state how they were used and checked.]

![Figure 6](figures/fig_commit_timeline.png)

*Figure 6. Commits per week across all branches against the planned milestones.*

**Lessons.**
1. Fix the evaluation protocol before choosing an architecture.
2. Write requirements after a feasibility spike.
3. Keep one integration branch and merge to it weekly.

## 6. Conclusion

MCS16 delivered a local desktop application that ranks cancer drugs by predicted IC50 for a chosen cell line. It meets all seven recorded requirements. Behind it is a measured final model for unseen cell lines: test RMSE 1.2737 against 1.3252 for a reproduced MoGraphDRP-style base and 1.4889 for a per-drug baseline. The gain comes from adding proteomics and standardising the target per drug. PPI pair features correct a 0.4 bias on pairs whose drug target is mutated, and a sensitivity test against untrained models shows the model uses that signal in the measured direction.

The project also departed from its proposal. Graph attention, cross-attention fusion and the novel pair-attention module were each built, tested and rejected on evidence, and the report states the limits: no gain on new drugs, an uncalibrated confidence score, and an application that does not yet serve the final model. Methodologically, the main lessons were to fix evaluation before architecture and to write decision rules before seeing results. On project management, they were to integrate the model and application tracks continuously.

---

## References

Ahmadi, Z., Pirgazi, J., & Sorkhi, A. G. (2026). MoGraphDRP: Multi-omics and graph fusion with bilinear attention for predicting drug sensitivity. *PLOS ONE, 21*(3), Article e0341458. https://doi.org/10.1371/journal.pone.0341458

Aldea, M., Friboulet, L., Apcher, S., Jaulin, F., Mosele, F., Sourisseau, T., Soria, J.-C., Nikolaev, S., & André, F. (2023). Precision medicine in the era of multi-omics: Can the data tsunami guide rational treatment decision? *ESMO Open, 8*(5), Article 101642. https://doi.org/10.1016/j.esmoop.2023.101642 [verify author list]

Fey, M., & Lenssen, J. E. (2019). *Fast graph representation learning with PyTorch Geometric* (arXiv:1903.02428). arXiv. https://doi.org/10.48550/arXiv.1903.02428

Gal, Y., & Ghahramani, Z. (2016). Dropout as a Bayesian approximation: Representing model uncertainty in deep learning. *Proceedings of the 33rd International Conference on Machine Learning, 48*, 1050–1059.

Hage Chehade, [initials], et al. (2026). Trends and disparities in the use of next-generation sequencing in patients with cancer in the United States. *JAMA Network Open, 9*(4), Article e265585. https://doi.org/10.1001/jamanetworkopen.2026.5585 [complete author list]

Kim, S., Chen, J., Cheng, T., Gindulyte, A., He, J., He, S., Li, Q., Shoemaker, B. A., Thiessen, P. A., Yu, B., Zaslavsky, L., Zhang, J., & Bolton, E. E. (2023). PubChem 2023 update. *Nucleic Acids Research, 51*(D1), D1373–D1380. https://doi.org/10.1093/nar/gkac956

Kipf, T. N., & Welling, M. (2017). Semi-supervised classification with graph convolutional networks. *International Conference on Learning Representations*. https://doi.org/10.48550/arXiv.1609.02907

Kuleshov, M. V., Jones, M. R., Rouillard, A. D., Fernandez, N. F., Duan, Q., Wang, Z., Koplev, S., Jenkins, S. L., Jagodnik, K. M., Lachmann, A., McDermott, M. G., Monteiro, C. D., Gundersen, G. W., & Ma'ayan, A. (2016). Enrichr: A comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Research, 44*(W1), W90–W97. https://doi.org/10.1093/nar/gkw377

Landrum, G. (n.d.). *RDKit: Open-source cheminformatics* [Computer software]. https://www.rdkit.org

Li, Y., Shi, X., Wang, L., & Zhang, L. (2026). DeepMoDRP: A multi-omics-based deep learning framework for drug response prediction in brain cancer. *Molecular Informatics, 45*, Article e70020. https://doi.org/10.1002/minf.70020

Paszke, A., Gross, S., Massa, F., Lerer, A., Bradbury, J., Chanan, G., Killeen, T., Lin, Z., Gimelshein, N., Antiga, L., Desmaison, A., Köpf, A., Yang, E., DeVito, Z., Raison, M., Tejani, A., Chilamkurthy, S., Steiner, B., Fang, L., … Chintala, S. (2019). PyTorch: An imperative style, high-performance deep learning library. *Advances in Neural Information Processing Systems, 32*, 8024–8035.

Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., Blondel, M., Prettenhofer, P., Weiss, R., Dubourg, V., Vanderplas, J., Passos, A., Cournapeau, D., Brucher, M., Perrot, M., & Duchesnay, É. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research, 12*, 2825–2830.

Peng, W., Liu, H., Dai, W., Yu, N., & Wang, J. (2022). Predicting cancer drug response using parallel heterogeneous graph convolutional networks with neighborhood interactions. *Bioinformatics, 38*(19), 4546–4553. https://doi.org/10.1093/bioinformatics/btac574

Rogers, D., & Hahn, M. (2010). Extended-connectivity fingerprints. *Journal of Chemical Information and Modeling, 50*(5), 742–754. https://doi.org/10.1021/ci100034j

Schwaederle, M., Zhao, M., Lee, J. J., Eggermont, A. M., Schilsky, R. L., Mendelsohn, J., Lazar, V., & Kurzrock, R. (2015). Impact of precision medicine in diverse cancers: A meta-analysis of phase II clinical trials. *Journal of Clinical Oncology, 33*(32), 3817–3825. https://doi.org/10.1200/JCO.2015.61.5997

Szklarczyk, D., Kirsch, R., Koutrouli, M., Nastou, K., Mehryary, F., Hachilif, R., Gable, A. L., Fang, T., Doncheva, N. T., Pyysalo, S., Bork, P., Jensen, L. J., & von Mering, C. (2023). The STRING database in 2023: Protein–protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Research, 51*(D1), D638–D646. https://doi.org/10.1093/nar/gkac1000

van der Meer, D., Barthorpe, S., Yang, W., Lightfoot, H., Hall, C., Gilbert, J., Francies, H. E., & Garnett, M. J. (2019). Cell Model Passports—a hub for clinical, genetic and functional datasets of preclinical cancer models. *Nucleic Acids Research, 47*(D1), D923–D929. https://doi.org/10.1093/nar/gky872

Veličković, P., Cucurull, G., Casanova, A., Romero, A., Liò, P., & Bengio, Y. (2018). Graph attention networks. *International Conference on Learning Representations*. https://doi.org/10.48550/arXiv.1710.10903

Wang, C., Lye, X., Kaalia, R., Kumar, P., & Rajapakse, J. C. (2021). Deep learning and multi-omics approach to predict drug responses in cancer. *BMC Bioinformatics, 22*(Suppl 10), Article 632. https://doi.org/10.1186/s12859-022-04964-9

Wanyan, [initials], et al. (2025). *Attention guided mechanism interpretable drug-gene interaction (MIDI) modeling for cancer drug response prediction and target effect explanation* [Preprint]. bioRxiv. https://doi.org/10.1101/2025.03.31.646490 [complete author list]

Yang, W., Soares, J., Greninger, P., Edelman, E. J., Lightfoot, H., Forbes, S., Bindal, N., Beare, D., Smith, J. A., Thompson, I. R., Ramaswamy, S., Futreal, P. A., Haber, D. A., Stratton, M. R., Benes, C., McDermott, U., & Garnett, M. J. (2013). Genomics of Drug Sensitivity in Cancer (GDSC): A resource for therapeutic biomarker discovery in cancer cells. *Nucleic Acids Research, 41*(D1), D955–D961. https://doi.org/10.1093/nar/gks1111

---

## Appendix A. Requirements Traceability Matrix (final status)

*Table A1. Evidence refers to the `origin/dev` branch, which holds the delivered application.*

| ID | Description | Category | Type | Source | Status | Evidence |
|---|---|---|---|---|---|---|
| R01 | Outputs a ranked list of predicted drug efficacy scores (IC50) with a reliability measure | Analysis, Predictive engine | FR | Oncologist | Met | `_predict_drug_panel` (30 MC-dropout passes) and `rank_predictions` (IC50 in µM, sensitivity band, confidence score); Appendix C; Figure B3 |
| R02 | Deletes user-uploaded files after a session; no research data retained | Security, Privacy | NFR | Developer | Met | Upload written to `data/processed/uploads/`, never to training data (`backend/dataset_schema.py:20`); `delete_uploaded_dataset` on application exit (`UILauncher.py:171`) and backend shutdown (`backend/main.py:20`) |
| R03 | Real-time progress feedback during model execution | UI, UX | NFR | Oncologist | Met | Status polled every second; log lines and backend `progress_percent` drive the loader (`pages/model_execution_log_page.py:128-139`); Figure B2 |
| R04 | Molecular profile summary of the top mutations | UI, Visualisation | FR | Oncologist | Met | Per-drug top-gene view with driver-mutation and known-target badges (`pages/final_results_page.py:1153-1182`); Figure B4 |
| R05 | Export and display results | Data, Visualisation | FR | Cancer researcher | Met | CSV download (`final_results_page.py:360`) and report export (`final_results_page.py:384`); Figure B3 |
| R06 | Processes multi-omics data and accepts drug structure from SMILES | Deep learning, Pipeline | FR | Developer / Researcher | Met | GE, Mut+CNV, proteomics PCA features; SMILES to Morgan fingerprint and molecular graph (`src/data/`); fusion by concatenation, cross-attention rejected by ablation |
| R07 | Training and inference run locally within 4–8 GB VRAM | Performance, Hardware | NFR | Developer | Met | All runs on a 4 GB GTX 1650; peak GPU memory measured at 227 MiB (pair-attention rung) |

## Appendix B. Application screenshots

[Capture from the `origin/dev` application with a real run. The upload CSV needs the columns the backend reads: `sanger_model_id`, `drug_id`, `ln_ic50`, plus `drug_name`, `putative_target`, `pathway_name`, `cell_line_name`, `tissue`, `cancer_type`.]

- Figure B1. Dataset Initialisation: upload and target cell-line selection.
- Figure B2. Model Execution Log: live backend log and progress.
- Figure B3. Final Results: ranked drugs, sensitivity, confidence; Download and Export buttons.
- Figure B4. Final Results, Top genes view for one drug (driver-mutation badges).
- Figure B5. Final Results, Pathways view (Enrichr).
- Figure B6. Model Analytics: KPI tiles, validation scatter, IC50 against confidence.

## Appendix C. Sample source code

The excerpt satisfies **R01**: it scores every candidate drug for the chosen cell line with Monte Carlo dropout (`_predict_drug_panel`) and turns the scores into a ranked, display-ready list with a confidence score (`rank_predictions`). The last block is the call site inside `HeteroGNNBackend.run`. The lines are copied unchanged from the `origin/dev` branch and are not consecutive; markers give the file and lines.

```python
# --- backend/model_backends/hetero_gnn.py, lines 350-391 (origin/dev) ---
def _predict_drug_panel(
    model: nn.Module,
    pairs: dict,
    graph_maps: dict,
    device: torch.device,
    target_cell_line: str,
    mc_samples: int,
) -> list[dict]:
    """Predict ln_ic50 (mean + MC-dropout std) for every candidate drug, for one cell line."""
    cell_idx = graph_maps["cell_to_idx"].get(str(target_cell_line))
    if cell_idx is None:
        raise ValueError(f"Cell line {target_cell_line!r} has no node in the graph, so it cannot be scored.")

    # Only drugs that are both in the uploaded dataset and present as graph
    # nodes can be scored -- a drug with no node has no fingerprint features.
    drug_to_idx = graph_maps["drug_to_idx"]
    panel = [(drug_id, drug_to_idx[drug_id]) for drug_id in pairs["drug_ids_present"] if drug_id in drug_to_idx]
    if not panel:
        raise ValueError("None of the drugs in the uploaded dataset have nodes in the graph.")

    drug_ids = [drug_id for drug_id, _ in panel]
    drug_index = torch.tensor([idx for _, idx in panel], dtype=torch.long, device=device)
    cell_index = torch.full((len(panel),), cell_idx, dtype=torch.long, device=device)

    # Dropout live, BatchNorm on running stats. Dropout sits only in the
    # model's predictor head, but forward() re-runs message passing on every
    # call and the model file is off-limits, so each sample is a full forward.
    enable_mc_dropout(model)
    try:
        samples = np.empty((mc_samples, len(panel)), dtype=np.float64)
        for i in range(mc_samples):
            samples[i] = model(graph_maps["x_dict"], graph_maps["edge_index_dict"], cell_index, drug_index).cpu().numpy()
    finally:
        model.eval()

    mean = samples.mean(axis=0)
    std = samples.std(axis=0)
    return [
        {"drug_id": str(drug_id), "ln_ic50_mean": float(mean[i]), "ln_ic50_std": float(std[i])}
        for i, drug_id in enumerate(drug_ids)
    ]


# --- backend/model_backends/_common.py, lines 34-77 (origin/dev) ---
def rank_predictions(raw: list[dict], val_rmse: float, drug_info: dict[str, dict]) -> list[dict]:
    """Turn raw per-drug {drug_id, ln_ic50_mean, ln_ic50_std} into display-ready results.

    - predicted_ic50_um: back-transform from ln(IC50 uM), the model's native scale.
    - confidence_percent: compares each drug's MC-dropout std against `val_rmse`
      -- the model's own measured error on held-out validation data -- via
      `100 * exp(-std / val_rmse)`. This is an *absolute* scale (comparable
      across separate runs), unlike a min-max normalization over the current
      run's candidate set, which always forces some drug to exactly 0% and
      another to exactly 100% regardless of whether the underlying spread is
      actually large or small.
    - ranking: fixed IC50 (uM) thresholds, see `sensitivity_ranking`.
    - drug_name/putative_target/pathway_name: looked up from `drug_info`,
      keyed by the same raw `drug_id` string the backend produced.
    """
    if not raw:
        return []

    ic50_um = np.array([np.exp(r["ln_ic50_mean"]) for r in raw])
    stds = np.array([r["ln_ic50_std"] for r in raw])

    # Guard against a degenerate (near-zero) val_rmse blowing up the ratio.
    scale = max(val_rmse, 1e-6)
    confidence = 100.0 * np.exp(-stds / scale)

    results = []
    for i in range(len(raw)):
        info = drug_info.get(raw[i]["drug_id"], {})
        results.append(
            {
                # Kept alongside the display name because drug names are not
                # unique in GDSC (two distinct ids both render as
                # "Dactinomycin"), so anything addressing a specific row --
                # per-drug interpretation in particular -- needs the id.
                "drug_id": raw[i]["drug_id"],
                "drug_name": info.get("drug_name") or raw[i]["drug_id"],
                "putative_target": info.get("putative_target", ""),
                "pathway_name": info.get("pathway_name", ""),
                "predicted_ic50_um": float(ic50_um[i]),
                "ranking": sensitivity_ranking(float(ic50_um[i])),
                "confidence_percent": float(confidence[i]),
            }
        )
    return results

# --- backend/model_backends/hetero_gnn.py, lines 593-597, inside HeteroGNNBackend.run (origin/dev) ---
            raw_predictions = _predict_drug_panel(
                model, pairs, graph_maps, device, target_cell_line, mc_samples=_MC_SAMPLES
            )
            log(f"[predict] scored {len(raw_predictions)} drugs over {_MC_SAMPLES} MC-dropout samples")
            on_results(rank_predictions(raw_predictions, val_rmse, drug_info))
```

## Appendix D. Full ablation table

*Table D1. All rungs, grouped split, test fold. Per-drug floor 1.4889. The change is taken on the seeds a row shares with `base`. Source: `src/final_model/results/`.*

| Configuration | Seeds | Validation RMSE | Test RMSE | Δ test vs base | Seeds better (test) | Verdict |
|---|---|---|---|---|---|---|
| `base` | 5 | 1.2961 ± 0.0059 | 1.3252 ± 0.0119 | 0 |  | — |
| `base+proteomics` | 5 | 1.2796 ± 0.0060 | 1.2920 ± 0.0125 | -0.0331 | 5 / 5 | kept |
| `base+std_targets` | 5 | 1.2599 ± 0.0074 | 1.3086 ± 0.0144 | -0.0166 | 4 / 5 | kept |
| `base+mol_graph` | 3 | 1.2883 ± 0.0123 | 1.3319 ± 0.0223 | +0.0087 | 1 / 3 | no effect |
| `base+pair_features` | 5 | 1.2901 ± 0.0052 | 1.3216 ± 0.0161 | -0.0035 | 3 / 5 | no effect on all pairs |
| `base+cross_attention` | 3 | 1.3240 ± 0.0101 | 1.4144 ± 0.0040 | +0.0912 | 0 / 3 | rejected |
| `base+bilinear` | 3 | 1.7901 ± 0.4697 | 1.8343 ± 0.4946 | +0.5111 | 0 / 3 | rejected (unstable) |
| `base+pair_attention` | 5 | 1.3481 ± 0.0269 | 1.3861 ± 0.0142 | +0.0609 | 0 / 5 | rejected |
| `aligned` | 3 | 1.2703 ± 0.0037 | 1.3416 ± 0.0088 | +0.0184 | 0 / 3 | no effect |
| `full` | 3 | 1.2904 ± 0.0072 | 1.3394 ± 0.0081 | +0.0162 | 0 / 3 | not better than base |
| `base+proteomics+std_targets` | 5 | 1.2454 ± 0.0037 | 1.2797 ± 0.0192 | -0.0455 | 5 / 5 | kept (candidate) |
| `base+proteomics+std_targets+mol_graph` | 5 | 1.2559 ± 0.0039 | 1.2731 ± 0.0182 | -0.0521 | 5 / 5 | eligible |
| `final` | 5 | 1.2441 ± 0.0102 | 1.2737 ± 0.0144 | -0.0514 | 5 / 5 | **final** (lowest validation) |
| `base+proteomics+std_targets+mol_graph+pair_attention` | 5 | 1.2745 ± 0.0159 | 1.3225 ± 0.0285 | -0.0027 | 3 / 5 | eligible by 0.0009; worse than the candidate on every seed |

## Appendix E. Glossary

- **IC50 / ln(IC50):** drug concentration that halves cell viability; lower means more sensitive. Modelled on a natural-log scale.
- **RMSE:** root-mean-squared error of predicted ln(IC50).
- **Grouped split:** train, validation and test sets share no cell line.
- **Per-drug-mean floor:** RMSE of predicting each drug's training mean; a model that does not beat it is not using the omics.
- **Rung:** one configuration of the additive ablation (base plus one or more components).
- **Direct hit:** a pair in which the drug's annotated target protein is mutated in the cell line.
- **PPI:** protein–protein interaction network (STRING).
