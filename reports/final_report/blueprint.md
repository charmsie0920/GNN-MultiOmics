# Final report: blueprint

The report is due **2 Nov 2026**. The word limit is 3,500, counting Sections 1–6 only; figures, tables, references and appendices are excluded. It must be drafted in an **edit-logged Google Doc**, and the authorship statement is mandatory (a missing statement costs 50%). Each member also clicks Submit on Moodle.

Files in this folder:
- `draft.md`: a complete first draft, currently about 3,280 words, leaving about 220 for the two team-written paragraphs in Section 5. Paste it into Google Docs with *Edit → Paste from Markdown*, then each owner rewrites their sections in their own words. The edit log is the evidence of authorship.
- `figures/`: Figures 1–4 and 6 as PNGs. Rebuild them with `.venv/bin/python reports/final_report/figures/make_report_figures.py`. Figure 5 is `docs/figures/final_model_architecture.png`, relabelled on 27 Oct for cross-attention fusion (edited in the `.svg` and exported with Chromium; the `.excalidraw` source still shows concatenation).

## Marks and where the draft spends words

| Section | Marks | Draft words | Owner (proposed) | What earns the marks (rubric) |
|---|---|---|---|---|
| 1 Introduction | 10 | ~270 | SM | Project, team, problem, outcome as achieved, a roadmap of the report |
| 2 Background | 10 | ~420 | WH (2.1), CC (2.2) | Two paragraphs; literature **updated since the proposal** (MIDI, the split-protocol finding); the gap we fill; competing tools |
| 3 Outcomes | 30 | ~1,480 | SM (3.1, 3.3), AA (3.2, 3.4), CC (3.5) | All results tied to requirements, the RTM, justified decisions, honest limits and future work |
| 4 Methodology | 15 | ~690 | CC (4.1), AA (4.2–4.4), SM (4.5–4.6) | Tools and versions, data steps with counts, development steps, replicable |
| 5 Critique | 15 | ~440 + team text | WH | Causes of limitations, project management (schedule, communication, workload, supervisor), lessons |
| 6 Conclusion | 5 | ~180 | WH | Recap of product, results, method, project management |
| Sample code (App. C) | 5 | – | SM | 50–100 lines, documented, a function and its call, satisfies an RTM item. Done: R01, 79 lines verbatim from `main` |
| References | 5 | – | CC | APA 7, cited items only |
| Style | 5 | – | CC (proof editor) | Numbered headings matching the TOC, numbered figures referenced in the text, front cover complete |

## Sources for every number

Do not quote numbers from anywhere else.

- **Ablation, final model, graph decision, protocols, Phase 6:** `CLAUDE.md` and `src/ABLATION.md`. The raw rows are in `src/final_model/results/*.csv`. All 60 decimals in the draft were checked verbatim against these files.
- **Data pipeline counts:** `docs/data_ingestion_and_alignment.md`, `docs/graph_construction.md`, `src/data/*`.
- **Early experiments:**
  - split protocol: `docs/09`
  - GAT vs GCN: `docs/07`
  - noise band: `docs/13`
  - hetero-GNN bug audit: `docs/hetero_gnn_test_bugfixes.md`
  - MIDI: `docs/17_midi_review.md`
- **Background statistics:** `reports/Precision oncology citation verification.md`, which also lists the slide citations that need correcting.
- **RTM evidence:** the `main` branch. The delivered app, serving the final model since `fbf49b3` (26 Oct), was merged there in `68372f8`. The RTM and Appendix C line numbers refer to `main`.
- **Final model since 27 Oct: `final_xattn`** (`final` with cross-attention fusion), `src/ABLATION.md` §5 update of 2026-10-27 and §8: validation 1.2652, test 1.3180 (5 seeds). `final` (concatenation, test 1.2737) is reported as the model it replaced.
- **Served model:** a five-seed ensemble of `final_xattn`, `src/ABLATION.md` §8: validation 1.2217, test 1.2764.

## Do not reuse

- **The "Our Model" slide of the final presentation (`docs/FINAL-3.pdf`, slide 10).** It labels pair attention between drug targets and mutated proteins as "Cross-Attention Fusion". That module was rejected (Section 3.4); the cross-attention in the final model is between omics layers, with every weight equal to one. Use Figure 5. The same applies to the older FINAL.pdf architecture slides (pair attention, bilinear head, GraphSAGE).
- **Headline numbers from the presentation's closing slide** (1.27 vs 1.33 grouped, 0.96 vs 1.02 random). They are the concatenation model's; the final model is 1.3180 vs 1.3252. Also not "Every part improves prediction" (slide 7).
- **The "MVP DEMO" slide image.** It is a mock-up with invented metrics (Pearson 0.842, an AUC axis). Use real screenshots (Appendix B).
- **`docs/architecture.md`.** It describes `HeteroIC50GNN`, and its data-source table names DepMap files the code does not use.
- **"0.66 vs 1.33".** The documented comparison is MoGraphDRP 0.6622 against our 1.2442 grouped and 0.8971 random (doc 09).
- **"41 days"** from the problem-statement slide. It is unsupported.
- **Experiment numbers (E01, E11, …).** They differ between docs, so name experiments instead.

## Claims to keep exactly as worded

These come from `CLAUDE.md`, "What to say and not say".
- The final model beats `base` on validation (−0.031, 5 of 5) but not on test (−0.007, 2 of 5). Do not say it beats the concatenation baseline on test.
- Cross-attention was adopted because it met a validation margin written before the run; it costs 0.044 test RMSE against concatenation. Every attention weight is 1.0, so no gain or behaviour may be credited to attention. Do not give "novelty" as the reason.
- The concatenation model's gain over `base` comes from proteomics and standardisation. The molecular graph and pair features do not move all-pairs RMSE beyond noise.
- In the final model the pair features only partly remove the direct-hit bias (+0.26 on test, `base` +0.40); its direct-hit RMSE is no better than without them.
- The only GNN in `final` is the GCN on the drug molecule. Do not call `final` a PPI GNN or a heterogeneous GNN.
- Interpretability is a sensitivity test of five inputs on 1.26% of pairs against untrained models. It is not SHAP and not a gene-level explanation.
- On unseen drugs the gain is not established. Do not claim generalisation to new compounds.
- Compare against our MoGraphDRP reproduction, never against the published 0.6622.
- The app serves `final_xattn` as a five-seed ensemble. Its gene view names only a drug's targets and a cell line's driver mutations, so pathway enrichment runs on about 4% of pairs.

## Team to-do before submission

1. **Confirm or change the authorship statement** in the draft. It is a proposal based on each member's work in git.
2. **Take screenshots B1–B6** from the `main` app with a real run. The upload needs `sanger_model_id, drug_id, ln_ic50`; add `drug_name, putative_target, pathway_name, cell_line_name, tissue, cancer_type` for readable names. For B4–B5, use KCL-22 (`SIDM01004`) with Olaparib (`1017`).
3. **Write the two placeholder paragraphs** in Section 5: workload and communication, and use of AI tools. The draft deliberately leaves both to the team.
4. **Complete the two references** marked `[verify]` or `[complete author list]`: Aldea 2023 and Wanyan 2025 (MIDI). Hage Chehade 2026 was filled from the presentation's reference slide. Check every DOI resolves.
5. **Insert the Google Docs TOC** and the final word count on the cover. Rename figure files if the doc needs them embedded rather than linked.
6. **Optional, if words allow:** one sentence on the growth-rate covariate decision once the supervisor answers.
