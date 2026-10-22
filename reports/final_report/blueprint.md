# Final report: blueprint

The report is due **2 Nov 2026**. The word limit is 3,500, counting Sections 1–6 only; figures, tables, references and appendices are excluded. It must be drafted in an **edit-logged Google Doc**, and the authorship statement is mandatory (a missing statement costs 50%). Each member also clicks Submit on Moodle.

Files in this folder:
- `draft.md`: a complete first draft, currently about 3,330 words. Paste it into Google Docs with *Edit → Paste from Markdown*, then each owner rewrites their sections in their own words. The edit log is the evidence of authorship.
- `figures/`: Figures 1–4 and 6 as PNGs. Rebuild them with `.venv/bin/python reports/final_report/figures/make_report_figures.py`. Figure 5 is `docs/figures/final_model_architecture.png`.

## Marks and where the draft spends words

| Section | Marks | Draft words | Owner (proposed) | What earns the marks (rubric) |
|---|---|---|---|---|
| 1 Introduction | 10 | ~270 | SM | Project, team, problem, outcome as achieved, a roadmap of the report |
| 2 Background | 10 | ~420 | WH (2.1), CC (2.2) | Two paragraphs; literature **updated since the proposal** (MIDI, the split-protocol finding); the gap we fill; competing tools |
| 3 Outcomes | 30 | ~1,330 | SM (3.1, 3.3), AA (3.2, 3.4), CC (3.5) | All results tied to requirements, the RTM, justified decisions, honest limits and future work |
| 4 Methodology | 15 | ~690 | CC (4.1), AA (4.2–4.4), SM (4.5–4.6) | Tools and versions, data steps with counts, development steps, replicable |
| 5 Critique | 15 | ~440 + team text | WH | Causes of limitations, project management (schedule, communication, workload, supervisor), lessons |
| 6 Conclusion | 5 | ~180 | WH | Recap of product, results, method, project management |
| Sample code (App. C) | 5 | – | SM | 50–100 lines, documented, a function and its call, satisfies an RTM item. Done: R01, 96 lines verbatim from `origin/dev` |
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
- **RTM evidence:** the **`origin/dev`** branch, the delivered app (commits `ea519cf`, `62b301a`, 13 Oct). Those commits are **not** on `aditya`, `main` or local `dev`, and the RTM line numbers refer to `origin/dev`.

## Do not reuse

- **FINAL.pdf "Our Model" and "Proposed final architecture" slides.** They show pair attention, a bilinear head and GraphSAGE over the PPI, none of which is in `final`. Use Figure 5.
- **The "MVP DEMO" slide image.** It is a mock-up with invented metrics (Pearson 0.842, an AUC axis). Use real screenshots (Appendix B).
- **`docs/architecture.md`.** It describes `HeteroIC50GNN`, and its data-source table names DepMap files the code does not use.
- **"0.66 vs 1.33".** The documented comparison is MoGraphDRP 0.6622 against our 1.2442 grouped and 0.8971 random (doc 09).
- **"41 days"** from the problem-statement slide. It is unsupported.
- **Experiment numbers (E01, E11, …).** They differ between docs, so name experiments instead.

## Claims to keep exactly as worded

These come from `CLAUDE.md`, "What to say and not say".
- The gain over `base` comes from proteomics and standardisation. The molecular graph and pair features do not move all-pairs RMSE beyond noise.
- The only GNN in `final` is the GCN on the drug molecule. Do not call `final` a PPI GNN or a heterogeneous GNN.
- Interpretability is a sensitivity test of five inputs on 1.26% of pairs against untrained models. It is not SHAP and not a gene-level explanation.
- On unseen drugs the gain is not established. Do not claim generalisation to new compounds.
- Compare against our MoGraphDRP reproduction, never against the published 0.6622.
- The app serves `HeteroIC50GNN` (1.3301), not `final`.

## Team to-do before submission

1. **Confirm or change the authorship statement** in the draft. It is a proposal based on each member's work in git.
2. **Take screenshots B1–B6** from the `origin/dev` app with a real run. The upload CSV needs these 9 columns: `sanger_model_id, drug_id, ln_ic50, drug_name, putative_target, pathway_name, cell_line_name, tissue, cancer_type`.
3. **Write the two placeholder paragraphs** in Section 5: workload and communication, and use of AI tools. The draft deliberately leaves both to the team.
4. **Complete the three references** marked `[verify]` or `[complete author list]`: Aldea 2023, Hage Chehade 2026, Wanyan 2025 (MIDI). Check every DOI resolves.
5. **Insert the Google Docs TOC** and the final word count on the cover. Rename figure files if the doc needs them embedded rather than linked.
6. **Optional, if words allow:** one sentence on the growth-rate covariate decision once the supervisor answers.
