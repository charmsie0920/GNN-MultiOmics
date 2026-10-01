# MIDI Review — Does It Already Do Our Novelty Component?

**Paper:** Wanyan et al., *Attention Guided Mechanism Interpretable Drug-Gene Interaction (MIDI) Modeling for Cancer Drug Response Prediction and Target Effect Explanation*, bioRxiv 2025.03.31.646490 (posted 2025-04-02, not peer reviewed)
**Date:** 2026-10-01 · **Why:** Phase 0 task in [`CLAUDE.md`](../CLAUDE.md)

## The question

Our planned novelty component is pair-specific cross-attention between a drug's target proteins and a cell line's mutated proteins over the STRING PPI graph. If MIDI's cross-attention already puts target proteins on the drug side, that claim would have to shrink to the evaluation alone.

**Answer: no, it does not.** The component is still ours, but the claim must be worded more narrowly, and two things we planned to present as new are already in MIDI.

## What MIDI does

| Part | MIDI |
|---|---|
| Drug input | Atoms, bond types and shortest-path distances from the SMILES, through a graph transformer. No target proteins. |
| Cell-line input | Per gene (6,144 genes): Geneformer identity embedding + expression embedding + mutation embedding, summed. |
| Cross-attention | The flattened drug embedding attends over the gene **identity** embeddings only (their eq. 11), with a temperature of 9. |
| Where targets enter | A supervised contrastive loss pulls a drug's embedding towards its known target genes' embeddings (eq. 16). Targets are a training signal, not an input. |
| Use of the attention | The scores scale each gene's embedding (eq. 12), which is then projected to a number and added to a drug-only number (eq. 13-14). |
| Graph prior | None. No PPI network. |
| Data | CCLE, 24 drugs, 471 cell lines, predicting activity area. |
| Split | By cell line, 80/20 within each cancer type (376 train, 95 test). |

The attention weights depend on the drug alone. The paper says so directly: each drug has one gene ranking, because its targets are fixed regardless of expression or mutation. The cell line only enters afterwards, through the embeddings being scaled.

## What is still ours

- **Target proteins as explicit tokens on the drug side.** MIDI learns a soft ranking of genes from drug structure; we feed the known targets in.
- **Pair-specific attention.** Our weights differ for every drug and cell line pair. MIDI's are per drug.
- **PPI topology** connecting mutated proteins to targets.
- **Untrained-model nulls** for the interpretability tests. MIDI compares against another model (TCRP), not against an untrained copy of itself.

## What we can no longer claim

- **Being first to use drug-target knowledge with attention over genes.** MIDI does this and claims to be first. Cite it as the closest prior work.
- **Cell-line-grouped splitting as a novelty.** MIDI also splits by cell line. It remains our headline protocol, but it is not new.
- **Testing predictions on target-mutated versus wild-type cell lines.** Their Fig. 6 does this for 19 drug-gene pairs (for example Nilotinib and ABL1, PLX-4720 and BRAF). Our direction-sign and BRAF tests in [experiment 15](./15_interpretability_validation_results.md) overlap with it.

## A limitation this exposes

MIDI works for any drug that has a structure. Our component needs known targets, and 130 of our drugs have no target edge. State this in the report, and report metrics separately for has-target and no-target pairs (already planned in Phases 3 and 5).

## Wording to use

> Pair-specific attention between a drug's known target proteins and a cell line's mutated proteins over the PPI graph, in contrast to MIDI's per-drug gene ranking learned from molecular structure.

## Caveats on the paper

- Preprint, not peer reviewed.
- The section "Adding Prior Drug-Gene Interaction Knowledge" is empty, and equations 9 and 10 are missing.
- The supplementary material was not read.
- Their scale (24 drugs) is far smaller than GDSC2, so their numbers are not comparable with ours.
