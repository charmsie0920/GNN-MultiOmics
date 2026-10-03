Pattern 1: every fusion result loses to proteomics alone. E04 (MLP, proteomics, fingerprint, 1.2843) beats every cross-attention configuration, every graph model, and the full architecture. More modalities have consistently made things worse.

Pattern 2: your cross-attention isn't doing what it claims. Each modality is a single pooled vector, so every attention call has one query and one key — softmax returns 1.0 by construction. Your own module docstring says so. The capacity comes from the projections and FFN, not from attention.

Those two facts are probably the same fact. Here's what I'd do about it.

Tier 1 — do these first
1. Seed sweep. Not an improvement, but nothing else is trustworthy without it. Your full-architecture run had val RMSE 1.2705 vs the no-PPI model's 1.2679 — 0.003 apart — yet their test scores differ by 0.038. That gap opening only at test time is exactly what unstable single-seed results look like. 3–5 seeds, mean ± std, for your top three configs. Small code change, decides whether your other findings are real.

2. Make the attention real: multi-token modalities. Instead of one 128-dim vector per modality, split it into k tokens (e.g. 4 × 32) so attention has something to distribute over. Then softmax is a genuine learned weighting rather than a no-op.

This is the highest-value architectural change available to you. It's cheap — a reshape in MultiOmicsCrossAttentionFusion plus pooling after — it makes "cross-attention fusion" an honest claim in your report, and it directly targets Pattern 1: real attention can learn to downweight a modality that's hurting, which the current degenerate version cannot.

3. Per-drug target standardization. ln(IC50) distributions differ enormously between drugs — some compounds are potent across the board, others inert. Standardize targets per drug using training-set statistics, predict in that space, then invert before scoring so RMSE stays comparable.

This is often a large win in drug-response prediction and it's a few lines. No leakage under your grouped split, since every drug appears in training. One caveat: it can't be used for leave-drugs-out, where unseen drugs have no statistics — so keep it as a separate arm, not a global change.

Tier 2 — benchmark alignment, good for the report
4. Bilinear attention head (MoGraphDRP §2.3). Replace concat → MLP with their multi-head bilinear module. ~50 lines, self-contained, isolates the component they call their key innovation. Your lecturer's "follow the benchmark's architecture" instruction points straight at this.

5. Fingerprint + molecular graph together. Worth noting you've been treating these as alternatives, but MoGraphDRP uses both — three fingerprint types and a graph GCN, concatenated. Your infrastructure already produces both representations on the same 111,799 pairs, so this is mostly plumbing. Their Table 4 claims the representations are complementary; you can test that claim directly.

6. COSMIC gene filtering instead of blanket PCA. They reduce each omics branch to ~600–700 cancer-relevant genes rather than PCA over everything. Given that proteomics alone beats every fusion, feature selection may matter more than fusion architecture. Bigger job — new preprocessing — but it's the most plausible explanation for the residual 0.235 gap you measured in 09_split_protocol_comparison.md.

Tier 3 — cheap experiments, run while waiting
Huber loss instead of MSE. Biological IC50 data is noisy with outliers; one flag, might buy a little.
Modality dropout during training — randomly zero a modality per batch. Directly attacks Pattern 1 by forcing robustness to any single modality dominating.
Seed ensembling. Averaging 5 seeds' predictions is a near-free RMSE improvement and is honest as long as you report it as an ensemble.
Light hyperparameter search. You've never tuned anything — LR, dropout and hidden dims were all inherited from the first script written. Tune on validation only.
What I'd skip
More GAT work. Three independent failures plus the benchmark's own Table 3. Write it up as a documented deviation and move on.

More PPI graph work, until Tier 1 lands. It currently costs 0.038 RMSE, and that branch has to climb back before any refinement to it is measurable.

AMOGEL's ARM subsystem. Still a semester of work for a branch that's behind.