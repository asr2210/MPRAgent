# 030 — Diversity selection v2 (4-pool, 200k candidates)

## Plan
Refine 029 by adding a 4th candidate pool (008-style: 5 motifs
no-overlap, random bg) for a richer candidate distribution. Same
greedy 6-mer selection of 50k. Hypothesis: more candidate diversity
→ better selection → higher mean_r.

## Result
**eval_01 = 0.4269.** K562 0.597, HepG2 0.627, SK-N-SH 0.056.
14-eval avg = 0.4251. **Below 029's 0.4288.**

Pool composition picked: A=14138, B=8600, **C=4622**, D=22640.
Pool D dominated because sparse-motif sequences have the largest
undisturbed-backbone region → highest unique-6-mer count. But this
PUSHED OUT the dinuc-bg contributors (C dropped 7696 → 4622),
exactly the seqs SK-N-SH needs.

## What this teaches
- **Diversity-selection is sensitive to candidate-pool composition.**
  Naive "more pools = better" fails. The selector greedily maximizes
  the local objective (k-mer novelty) and can systematically
  under-pick the slices that carry per-cell-type signal.
- **The 029 3-pool composition was lucky/optimal**: random-bg dense
  seqs have decent diversity AND carry motif signal; dinuc-bg has
  enough rare 6-mers to still get selected at ~15% rate.
- **SK-N-SH signal is fragile under diversity-selection.** Going from
  3 to 4 pools dropped SK-N-SH from 0.0653 → 0.0563 (−0.009).
- The 008-style sparse-motif design (which dominated selection here)
  is actually a WORSE training signal than the dense designs — it
  has high k-mer count but low motif content. The selector got
  fooled by pure backbone novelty.

## Theory T21 (refinement)
Greedy k-mer diversity isn't a magic objective — it's optimal only
when candidate pools all carry similar PER-SEQUENCE signal density.
Adding a high-diversity but low-signal pool causes selection to
drift toward that pool, reducing total signal in the selected library.

For optimal diversity-selection libraries:
- All candidate pools should have COMPARABLE k-mer entropy per seq
- Otherwise the selector over-picks the high-entropy ones regardless
  of signal content
- Solution: weight selection by signal density, OR limit pool sizes
  proportional to desired final composition, OR use a more nuanced
  score that balances novelty against signal proxies

## Final library decision
**029 (eval_01 = 0.4288) is the best library produced** in this 30-exp
run, beating the prior plateau of 0.4283-0.4284 (008/020). The exact
configuration was:
  - 150k candidates: 50k @ 020-design + 50k @ 022-design + 50k @ 026-design
  - Greedy 6-mer diversity selection of 50k
  - Final composition: A=24832, B=17472, C=7696
  - eval_01 = 0.4288, 14-eval avg = 0.4265
  - SK-N-SH means 0.0653, with three evals above 0.07
