# 017 — Dirichlet(0.3) with mid-GC rejection

## What I built
Sampled Dirichlet(0.3) compositions, kept only those with GC ∈ [0.3, 0.7].
Acceptance rate 29.4%. 50k accepted compositions, one sequence each.

## Result
- eval_01 = 0.0773 (vs Dirichlet(0.3) 0.0786). -0.0013.
- mean ≈ 0.0951 (vs 0.0954). approximately neutral.
- eval_04, eval_09 went UP (0.0902 → 0.0924 each).

## Interpretation
Concentrating Dirichlet(0.3) compositions to mid-GC range HURT eval_01. The extreme-GC
sequences in the natural Dirichlet(0.3) distribution were apparently CONTRIBUTING to
eval_01 performance, not just diluting it.

Both directions (014: more uniform GC, 017: more concentrated GC) hurt vs unaltered
Dirichlet(0.3). **Natural Dirichlet(0.3) GC distribution is the precise sweet spot
for eval_01.**

But eval_04/09 prefer mid-GC concentration. **Different evals have different preferred
composition distributions.** This is an irreducible tradeoff.

## What to try next
The Dirichlet ceiling is firmly at ~0.078 eval_01 for any single-strategy library.
Time to test a NEW source: random hg38 regions (not regulatory DHS). Maybe non-DHS
genomic sequences provide different compositional coverage that adds info.

Or test the alpha extreme: alpha=0.05 completes the sweep.
