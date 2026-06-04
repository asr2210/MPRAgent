# 005 — Dirichlet(0.1) very-extreme composition

## What I built
50k sequences with per-sequence base composition ~ Dirichlet(0.1)^4. Even more concentrated near corners (single-base dominance) than exp 004.

## Result
- eval_01 = 0.0752 — WORSE than Dirichlet(0.3) (0.0786) but better than DHS-stratified (0.0739).
- mean across evals: 0.0935 (vs 0.0954 in exp 004).

## Interpretation
Non-monotonic in alpha. Sweet spot around 0.3 (or maybe 0.2-0.5). Going more extreme (toward homopolymers) hurts — consistent with homopolymer_rich baseline (0.0570).

Composition variance lever has DIMINISHING RETURNS. Likely we're near its plateau (~0.078-0.079 on eval_01 from composition alone).

## Theory updates
- T1.2 holds, but the sweet spot is INTERIOR (alpha ~ 0.3), not at the limit.
- To break above ~0.08 on eval_01, I likely need a DIFFERENT axis of learnable variance, not more composition.

## Sequence of dirichlet alpha results:
- alpha=0.1 → 0.0752
- alpha=0.3 → 0.0786 (best so far)
- alpha=1.0 → 0.0768 (strategies.md baseline)

So a slight peak at alpha ~ 0.3. Pure composition ceiling appears to be ~0.079.

## What to try next
Need to test a different axis. Most informative next:
(a) DHS biology PLUS Dirichlet(0.3): does biology add value on top of the winning composition? Or is composition all that matters?
(b) Sequences explicitly maximizing dinucleotide diversity (16-D vs 4-D Dirichlet) — captures CpG, TpA, etc. variation.
(c) Test mixes of Dirichlet at different alphas (some pathological, some uniform) — wider variance regime.

(a) is the most direct test of "is biology useful at all in this harness" — important before abandoning biology entirely.
