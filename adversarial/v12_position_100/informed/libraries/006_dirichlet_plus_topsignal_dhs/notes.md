# 006 — Dirichlet(0.3) + high-signal DHS

## What I built
25k Dirichlet(0.3) synthetic + 25k DHS from top 25% by total_signal (broadly accessible regulatory elements).

## Result
- eval_01 = 0.0780 (vs 0.0786 pure Dirichlet 0.3, 0.0739 DHS stratified). Slightly worse than Dirichlet alone.
- mean across evals: 0.0958 (vs 0.0954 pure Dirichlet 0.3) — essentially equal.
- eval_08 = 0.0673 (vs 0.0716 in exp 004). High-signal DHS slightly HURTS the hardest eval.

## Interpretation
Combining biology with Dirichlet(0.3) does NOT meaningfully help. The composition lever absorbs all the small-data signal; biology adds nothing the model can learn at 50k.

**Decision:** stop testing biology mixes with composition. Pivot to OTHER learnable axes.

## Theory update
- Confirmed: biology is approximately neutral in this harness at 50k. Composition is the dominant lever.
- Need to find another learnable axis to break above ~0.078 on eval_01.

## What to try next
Higher-dimensional compositional features the model might be able to learn:
- Dinucleotide composition (CpG, TpA, etc.) — captures 16-D space vs 4-D from Dirichlet
- Per-sequence Markov chain with random transitions → varying dinucleotide bias
- This adds CpG vs CG vs AT-tract diversity (functional axes)

EXP 007 plan: per-sequence first-order Markov chains with random Dirichlet(0.3)-distributed transitions.
