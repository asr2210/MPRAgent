# 026 — Chaotic 5-way mix

## What I built
10k each from: Dirichlet(0.3), Dirichlet(0.1), shuffled DHS, hg38 random,
Markov-Dirichlet(0.3). Total 50k.

## Result
- eval_01 = 0.0784 (vs Dirichlet(0.3) seed=0 0.0786). **Within seed noise of best.**
- mean ≈ 0.0950.

## Interpretation
5-way chaotic mix matched the pure Dirichlet(0.3) ceiling within noise. This is the
FIRST mixture that didn't underperform pure Dirichlet(0.3) by a clear margin.

Key insight: as long as ~20%+ of sequences are extreme-composition Dirichlet, the
model gets enough composition signal. The other 80% (lower-composition variants, bio,
hg38) doesn't add nor subtract significantly — within the ceiling.

This suggests there's a **robust plateau** of ~0.078-0.079 for ANY library with sufficient
extreme-composition coverage.

## What to try next
Test if multi-seed Dirichlet(0.3) ensemble (16.7k each from seeds 0, 1, 2) gives
better uniform coverage and possibly slightly higher result (or reduced seed variance).
