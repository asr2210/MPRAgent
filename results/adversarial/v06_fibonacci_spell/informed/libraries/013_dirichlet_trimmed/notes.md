# Experiment 013 — Trimmed Dirichlet(0.5)

## Result
eval_01 = **0.1379** (vs Dirichlet(0.5) 0.1395 — slightly worse)

Acceptance rate: 33.5% (kept compositions where min(p) >= 0.03).
Trimmed distribution: GC std=0.222 (vs original 0.286), max-base
mean=0.527 (vs original 0.616).

Per-cell:
- K562 = 0.0424 (vs 002 = 0.0473) — LOWER, not higher
- HepG2 = 0.1699 (vs 002 = 0.1697) — tied
- SK-N-SH = 0.2014 (vs 002 = 0.2015) — tied

## Interpretation
Trimming hurts K562 specifically, while leaving HepG2/SK-N-SH unchanged.
The EXTREME tail of Dirichlet(0.5) (where one base dominates >70%)
is actually beneficial to the K562 head, not harmful.

This is the opposite of what I expected. The near-homopolymer
baseline (50k pure homopolymer) was 0.058, but a few thousand
near-homopolymers MIXED with normal Dirichlet samples actually
HELPS K562.

## Lesson
Dirichlet(0.5)'s tail is informative; don't trim it. Aggressive
extremes work IN MIXTURE with moderate compositions. Pure homopolymers
fail because there's nothing else to contrast with.

## Next
Measure single-seed noise floor: re-run Dirichlet(0.5) with seed=99.
If 0.1395 ± 0.005, plateau is real and I should focus remaining
experiments on diagnostics. If 0.1395 is just lucky seed, more
experiments could break through.
