# Experiment 012 — Topic-uniform SynthSeqs

## Result
eval_01 = **0.1314** (vs exp 001 topic-weighted 0.1319 — within noise)

Per-cell:
- K562 = 0.0311 (vs 001: 0.0289 — slight bump)
- HepG2 = 0.1688 (vs 001: 0.1712 — slight drop)
- SK-N-SH = 0.1941 (vs 001: 0.1957 — slight drop)

## Interpretation
Topic balance shifts losses slightly between cell types but doesn't
meaningfully change overall correlation. The K562 lift (+0.002) is
within single-seed noise (~0.005).

Both synthseqs variants underperform Dirichlet(0.5) by ~0.008.
Real-DHS sequence structure doesn't help the evaluator regardless of
how topic representation is balanced.

## Implication for K562 strategy
Cell-type-targeted real sequences don't lift K562 meaningfully. The
K562 head's poor correlation (~0.04) appears to be a property of the
evaluator's K562 model itself, not something library composition can
easily fix via topic balancing.

## Next
Pivot back to Dirichlet (the winning regime). Test if trimming the
near-homopolymer tail (where any base < 3%) improves Dirichlet(0.5)
by removing the unproductive extreme compositions.

**Exp 013**: trimmed Dirichlet(0.5) — reject samples where min base
proportion < 0.03. Keeps most of the distribution (~90%) but removes
the very-skewed extreme tail.
