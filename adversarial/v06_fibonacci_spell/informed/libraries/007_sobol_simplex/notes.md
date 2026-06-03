# Experiment 007 — Sobol-uniform simplex coverage

## What I did
Sample 50,000 base compositions using a Sobol sequence on (0,1)^4
mapped to the 4-simplex via -log(u) normalization. This gives
quasi-uniform coverage of the composition simplex. Then iid 200bp per
composition.

## Result
- eval_01 = **0.1371** (vs dirichlet(0.5) 0.1395 — WORSE by 0.002)

## Interpretation
Uniform simplex coverage is *worse* than random Dirichlet(0.5). Why?
Dirichlet(0.5) biases compositions toward simplex EDGES (one or two
bases dominant). Sobol-uniform gives more compositions near the
center (0.25, 0.25, 0.25, 0.25). The model gets more signal from
EXTREME compositions than from balanced ones.

This is consistent with v06 baseline: random_uniform (all ~25% each)
scores 0.1158, lower than any composition-diverse strategy. The model
learns from variation in composition; balanced compositions don't
provide that.

## Sweet-spot finding
- alpha=0.5 → 0.1395 (best so far)
- mixed (0.3..2.0) → 0.1377
- Sobol (uniform on simplex) → 0.1371

Pure dirichlet(0.5) hits the optimum for this evaluator.

## Next
Test the "composition-matched random" idea: take real DHS sequences,
keep ONLY their per-sequence base composition, regenerate as iid.
This decomposes the real-DHS underperformance (0.1319 vs 0.1395)
into "limited composition diversity" vs "structural patterns hurt".
