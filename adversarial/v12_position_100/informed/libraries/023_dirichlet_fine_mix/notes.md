# 023 — Fine Dirichlet alpha mix {0.2, 0.3, 0.4}

## What I built
~16.7k sequences each at alpha ∈ {0.2, 0.3, 0.4} = 50k total.

## Result
- eval_01 = 0.0779 (vs Dirichlet(0.3) 0.0786). -0.0007 loss.
- eval_08 = 0.0728 (vs 0.0716). +0.0012 gain.
- mean ≈ 0.0952 (vs 0.0954). Essentially tied.

## Interpretation
Fine alpha mix is approximately neutral. Confirms peak is sharp at exactly 0.3, and
mixing neighboring alphas dilutes slightly. eval_08 prefers slightly higher alpha, so
mixing helps it; eval_01 prefers exactly 0.3, so mixing hurts it slightly.

## What to try next
Test if SMALL DOSE of biology can supplement without diluting too much. Exp 024:
40k Dirichlet(0.3) + 10k random hg38.
