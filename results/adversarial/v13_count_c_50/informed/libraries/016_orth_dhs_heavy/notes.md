# Experiment 016 — orth-DHS heavy (35k orth-DHS + 15k cCRE)

## Design
Push 015's orthogonal-DHS direction further: 35 k DHS-non-cCRE +
15 k cCRE class-balanced.

## Result
| eval | 015 (25/25) | **016 (35/15)** | Δ |
|------|------------:|----------------:|--:|
| 01   | **0.5736** | 0.5689 | −0.005 |
| 07   | 0.6131 | **0.6250** | **+0.012** (new best) |
| 13   | 0.5924 | **0.6058** | **+0.013** (new best) |
| 04/09| 0.5535 | 0.5279 | −0.026 |
| 08   | 0.1492 | 0.0972 | −0.052 |
| mean | 0.560  | 0.555  | −0.005 |

## Interpretation
**The 25/25 split (015) is near-optimal for eval_01.** Pushing
toward more orth-DHS lifts eval_07/eval_13 (which reward broad
regulatory diversity) but loses eval_01 and crashes eval_04/eval_08.

**Different evals continue to want different priors**:
- eval_01: 25 k orth-DHS + 25 k cCRE balanced (015)
- eval_07/13: more orth-DHS, less cCRE (016)
- eval_04/09: more cCRE (012)
- eval_08: cCRE-PLS/pELS heavy (009 only)

Strong K562 drop (0.61 → 0.60 on eval_01) when going orth-DHS heavy.
This further confirms K562 is well-covered by cCRE; removing cCRE
budget loses K562 directly.

## Numbers
mean_r: 0.555
eval_01: 0.5689
eval_07: 0.6250 (best ever)
eval_13: 0.6058 (best ever)
eval_08: 0.0972
time_s: 23

## Next
- 25k orth-DHS COMPONENT-TARGETED + 25k cCRE (combine 011's
  component-weighting with 015's orthogonality)
- And try 20/30 (cCRE-heavy with orth-DHS) for the opposite direction.
