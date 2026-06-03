# Experiment 017 — orth-DHS component-targeted

## Design
Combine 015's orthogonal-DHS pool with 011's component-targeted
sampling (Neural, Cancer/epi, Digestive heavy). 25 k + 25 k cCRE.

## Result
| eval | 015 (uniform orth) | **017 (CT orth)** | Δ |
|------|-------------------:|------------------:|--:|
| 01   | **0.5736** | 0.5696 | −0.004 |
| 07   | 0.6131 | 0.6246 | +0.012 |
| 13   | 0.5924 | 0.6058 | +0.013 |
| 04/09| 0.5535 | 0.5271 | −0.026 |
| 08   | 0.1492 | 0.0942 | −0.055 |

## Interpretation
**Component-targeting does NOT compose with orthogonality on eval_01.**
In the full DHS pool, component-targeting helped (011: +0.002 over
008). In the orthogonal pool, it hurts (-0.004 vs 015).

**Why**: The orthogonal pool (DHS-non-cCRE) is already enriched for
non-K562 cell types — K562-relevant elements live in cCRE. So the
orthogonal pool implicitly does what 011's component targeting tried
to do. Stacking both over-corrects and trades diversity.

Diagnostic: orthogonal pool component distribution
  Primitive/embryonic: 233k (biggest, NOT in 011's targets)
  Neural: 124k
  Lymphoid: 89k (K562-relevant, NOT in 011's targets)
  Cancer/epi: 66k

Uniform sampling from the orth pool gives natural diversity weighted
toward "everything except cCRE-covered" — which is the right thing.

eval_07/eval_13 *did* lift further (+0.012/+0.013) — they continue
to reward cell-type-diversifying priors regardless of source.

## Theory update
- **Orthogonality and component-targeting are not independent levers.**
  They address the same problem (HepG2/SK-N-SH coverage) via
  different mechanisms; combining them is double-dipping.
- **Uniform sampling from orth-DHS captures the cell-type lift
  automatically.** Better than explicit component weighting.

## Numbers
mean_r: 0.555
eval_01: 0.5696
eval_07: 0.6246
eval_13: 0.6058
time_s: 20

## Next
Try the OTHER ratio direction: 15 k orth-DHS + 35 k cCRE class-bal.
And test pure 50k orth-DHS to isolate the orth effect.
