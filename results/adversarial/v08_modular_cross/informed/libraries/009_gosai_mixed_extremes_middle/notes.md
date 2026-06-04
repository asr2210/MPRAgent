# Experiment 009 — Mixed extremes + stratified middle

## Method
- 25K activity extremes (12.5K bottom + 12.5K top, lfcSE<0.3)
- 25K middle quintile-stratified (within middle 60% of activity, 5 bins x 5K)

## Result
**eval_01 = 0.0091** (vs plateau 0.018; pure extremes 0.014).
eval_04/09 = 0.0184 (vs pure extremes 0.022).

WORSE on both axes than either pure strategy.

## Interpretation
Combination is sub-additive — each component dilutes the other. Mixing the
"shape" of two libraries doesn't produce a model that satisfies both objectives.
Probably because:
- Extremes alone teaches the model to predict strong activity differences
- Stratified middle teaches graded response
- Combined, the model splits capacity between these regimes and is worse at
  both than a specialist

## Theory update
A single 50K library encodes ONE training distribution. The model learns the
mapping for that distribution. To match an eval set, the library needs to match
that eval's distribution — not a weighted average of multiple distributions.

For eval_01 the best library is "broad stratified Gosai with quality filter."
For eval_04/09 the best library is "Gosai activity extremes."
These are different libraries; cannot combine.

Next direction: try sub-source selection. Gosai = 14K CRE + 446K GTEX + 338K
UKBB. Sub-source may matter for eval composition.
