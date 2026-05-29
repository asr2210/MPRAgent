# 018 — GC-stratified cCRE selection

## Goal
Test if forcing uniform GC distribution (vs natural peaked distribution
~45-50% GC) helps generalisation.

## Method
- Compute GC for each of 745k cCRE windows
- Bin into 30 GC bins (1-1.5% each)
- Sample 4500 per bin (with replacement for sparse bins)
- 67.5k for fwd, 67.5k for RC of different
- 15k motif-embedded synthetic

Natural GC distribution: very peaked, mode at GC ~50%, range 15-94%.
Extreme bins (GC <23%, >88%) have <50 cCREs each and required heavy
replacement sampling.

## Result: regression. Mean = 0.8791 (vs 015 0.8905, **−0.011**).

| eval | 015 | 018 | Δ |
|------|------|------|---|
| 01 | 0.8330 | 0.8219 | −0.011 |
| 02 | 0.9303 | 0.9202 | −0.010 |
| 03 | 0.9230 | 0.9118 | −0.011 |
| 04 | 0.8657 | 0.8629 | −0.003 |
| 05 | 0.8328 | 0.8216 | −0.011 |
| 06 | 0.9307 | 0.9205 | −0.010 |
| 07 | 0.9063 | 0.8897 | −0.017 |
| 08 | 0.9115 | 0.8985 | −0.013 |
| 09 | 0.9467 | 0.9423 | −0.004 |
| 10 | 0.9313 | 0.9152 | −0.016 |
| 11 | 0.8187 | 0.8074 | −0.011 |
| 12 | 0.8011 | 0.7883 | −0.013 |
| 13 | 0.9049 | 0.8870 | −0.018 |
| 14 | 0.9308 | 0.9208 | −0.010 |

## Key observations
1. **Regression on every eval** (~−0.01 each). Forcing uniform GC hurts
   the same way forcing equal classes (007) and dELS-only (017) did.
2. **Worst losses on eval_13 (−0.018) and eval_07 (−0.017)**. These
   evals seem most sensitive to having natural cCRE distribution.
3. **Sample-with-replacement in sparse GC bins (very low or very high
   GC) added duplicates**, which contributes to the loss — over-representing
   rare regions effectively decreases unique-cCRE coverage.

## Theory update (v17 → v18)
**Three independent strong negative results now confirm the same finding:**
- 007 class-balanced: −0.004 vs natural
- 017 dELS-only:      −0.015 vs natural
- 018 GC-stratified:  −0.011 vs natural
The natural cCRE distribution is OPTIMAL across multiple axes. Don't
restrict, don't rebalance, don't stratify — the model wants the unbiased
sample from the cCRE catalogue.

This is one of the most robust findings of the project so far.

## Next
Time to try a different lever entirely: **engineered synthetic with
designed regulatory architecture**. Specifically, homotypic motif
clusters (3 copies of the same JASPAR motif spaced 30bp apart, mimicking
real enhancer cooperativity). If structured synthetic > random-placement
synthetic, we have a new lever to push the synthetic component.

EXPERIMENT 019 = 015 base + 15k homotypic-cluster synthetic (3 copies of
same motif spaced 30bp apart) instead of 015's random-placement
synthetic.
