# 017 — dELS-only cCREs (winning recipe)

## Goal
Test whether restricting to the dominant cCRE class (dELS, 74% of natural
distribution) lifts mean above 0.890.

## Method
- 67,500 dELS-only cCREs (forward), seed=16
- 67,500 different dELS-only cCREs (reverse complement), seed=16
- 15,000 motif-embedded synthetic, seed=16
- Total 150k.

dELS pool (class = `dELS` or `dELS,CTCF-bound`): 536,243 after filter.

## Result: BIG regression. Mean = 0.8753 (vs 015 0.8905, **−0.015**).

| eval | 015 | 017 | Δ |
|------|------|------|---|
| 01 | 0.8330 | 0.8187 | −0.014 |
| 02 | 0.9303 | 0.9171 | −0.013 |
| 03 | 0.9230 | 0.9080 | −0.015 |
| 04 | 0.8657 | 0.8499 | −0.016 |
| 05 | 0.8328 | 0.8186 | −0.014 |
| 06 | 0.9307 | 0.9174 | −0.013 |
| 07 | 0.9063 | 0.8929 | −0.013 |
| 08 | 0.9115 | 0.8868 | **−0.025** |
| 09 | 0.9467 | 0.9238 | **−0.023** |
| 10 | 0.9313 | 0.9185 | −0.013 |
| 11 | 0.8187 | 0.8049 | −0.014 |
| 12 | 0.8011 | 0.7856 | −0.015 |
| 13 | 0.9049 | 0.8941 | −0.011 |
| 14 | 0.9308 | 0.9174 | −0.013 |

## Key observations
1. **Both extremes hurt.** Class-balanced (007, −0.004) hurt; dELS-only
   (017, −0.015) hurts WORSE. The natural cCRE class distribution is the
   sweet spot.
2. **Every eval lost** by 0.011–0.025. The biggest losses are on eval_08
   (−0.025) and eval_09 (−0.023). Other cCRE classes (pELS, PLS,
   CTCF-only) are NOT redundant — they contribute information the dELS
   pool lacks.
3. **Hard evals 11/12 BOTH lost** (−0.014, −0.015). So they're NOT
   specifically about distal enhancers — they care about regulatory
   class diversity too.

## Theory update (v16 → v17)
- The natural cCRE class distribution (74% dELS, 13% pELS, 4% PLS, 9%
  other) is OPTIMAL. Both over-balancing (007) and over-focusing (017)
  hurt.
- Each class contributes unique regulatory grammar:
  - dELS: bulk distal enhancer grammar
  - pELS: proximal enhancer (TSS-near) grammar
  - PLS: promoter grammar (TATA, CCAAT, CpG)
  - CTCF: insulator boundary grammar
- Removing or over-representing any class breaks the model's
  generalisation.
- "Diversity of biological annotation type" is a real lever — at least
  *within* the cCRE catalogue.

## Next
The natural distribution is optimal, so filter-by-class is dead. Next:
test if STRATIFICATION by a DIFFERENT property (GC content) helps or
hurts. If the model wants natural-distribution-of-GC, GC-stratification
will hurt like class-stratification did. If GC is an axis the natural
distribution doesn't sample well, GC-stratification could help.

EXPERIMENT 018 = cCRE selection stratified by GC content (force uniform
GC distribution) + 015 recipe.
