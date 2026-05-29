# 023 — Kitchen sink (all winning levers combined)

## Goal
Test if stacking all winning ingredients (cCRE-RC bulk + FANTOM5-RC +
multi-source synthetic) gives any compound benefit above noise.

## Method
- 60k cCRE-fwd + 60k different cCRE-RC (120k cCRE, mixed strand)
- 7.5k FANTOM5-fwd + 7.5k different FANTOM5-RC (15k, mixed strand)
- 10k motif-embedded + 5k uniform random (15k synthetic, 2 sources)
- Total 150k.

## Result: 0.8886. Within noise of 0.889 cluster.

| design | mean |
|--------|------|
| 015 (cCRE + motif) | 0.8905 |
| 022 (015 replicate, different seed) | 0.8875 |
| 016 (cCRE + FANTOM + motif) | 0.8907 |
| 023 (kitchen sink) | 0.8886 |

Spread is 0.8875–0.8907 = 0.003 range, which is exactly the measured
noise floor. So all four are statistically equivalent at 0.889 ± 0.003.

## Key observations
1. **No compound benefit from stacking levers.** Combining cCRE-RC,
   FANTOM-RC, and 2-source synthetic gives the same result as plain
   cCRE+motif.
2. **The design space within "cCRE-derived bulk + synthetic" is
   saturated.** Multiple paths lead to 0.889 ± 0.003.
3. **To break the ceiling, need genuinely new data/lever:**
   - Cell-type-specific cCRE selection (needs per-cell-type signal data)
   - Conservation filter (needs phyloP/phastCons)
   - Multi-organism conserved sequences
   - Activity-stratified selection (would require pilot)

## Next
Given saturation, remaining experiments should focus on:
- One more replicate of best recipe for confidence
- Variations that might give >0.005 effects
- Final library submission

Going with experiment 024 = test if MULTI-SOURCE synthetic (3 different
diversity sources in equal small chunks) helps. Last variation in
synthetic component design.
