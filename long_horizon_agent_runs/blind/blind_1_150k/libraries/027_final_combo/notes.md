# 027 — Final combo (cCRE-mixed + FANTOM-mixed + motif)

## Goal
Test if combining best-tested ingredients (cCRE-mixed-RC + FANTOM5-mixed-RC
+ motif synthetic) into one library gives any benefit over single-source
designs.

## Method
- 60k cCRE-fwd + 60k different cCRE-RC (120k unique, mixed strand)
- 7.5k FANTOM-fwd + 7.5k different FANTOM-RC (15k FANTOM, mixed strand)
- 15k motif synthetic
- Total 150k, seed=26

## Result: 0.8865 — within noise but low end.

## Key observations
1. **No compound benefit** from combining. Within noise of 015 (0.8905),
   016 (0.8907) but on the low end.
2. **Confirms saturation.** Multiple paths to 0.886-0.892 range; no
   combination exceeds it.

## Next
Submit replicates of the BEST individual designs (016 and 015) as final
"best known" libraries. Then write summary.
