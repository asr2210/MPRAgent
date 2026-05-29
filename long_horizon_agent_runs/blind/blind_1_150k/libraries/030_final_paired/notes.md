# 030 — FINAL submission (012 recipe, paired-RC + 15k motif, seed=29)

## Goal
Final library submission. Uses the cleanest single-source winning recipe
(012-style paired RC + motif), with a fresh seed.

## Method
- 67,500 unique cCRE-centered windows (forward), seed=29
- 67,500 reverse complements of those same 67,500 cCREs (paired RC)
- 15,000 motif-embedded synthetic (2 JASPAR consensus motifs/seq)
- Total 150k.

## Result: 0.8883. Within noise of the cluster.

| recipe family | replicates | mean across | best individual |
|---------------|-----------|-------------|-----------------|
| 015 (cCRE-mixed + motif) | 015, 022, 029 | 0.8889 | 0.8905 |
| 012 (cCRE-paired-RC + motif) | 012, 030 | 0.8894 | 0.8905 |
| 016 (+ FANTOM half) | 016, 028 | 0.8889 | 0.8907 |

## Pick for final
The three best recipes all give the same population mean (≈ 0.889) within
the ±0.003 noise floor. Choosing 012-style for final because it is the
single cleanest "pure cCRE + RC + 15k motif" recipe.
