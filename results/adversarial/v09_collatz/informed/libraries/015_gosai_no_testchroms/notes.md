# E015 — Gosai EXCLUDING test chroms (no chr 7, 9, 13, 21, X)

## Result
eval_01 = 0.3191. Per-cell K562=0.146, HepG2=0.204, SKNSH=0.608.

## Interpretation
Excluding the 5 "test" chroms drops the score by only -0.004 vs random
Gosai. So the test-chrom boost (E008 +0.013) is real but small — not
catastrophic when removed.

Symmetric picture: those chroms ADD a little, removing them costs a little.
The pipeline isn't dominated by any single chrom-leak. Most signal is
chrom-agnostic.

## Conclusion
The test-chrom lever (≤+0.013) is the BIGGEST library-side effect found.
The 0.34 ceiling in this pipeline is unlikely to be breakable by any
single-source library re-weighting. Time to try fundamentally different sources.
