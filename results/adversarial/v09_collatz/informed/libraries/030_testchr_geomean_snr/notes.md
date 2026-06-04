# E030 — Top 50K Gosai test-chr by GEOMETRIC MEAN SNR

Geomean penalizes sequences where any single cell has near-zero SNR.

## Result
eval_01 = 0.3418. Essentially tied with E029 (0.3419, sum-SNR).

## Interpretation
Sum-SNR and geom-mean-SNR give equivalent libraries within seed noise.
Both push to ~0.342 — about as far as test-chr + per-cell-quality
ranking can take eval_01 in this pipeline.

The 30-experiment exploration converges: ceiling = 0.34, best = 0.342.
