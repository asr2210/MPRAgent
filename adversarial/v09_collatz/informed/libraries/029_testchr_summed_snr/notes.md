# E029 — Top 50K Gosai test-chr by SUMMED SNR (across all 3 cells)

Combined score = SNR_K562 + SNR_HepG2 + SNR_SKNSH (sum, not max).

## Result — NEW BEST
eval_01 = 0.3419. Per-cell K562=0.172, HepG2=0.222, SKNSH=0.632.

vs E022 max-SNR: 0.3410 → +0.001
vs E008 random chr: 0.3359 → +0.006

## Interpretation
Summed SNR slightly beats max SNR. The sum-based ranking selects
sequences with consistently good signal in multiple cells, rather than
outliers loud in just one cell. The model gets richer cross-cell signal
per sequence.

Within-seed noise: ±0.001-0.002, so +0.001 over E022 is borderline.
But the per-cell improvements are consistent (all three cells equal or
slightly above E022 levels). Suggests a real (small) effect.

## Final best library: E029
chr 7/9/13/21/X + top 50K by summed per-cell SNR.
