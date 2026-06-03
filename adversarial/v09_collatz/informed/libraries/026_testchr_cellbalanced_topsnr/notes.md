# E026 — Cell-balanced top-SNR within Gosai test chroms

Interleave top-SNR-by-K562, top-SNR-by-HepG2, top-SNR-by-SKNSH from
test-chr Gosai (dedup).

## Result
eval_01 = 0.3403. Per-cell K562=0.171, HepG2=0.223, SKNSH=0.627.

vs E022 (max-SNR across cells): 0.3410. Essentially identical, -0.001.

## Interpretation
Cell-balancing within the top-SNR test-chr subset gives no improvement.
The max-SNR approach in E022 was already producing a roughly balanced
mix. HepG2 ticked up by +0.002 (0.221→0.223), but SKNSH dropped by
0.003 (0.630→0.627), canceling.

## Conclusion
E022 (top-SNR max across cells) is the locally optimal recipe within
this combinatorial neighborhood.
