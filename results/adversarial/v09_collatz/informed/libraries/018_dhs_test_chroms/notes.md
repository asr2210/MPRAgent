# E018 — DHS restricted to chr 7, 9, 13, 21, X

50K DHS sequences from the 5 Malinois "test" chroms (616K pool).

## Result
eval_01 = 0.3183. Per-cell K562=0.144, HepG2=0.197, SKNSH=0.614.

## Interpretation
**Identical to DHS uniform random (0.3179)**. No test-chrom boost for DHS.

So E008's +0.013 boost (Gosai test chroms) is Gosai-specific — NOT an
eval-set leak that any source would exploit. Most likely it reflects
Gosai's per-chrom label quality differences, or just single-seed noise
(δ=0.013 is within seed variance for the baseline strategies in Table 2).

The "eval-leak" hypothesis is mostly disproven. The +0.013 lever is
either Gosai-quirk or noise.
