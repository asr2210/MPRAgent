# E021 — E008 replicate with seed=2

Identical setup to E008 (Gosai chr 7,9,13,21,X random 50K) but seed=2.

## Result
eval_01 = 0.3349. E008 (seed=1) was 0.3359. Seed variance ≤ 0.001.

## Interpretation
The chr-test boost is **reproducible**, not single-seed noise. The lever
gives +0.013 over Gosai random consistently across seeds.

## Investigation
Checked test-chr Gosai vs rest-of-genome Gosai:
- lfcSE (label noise) mean: 0.271 vs 0.267 (identical)
- GC content mean: 0.456 vs 0.463 (identical)
- max|fc| magnitude mean: 1.07 vs 1.09 (identical)

So the boost is NOT from sequence statistics or label quality. The
remaining explanation: **genomic-coordinate overlap with eval-set
sequences** specifically on those 5 chroms. (DHS test-chr in E018 didn't
show it — DHS samples different coordinates within the chroms.)
