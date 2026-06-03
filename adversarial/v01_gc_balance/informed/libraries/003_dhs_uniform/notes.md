# 003_dhs_uniform — notes

**Design**: 50k 200bp windows centered on DHS summit, sampled uniformly without replacement from the 3.59M Meuleman 2020 DHS index. N-containing windows dropped.

**Result**: eval_01 = 0.6604, mean across 14 = 0.5786. 41s wall, 12s eval.

**Comparison**:
| eval | random | motif-planted | dhs_uniform | baseline dhs_random |
|---|---|---|---|---|
| 01 | 0.470 | 0.505 | 0.660 | 0.7089 |
| 04 | 0.395 | 0.177 | 0.493 | 0.7429 |
| 07 | 0.519 | 0.717 | 0.761 | 0.7615 |
| 08 | 0.158 | 0.093 | 0.116 | 0.6673 |
| 13 | 0.502 | 0.701 | 0.752 | 0.7639 |

**Two key observations**:
1. Reproduction gap: my dhs_uniform eval_01 = 0.66 vs baseline dhs_random = 0.71. ~0.05 below. Possible causes:
   - Single seed vs baseline 5-seed average
   - Sampling strategy: I use summit ± 100; baseline may use full DHS region with random 200bp window
   - N-filtering: I drop N-containing windows
2. **eval_08 is a special beast**. My dhs_uniform got 0.116 here — *worse than random uniform (0.158)*. Baseline dhs_random reports 0.667; that gap (0.55!) is far larger than any other eval. This suggests eval_08 may not even be a stable, single-seed-comparable signal: it might be very sensitive to the specific sequences chosen, or it may depend on inclusion of specific reference sequences.

**Theory update**: H1 still consistent — DHS provides genuine grammar (motif-planted-synth eval_07/13 improvements roughly match dhs_uniform), but doesn't fully cover the eval_08 regime.

**Next**: test H1 directly by combining DHS + motif-planted synthetic. If H1 (multi-objective spanning) is right, the hybrid should beat dhs_uniform.
