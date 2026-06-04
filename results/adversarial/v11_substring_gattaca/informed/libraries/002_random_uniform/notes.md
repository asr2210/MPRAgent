# E2: random_uniform (calibration)

## Design
- 50,000 IID random 200bp sequences. Uniform over A/C/G/T.
- Seed 0. No additional structure.

## Result (single seed)
- eval_01 = **0.8565**  (published 5-seed baseline: 0.8566 — match!)
- mean across 14 evals ≈ **0.841**
- Published random_uniform 5-seed mean across 14 evals: 0.840

**Pipeline is calibrated.** Single-seed variance is tiny — published
5-seed means are within 0.001 of my single-seed run.

## Implications
- Baseline ceiling for "pure random" approaches is ~0.857 on eval_01 / ~0.841 averaged.
- The published gc_50 baseline (0.8591, exactly-50%GC sequences) sits ~0.003 above.
- My single-seed E2 effectively reproduces random_uniform → I trust deltas of ~0.005+ as meaningful.

## Reference for future experiments
| eval | E2 (my random_uniform) |
|------|------------------------|
| 01 | 0.8565 |
| 02 | 0.8565 |
| 03 | 0.8514 |
| 04 | 0.8604 |
| 05 | 0.8565 |
| 06 | 0.8585 |
| 07 | 0.8007 (hard) |
| 08 | 0.7723 (hardest) |
| 09 | 0.8604 |
| 10 | 0.8058 (hard) |
| 11 | 0.8585 |
| 12 | 0.8514 |
| 13 | 0.8264 (hard-ish) |
| 14 | 0.8565 |

Evals 01,02,05,14 are identical (likely same eval, different scoring).
Evals 04,09 are identical. 06,11 are identical. 03,12 are identical.
Distinct evals: 01-group, 03-group, 04-group, 06-group, 07, 08, 10, 13.
So functionally ~8 distinct evals.

Evals 07, 08, 10, 13 are "harder" — random_uniform scores 0.77-0.83 on these.
These probably contain structured/biological sequences that pure random can't fit perfectly.
