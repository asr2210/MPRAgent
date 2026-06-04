# Experiment 014 — ccre_heavy + celltype-targeted DHS (mosaic)

## Design
35 k cCRE class-balanced (eval_04/09 lifter, from 012)
+ 15 k DHS component-targeted (eval_07/13 lifter, from 011).

Hypothesis: priors compose — combining cCRE-heavy + CT-DHS should
lift eval_01 if effects are additive.

## Result
| eval | 008 (25/25) | 011 (CT, 25/25) | 012 (15/35) | **014** |
|------|------------:|----------------:|------------:|--------:|
| 01   | 0.5671 | 0.5688 | 0.5678 | 0.5664 |
| 04/09| 0.5723 | 0.5579 | **0.5772** | 0.5677 |
| 07   | 0.5833 | **0.5989** | 0.5828 | 0.5881 |
| 13   | 0.5603 | **0.5769** | 0.5589 | 0.5647 |
| 08   | 0.2111 | 0.1813 | 0.2169 | 0.2008 |
| mean | 0.554  | 0.555 | 0.554  | 0.554 |

## Interpretation
**Mosaic strategy gives partial credit on each direction**:
- eval_04/09: 0.5677 (between 011's 0.5579 and 012's 0.5772 — closer to 012)
- eval_07: 0.5881 (between 008's 0.5833 and 011's 0.5989 — closer to 008)
- eval_13: 0.5647 (similar logic)
- eval_01: 0.5664 (tie with all)

**The priors don't compose additively on eval_01.** Each design lifts
its specific evals at the *cost* of mildly losing others. The mosaic
averages those — no peak.

## Theory update
- **eval_01 hard ceiling around 0.568.** Five element-level designs
  (008, 011, 012, 013, 014) clustered between 0.562–0.569.
- Different evals reward different priors → no single mix can win
  everything. eval_01 specifically seems insensitive to the priors I
  can vary with element-level sampling.
- The lever I haven't tried: **data quality filtering**. Maybe a
  small fraction of low-complexity / poly-A / repeat-rich sequences
  is dragging the model's training. 0.5% of DHS have 3-mer diversity
  < 40, 1.6% have homopolymers ≥ 15.

## Numbers
mean_r: 0.554
eval_01: 0.5664
eval_07: 0.5881
eval_13: 0.5647
eval_08: 0.2008
time_s: 26

## Next
Test data-quality filtering: 008 design but exclude low-complexity
sequences (3-mer diversity < 45 or homopolymer ≥ 12) from both
DHS and cCRE pools before sampling.
