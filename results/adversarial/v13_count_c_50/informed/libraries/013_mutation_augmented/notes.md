# Experiment 013 — mutation_augmented

## Design
25 k unique elements (12.5 k DHS uniform + 12.5 k cCRE class-bal)
+ 25 k mutated copies (5 random SNPs each). Total 50 k.

## Hypothesis
Local sequence-space densification (paired original + 5-SNP mutant)
might teach the model robustness to small variants and tighten its
sequence-feature mapping. Could lift eval_01 beyond the 0.568 plateau.

## Result
| eval | 008 (50 k unique) | 010 (RC aug)  | **013 (mut aug)** |
|------|------------------:|-------------:|------------------:|
| 01   | 0.5671 | 0.5659 | 0.5624 (−0.005) |
| 07   | 0.5833 | 0.5852 | 0.5777 (−0.006) |
| 13   | 0.5603 | 0.5621 | 0.5538 (−0.007) |
| 04/09| 0.5723 | 0.5674 | 0.5676 |
| 08   | 0.2111 | 0.2083 | 0.2161 |
| mean | 0.554  | 0.553  | 0.551 |

## Interpretation
**Mutation augmentation slightly hurts**, very similar to RC
augmentation (010) but a touch worse. Confirmation: **unique element
diversity > local sequence-space density**. The model already learns
local sequence invariances implicitly; adding paired near-duplicates
trades diversity for a redundant signal.

Combined with 010's result:
- RC aug: −0.001 (wash)
- Mut aug: −0.005 (slight loss)
- Both: trade real unique elements for redundant pairs

**The 50 k budget is best spent on 50 k unique elements.**

## Theory update
- **Augmentation strategies are net-negative on this task.** Both RC
  and SNP mutation hurt or tie. The bottleneck is *information per
  element*, not local variability.
- Anything that consumes the 50 k budget on near-duplicates is going
  to underperform. Skip augmentation entirely from now on.

## Numbers
mean_r: 0.551
eval_01: 0.5624
eval_07: 0.5777
eval_08: 0.2161
time_s: 22

## Next
Augmentation is a dead end. Element-level design is plateaued at 0.568.
Try a **mosaic** combining priors that each lift different evals:
- 011's component-targeted DHS lifts eval_07/eval_13 (+0.016/+0.017)
- 012's cCRE-heavy lifts eval_04/eval_09 (+0.005)
- 008's 50/50 mix is the eval_01 baseline
Test whether composing these priors (e.g., cCRE-heavy + component-
targeted DHS) gives additive lift on eval_01 too.
