# Experiment 007 — ccre_uniform

## Design
Uniform random 50 000 of the 2.35 M cCRE pool (no class balancing).

## Hypothesis
Diagnostic for exp 006. Tells me whether the +0.031 mean improvement
in 006 came from (a) using cCRE rather than DHS, or (b) forcing class
balance.

## Result
| variant | eval_01 | mean(14) |
|---------|--------:|---------:|
| 002 uniform DHS        | 0.5627 | 0.516 |
| 007 uniform cCRE       | **0.5532** | 0.534 |
| 006 cCRE class-balanced| 0.5637 | 0.547 |

## Interpretation
The improvement is **mostly class balance, not source**.
- Uniform cCRE actually scores **lower than uniform DHS** on eval_01
  (−0.009). Source-only effect: cCRE is slightly *worse* than DHS at
  the standard random-draw level.
- Mean across all 14 evals: cCRE uniform 0.534 vs DHS uniform 0.516
  → cCRE source is +0.018 better on average (helps non-eval_01 evals).
- Class balance on cCRE adds another +0.011 on eval_01 and +0.013 on
  mean.

The class-balance gain is consistent with the theory: forcing equal
representation of PLS / pELS / dELS / CTCF / TF /etc. injects functional
diversity that uniform sampling misses (because dELS would dominate at
62 %).

## Why uniform cCRE underperforms uniform DHS on eval_01
Uniform cCRE = 62 % distal enhancers + 10 % proximal enhancers + a few
% each of other classes. That's a much narrower functional diversity
than a uniform draw from the open-chromatin universe (DHS), even though
each individual cCRE element is more reliably regulatory. Quality but
not enough range.

## Theory update
- Combining "more functional certainty per element" (cCRE) with
  "broader functional range" (class balance) wins.
- Pure DHS is good *because* it's wider, not because it's cleaner.
- For mixing sources: DHS likely helps on evals where uniform DHS
  excels (eval_07/13) and cCRE on evals where uniform cCRE excels
  (eval_04/08/09). A blended library should beat either alone.

## Numbers
mean_r averaged across 14 evals: 0.534
eval_01: 0.5532
eval_08: 0.2426
time_s: 32.3

## Next
008 = 50/50 DHS uniform + cCRE class-balanced mix. Tests complementarity.
