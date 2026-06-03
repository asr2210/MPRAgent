# Experiment 007 — activity extremes

## Method
25K most-active + 25K most-inactive from Gosai with lfcSE<0.3.

## Result
**eval_01 = 0.0142** (vs 005 plateau 0.0180). Lower mean but mixed pattern:
- eval_04/09 jumped to 0.0224 (BEST so far for these)
- eval_07/08/13 dropped to ~0 or negative

## Interpretation
Different eval sets prefer different sequence distributions. Extremes train
the model to predict extreme activity but lose middle-range information.
eval_04/09 may measure binary active/inactive classification; others measure
graded activity.

For maximizing eval_01 (primary), broad stratification (exp 004/005) still wins.
