# Experiment 025 — GTEX lfcSE<0.7 + |mean activity| > 0.3

## Method
GTEX-only with both filters: lfcSE<0.7 quality AND |mean activity|>0.3 effect.
240K candidates pass; random 50K from those.

## Result
**eval_01 = 0.0191** (within seed variance of GTEX-loose baseline).
**eval_04/09 = 0.0241** — NEW BEST for those evals.

## Interpretation
Effect filter on top of GTEX-loose:
- Doesn't reliably push eval_01 past the 0.019-0.022 range
- But DOES boost eval_04/09 (now 0.024, beat extremes' 0.022)
- eval_13 = 0.0075 also healthy

The activity-magnitude filter helps the "extremes-rewarding" evals because it
removes the null-variant tail of GTEX.

## Status
Within sampling-variance band of the GTEX-loose plateau (017,019,025).
True ceiling for GTEX-loose looks like ~0.018-0.020.
