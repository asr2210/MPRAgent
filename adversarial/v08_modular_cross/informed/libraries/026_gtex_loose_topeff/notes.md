# Experiment 026 — GTEX lfcSE<0.7, top 50K by |mean activity|

## Method
Deterministic: rank GTEX lfcSE<0.7 sequences by |mean activity|, take top 50K
(cutoff |mean|>1.626). Removes seed variance.

## Result
**eval_01 = 0.0131** — well below GTEX-loose baseline.
**eval_04/09 = 0.0292** — NEW HIGH for these evals.

## Interpretation
Top-effect selection (deterministic) is a stronger version of activity extremes.
- Massively boosts eval_04/09 (now 0.029, up from random's 0.022)
- Hurts eval_01 (drops to 0.013 from random's 0.019)

Consistent with earlier finding: eval_01 prefers broad activity, eval_04/09
prefers high-magnitude effects.

## What this tells me
The two main eval clusters genuinely want different sequences. Cannot satisfy
both with one 50K library — selection for one degrades the other.

## Final direction
With eval_01 as primary metric, stick with GTEX-loose random and accept the
0.018-0.022 band. Try a few last experiments to:
1. Push eval_01 specifically (cell-type-targeted selection?)
2. Validate the best configuration with replicate seeds
