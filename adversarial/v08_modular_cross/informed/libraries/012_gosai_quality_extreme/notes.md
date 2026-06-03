# Experiment 012 — Extreme quality filter (lfcSE<0.15)

## Method
Gosai with lfcSE<0.15 (vs 0.3 in 005, 0.5 in 004). 235K sequences pass.
Quintile-stratified by mean activity.

## Result
**eval_01 = 0.0108** — substantially below 0.018 plateau.

## Interpretation
Quality has a peak around lfcSE<0.5; going tighter HURTS. Why?
- Very tight SE biases toward sequences with: (a) high replicate count
  (oversampled), (b) low intrinsic activity variability, (c) "easy to measure"
- This narrows the distribution and removes hard-but-real examples that the
  model needs to learn generalizable features

So quality is not the limiting axis. The plateau is set by something else.

## Theory
- Composition (sub-source mix) matters: UKBB lifts eval_01, GTEX lifts eval_04/09
- Quality has a U-curve: too noisy fails, too clean loses diversity
- Selection on extremes/variance specializes for some evals at cost of others

Next probe: signal-to-noise ratio (|effect|/SE). Selects sequences with both
high-magnitude AND reliable effects. Different from low-SE alone or
high-|effect| alone.
