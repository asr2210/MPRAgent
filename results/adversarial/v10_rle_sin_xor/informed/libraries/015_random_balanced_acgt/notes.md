# Experiment 015: random with exact 50/50/50/50 ACGT per sequence — BIG FINDING

## Result
- eval_01: **0.2697** (HUGE DROP from random's 0.521)
- K562=0.9928 (unchanged)
- **HepG2=-0.1824 (NEGATIVE — inverted prediction!)**
- SKNSH=-0.0013

## Interpretation — major finding
Removing per-sequence compositional variance INVERTS HepG2 prediction.
- K562 oracle unaffected (still saturated at 0.99)
- HepG2 model learned a compositional cue (likely GC% related) that, when
  removed from training, makes predictions anti-correlate with truth.

**The HepG2 model's signal is largely driven by per-sequence composition variance.**
When every sequence has identical (A,C,G,T) counts, the model has no input
signal to differentiate sequences for HepG2 prediction, and the residual
prediction direction is wrong.

## Implication for theory
Per-sequence COMPOSITION DIVERSITY is critical for HepG2.
- Removing it: HepG2 → -0.18
- Standard random uniform: HepG2 → 0.56
- Question: could AMPLIFYING composition variance (more extreme per-sequence
  GC differences) PUSH HepG2 above 0.56?

This is the strongest lead so far. Next experiment: amplify composition variance.
