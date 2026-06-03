# Experiment 004 — DHS Neural-only (diagnostic)

## Hypothesis
Library of only neural-component DHS should move SK-N-SH r off zero.

## Result
eval_01 = 0.4850. K562=0.90, HepG2=0.56, SK-N-SH=0.0005 (still ≈ 0).
Slightly WORSE than synth_random on K562.

## Verdict
**SK-N-SH cannot be moved by changing the source distribution to neural-specific.**
Confirms SK-N-SH ≈ 0 is a pipeline-level constraint. Library design CANNOT impact
SK-N-SH r directly.

## Implications for theory
- Mean_r upper bound ≈ (0.99 + ~0.6 + 0) / 3 = 0.53 with random sequences
- To beat that I need to push HepG2 above 0.56 without losing K562
- Or find specific sequence properties the K562/HepG2 oracles strongly respond to
