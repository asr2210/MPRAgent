# Experiment 013: 25k random + 25k revcomp

## Result
- eval_01: **0.5190** mean_r. K562=0.9947, HepG2=0.5650, SKNSH=-0.0027

## Comparison
- Pure random seed 1: 0.5210
- 25k+revcomp seed 1: 0.5190
- Slight drop (-0.002)

## Interpretation
Revcomp augmentation is approximately neutral. The marginal information from
25k unique random sequences with both strands ≈ 50k unique random sequences.
No strand-symmetry inductive bias benefit.
