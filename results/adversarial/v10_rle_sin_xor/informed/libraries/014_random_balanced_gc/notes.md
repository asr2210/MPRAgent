# Experiment 014: random with GC=50% per sequence (exact 100 GC, 100 AT)

## Result
- eval_01: **0.5195** mean_r. K562=0.9943, HepG2=0.5641, SKNSH=0.0001

## Interpretation
Forcing GC=50% per sequence is neutral compared to fully random (0.521).
Per-sequence GC variance is not load-bearing for the model.
SKNSH crept up to ~0 (vs slightly negative for fully random) — possibly noise.

Note: this still has free A↔T and C↔G choices within positions, only constrains
the AT/GC counts. So the entropy reduction is modest.
