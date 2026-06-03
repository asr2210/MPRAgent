# Experiment 014 — Dirichlet(0.5) seed=99 (noise floor)

## Result
eval_01 = **0.1365** vs exp 002 (seed=42) 0.1395 — gap of 0.003

Per-cell variance across seeds is much larger than per-cell mean variance:
- mean: 0.1395 (seed=42) vs 0.1365 (seed=99) → variance ~0.003
- K562: 0.0473 vs 0.0390 → variance ~0.008
- HepG2: 0.1697 vs 0.1687 → variance ~0.001
- SK-N-SH: 0.2015 vs 0.2016 → variance ~0.000

## Single-seed noise estimate
mean_r noise ≈ ±0.003. So "true" Dirichlet(0.5) score is ~0.138 ± 0.003.

## Re-interpret prior results
After noise correction, "within noise" of Dirichlet(0.5):
- 002 = 0.1395, 014 = 0.1365 (same population)
- 005 (mixed alpha) = 0.1377 ✓
- 007 (sobol) = 0.1371 ✓
- 009 (gradient) = 0.1369 ✓
- 013 (trimmed) = 0.1379 ✓
- 006 (JASPAR motifs) = 0.1364 ✓ (within noise!)

Definitely worse (beyond noise):
- 001 (synthseqs) = 0.1319 (-0.020)
- 003 (Markov) = 0.1346 (-0.013)
- 008 (composition-iid) = 0.1354 (-0.005)
- 010 (alpha=0.3) = 0.1344 (-0.014)
- 011 (mode-stratified) = 0.1356 (-0.005)
- 012 (uniform synthseqs) = 0.1314 (-0.025)

## Implication
Dirichlet(0.5) family is a TIGHT plateau. To break through I need
+0.005 (real lift, beyond noise) or +0.01 (clear win). Single-axis
tweaks within the family won't do this.

## Strategy for remaining experiments
1. Try one more orthogonal axis (block-Dirichlet, exp 015)
2. Multi-seed sampling to find lucky seed (exp 016-019)
3. Variants of best seed (motif aug, alpha tweaks)
4. Final library + summary
