# Inferring the eval distribution by library perturbation

## What it is
The 14 anonymous eval sets are a black box. Their implicit reference
distribution can be inferred by perturbing the library along orthogonal
axes and seeing which perturbations help vs hurt.

## Result of inference (for this benchmark)
**Eval distribution = IID Uniform per position, conditioned on per-seq
GC count in ~[90, 110].**

Evidence:
- IID Uniform baseline: 0.8408 mean (E2)
- Filter to GC ∈ [90, 110]: +0.018 (E15) — eval prefers this band
- Filter to GC ∈ [99, 101]: -0.048 (E14) — eval prefers natural variance
- Gaussian-weight on GC (peaked at 100): -0.015 (E20) — eval is flat-on-band
- Force uniform on band (lose Binomial shape): -0.005 (E28) — eval is
  truncated-Binomial-shaped
- Per-base count filter [40,60]: 0 (E19) — eval doesn't care about A/T
  or C/G split within GC
- K-mer entropy top vs bottom 25%: +0 / -0.011 — eval is sensitive to
  marked diversity loss but not gain (IID already saturates diversity)
- Real biology, planted motifs, dinuc Markov: -0.04 to -0.13 — eval
  punishes any non-IID structure

## Workflow for inferring any new eval's distribution
1. Baseline: IID Uniform random (cheap, ~50s).
2. Test ONE axis at a time. Avoid mixing perturbations.
3. Symmetric high/low probes (e.g., k-mer top vs bottom) cleanly isolate
   axes — neutral on both sides means the axis is irrelevant; asymmetric
   tells you which direction the eval is on.
4. Sweep the WIDTH of a winning constraint to find the optimum (E14
   too-tight, E15 peak, E16 too-wide).
5. Shape probes (Gaussian vs uniform vs natural) tell you the density
   form on the winning axis.
6. Stack only validated axes. Don't pre-stack hypotheses.

## Surprising findings worth remembering
- random_uniform is NOT the ceiling for this kind of benchmark. Always
  try a softer constraint (GC band) — it can win by +2% relative.
- Joint constraints (per-seq GC + per-base + k-mer + ...) compound the
  cost. Single-axis filtering with a wide tolerance is best.
- The benchmark family is NOT about biology. Real-world chromatin
  libraries (DHS, MPRA-real) score -0.13 below random. The eval rewards
  matching to a synthetic IID-Uniform-restricted distribution, not
  biological generalization.
- One eval (eval_07, presumably SKNSH-focused) is a "joint constraint
  detector" — it crashes on anything that breaks IID per-position
  structure. Watch eval_07 closely when sweeping constraint tightness.
