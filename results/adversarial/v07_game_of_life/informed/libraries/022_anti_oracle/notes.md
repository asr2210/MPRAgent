# Experiment 022 — Anti-oracle (bottom 50k by max Malinois)

## Design
500k GC-50 random, score with Malinois, take BOTTOM 50k by max(K562, HepG2, SKNSH).
These are sequences where ALL THREE cells predict LOW or NEGATIVE activity.

## Selected stats
- Max-cell score: mean=-0.135 (vs 007's +2.6, random's +0.7)
- K562 mean=-0.166, HepG2=-0.264, SKNSH=-0.317 — ALL NEGATIVE
- GC: mean=0.489, std=0.035 (slightly DECREASED — low activity correlates with low GC)

## Result
- eval_01 = **0.3933** (Δ vs 007 = -0.0036, near 014 noise)
- mean_r = **0.3830** (Δ vs 007 = -0.0031, near random GC-50)
- Per-cell: K562=0.618, HepG2=0.432, SK-N-SH=0.130

## STUNNING finding
**Anti-oracle ≈ random ≈ 014.** Training on sequences with NEGATIVE
predicted activity gives essentially the SAME mean_r as random GC-50 or
the lucky 007 oracle. The "+0.005 oracle benefit" we thought we had
ALMOST DISAPPEARS when we use anti-direction.

This means:
1. The ORACLE DIRECTION is essentially irrelevant
2. What matters is the SUBSPACE (GC-50 random, independent samples)
3. The 007/017/018 mean_r of 0.386 is at most marginally above random
   (and may itself be near the upper-noise envelope of random)

## Theory v15 (major consolidation)
The 0.397/0.386 ceiling is achieved by ANY library satisfying:
1. ~25k+ independent samples (021 fails this)
2. Per-sequence composition in natural-like range (016 fails this)
3. Natural-like local structure (015 fails this)
4. From a sequence subspace that resembles eval distribution
   (003 PLS fails — wrong composition; 009 cCRE-oracle fails — biased)

GIVEN these constraints, the SELECTION CRITERION is irrelevant:
- Random | top | bottom | pan-active | discrim — all hit the same ceiling
- 022 (bottom) and 007 (top) differ by only 0.003 mean_r — noise

This dramatically reshapes the theory: the 0.397 ceiling is fundamentally
about the (model, evaluator, sequence subspace) tuple, not about training
data design within the subspace.

## What this rules out
- Oracle activity-based selection as a meaningful lever
- Almost any single-axis training-data variation

## What this opens
The remaining levers MUST operate on:
- Subspace fundamentals (composition, structure, independence)
- Cross-modal information (e.g., model ensemble agreement)
- Trainer-level features (which we cannot modify)

We have essentially mapped the plateau. Remaining experiments will
either find a NEW lever or definitively confirm the structural ceiling.
