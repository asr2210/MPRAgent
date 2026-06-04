# Experiment 016 — Bimodal per-sequence GC, library mean GC=0.50

## Design
- 25k sequences sampled with per-base GC prob 0.30 → per-seq GC ≈ 0.30
- 25k sequences sampled with per-base GC prob 0.70 → per-seq GC ≈ 0.70
- Library MEAN GC = 0.500, STD = 0.202 (bimodal)

## Result
- eval_01 = **0.3760** (Δ vs 007 = -0.021)
- mean_r = **0.3674** (Δ vs 007 = -0.019)
- Per-cell: K562=0.584, HepG2=0.415, SK-N-SH=0.130

## Interpretation
**Per-sequence composition matters; library-mean is NOT sufficient.**
Even though the library averages to GC=0.50 (matching 007), each
individual sequence is at the composition extremes (0.30 or 0.70),
and the trained model performs as poorly as PLS-only (mean GC 0.62)
or motif-tiled libraries.

This says the constraint is per-sequence natural-likeness. The
evaluator's good performance on GC≈0.50 training data is because each
TRAINING sequence individually resembles eval sequences in composition,
not because the LIBRARY does.

## Theory v10 (refined)
- The model is sensitive to per-sequence composition, not just library
  composition distribution.
- The natural-like subspace requires per-sequence GC ∈ ~[0.45, 0.55]
  (broad but not arbitrary).
- Combined with v9: the natural-like subspace is BOTH per-sequence
  composition matched AND has natural-style local structure (no dense
  motif tiling).

## What this rules out
- Library-aggregate composition design (it must be per-sequence)
- "Extreme compositions cancel out" hypothesis
- Bimodal training as a diversity-providing strategy

## What to try next
1. **Multi-seed oracle pool** (017): combine top selections from 5 different
   seed runs of 007 — does noise averaging help?
2. **GC=0.40 / GC=0.60 uniform** (020/021): map composition lever precisely
   between the obvious extremes
3. **Natural cCRE GC distribution** (019): cCREs with no GC filter —
   natural composition variation per sequence
4. **Pan-active oracle** (018): select sequences with high MIN cross-cell
   activity — biases toward cell-type-transferable
