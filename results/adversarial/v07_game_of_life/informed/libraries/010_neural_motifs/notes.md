# Experiment 010 — Neural TF motif implants (targeted SK-N-SH lift)

## Design
50k random GC=0.50 200bp scaffolds. For each, implant 3-5 JASPAR PFM
instances drawn from 128 curated neural-lineage TF motifs (NEUROG, NEUROD,
ASCL, OLIG, SOX, POU3F, LHX, ISL, ZIC, PAX, OTX, FOX, NRSF, NHLH, DLX,
TFAP2, HAND, PHOX, GATA3, RFX, etc.).

GC: mean 0.485, std 0.039.

## Hypothesis
SK-N-SH (stuck at 0.13-0.15 across all 9 libraries) might lift with
targeted neural motif augmentation. If it does, motif targeting is a new
lever. If not, the SK-N-SH ceiling is structural.

## Result
- eval_01 = **0.3927** (vs 002 random JASPAR motif = 0.3939, gc_50 = 0.397)
- Mean across 14 evals = **0.3818** (vs 002 = 0.3827)
- Per-cell-type on eval_01: K562 0.612, HepG2 0.428, **SK-N-SH 0.138** (unchanged)
- Runtime: 1194s

## Interpretation
**Strong negative result for SK-N-SH liftability.** Neural-targeted motif
augmentation scored essentially identical to random JASPAR motif augmentation
(exp 002), with SK-N-SH unchanged at 0.138. The 0.13-0.15 SK-N-SH ceiling
is **structural in the evaluator, not training-content-dependent**.

This rules out:
- Cell-type-specific motif targeting as a SK-N-SH lever
- "SK-N-SH needs neural TFs" hypothesis
- The possibility of differential per-cell-type improvement via motif
  selection (at this scale)

## Theory v6 update — confirmed
The per-cell-type ceilings (K562 ~0.62, HepG2 ~0.43, SK-N-SH ~0.15) are
properties of the evaluator's held-out sequences, not of training data.
The mean_r ceiling at ~0.397 is correspondingly structural.

To break it, we need a qualitatively different mechanism — most likely
training-distribution shift towards eval-distribution via:
1. Hyperactive sequence generation (random scaffolds + targeted SA design)
2. Mixed-strategy library (combining multiple known-good approaches)
3. Larger candidate pool to find more extreme activity examples

Targeted motif augmentation is no longer worth more experiments.
