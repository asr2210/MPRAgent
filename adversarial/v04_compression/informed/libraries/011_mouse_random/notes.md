# 011_mouse_random

## Design
50k random 200bp tiles from mouse mm10 chr1/11/17/19 (gene-rich,
~474Mb total), length-weighted. N-rejection. seed=0.

Diagnostic: is the real-DNA prior vertebrate-universal or human-specific?

## Result — strongly positive
- eval_01 = **0.4485** (vs E2 human 0.4992; Δ -0.051)
- mean over 14 ≈ 0.466
- eval_07 = 0.5438 (vs E2 0.5985; Δ -0.055)
- eval_13 = 0.5401 (vs E2 0.6020; Δ -0.062)
- eval_08 = 0.0866 (vs E2 0.0916; Δ -0.005)

Mouse DNA, with ZERO human-specific training, gets within 0.05 r-points
of human DNA. Δ is small and uniform across evals.

## Where this lands
```
synth uniform      0.307
synth Markov k4    0.268
motif implant      0.347
promoter TSS       0.370
DHS stratified     0.438
mouse random       0.449  ← cross-species real DNA
DHS uniform        0.470
all-chrom human    0.482
gene-rich human    0.499  ← my ceiling
```

Mouse is closer to DHS-uniform than to human random. Most of the
real-DNA prior IS vertebrate-universal. The -0.05 gap is the human-
specific component (likely diverged TFBS/repeat content).

## Theory update v5
The "+0.19 from real DNA" prior decomposes:
- Vertebrate-universal component (works in mouse): ~0.14
- Human-specific refinement (TFBS, ALU repeats, gene grammar): ~0.05

The human-specific component is small but real. This means:
- Mixing human + mouse could broaden effective coverage WITHOUT
  losing much human-specificity (if dilution is gentle).
- OR mixing could simply dilute the human signal and underperform pure
  human (if the model can't use both).

## Implications
1. Cross-species mixing is a viable direction for breaking 0.50.
2. The eval set is testing primarily vertebrate-universal grammar
   with a small human-specific component.
3. If mouse works this well, other mammals (zebrafish, dog) likely
   work too — the prior is broadly vertebrate.

## Plan for E12
**human + mouse 50/50 mix**: 25k from E2's human chr1/17/19/22 +
25k from this mouse design. Tests:
- Mix ≥ 0.50: cross-species diversity ADDS information beyond pure
  human. Continue this axis.
- Mix ~ 0.48: dilution = diversity gain. Mouse helps modestly or not
  at all when blended.
- Mix < 0.48: dilution dominates. Drop cross-species mixing.

Predicted: ~0.49 (between pure human 0.499 and pure mouse 0.449,
weighted to human's higher per-slot signal). Hopeful: ~0.51 if
diversity adds modeling regularization.

## What I'm NOT sure about
The "vertebrate-universal" hypothesis predicts mouse should generalize
well in EITHER direction. But this experiment trains on mouse and
evaluates on (presumably) human cell lines. The reverse — training on
human, evaluating on mouse — isn't testable in this pipeline. So the
asymmetry of the result could reflect the eval being human-biased.

Even so: a 0.45 score with cross-species training is a strong endorsement.
