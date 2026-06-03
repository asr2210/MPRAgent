# 023 — 10 motifs + 60/40 weighted pool (combined SK-N-SH levers)

## Plan
Stack two SK-N-SH-friendly levers from T13 + T16:
- 10 motifs/seq with overlap (T16: SK-N-SH peaks at density 10)
- 60% cell-type pool + 40% full JASPAR (T13: SK-N-SH likes broader pool)

Hope: SK-N-SH at 0.07+ while K562/HepG2 hold near 0.59/0.62.

## Result
**eval_01 = 0.4249.** K562: 0.598 (↑), HepG2: 0.631 (↑), **SK-N-SH: 0.045**
(↓ from 022's 0.061).

Eval_04 = 0.4334 (high), Eval_09 = 0.4334 (high) — top single-eval values.

## What this teaches
- SK-N-SH levers don't STACK: density 10 alone (022) gives 0.061, weighted
  alone (017) gives 0.057, both together gives 0.045 (worse than either!).
- Combining the levers reduces effective density of CELL-TYPE motifs to
  ~6 (60% of 10), which is below the SK-N-SH-sweet-spot of 10. So the
  combination is actually worse for SK-N-SH.
- K562/HepG2 BENEFITED from the weighting (higher than 020/022). Broader
  pool adds variety they can use.

## Theory T17
SK-N-SH wants HIGH DENSITY OF NEURAL MOTIFS specifically — not just dense
sequences, and not just broader pools. Diluting cell-type density via
broader pool weighting cancels the density effect.

K562/HepG2 are more flexible — they like either narrow density-3 (008)
or weighted density-10 (023).

## Next
Test if stacking SK-N-SH levers in a different way works:
- Exp 024: 12 motifs/seq narrow pool (just 289), pure density push for
  SK-N-SH. Tests if 12 narrow > 10 narrow.
- Or hybrid: 50% exp 020 + 50% exp 022 designs.
