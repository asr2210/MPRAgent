# Experiment 011 — celltype_targeted_mix

## Design
25 k cCRE class-balanced + 25 k DHS with components weighted for
HepG2/SK-N-SH relevance (Neural, Cancer/epithelial, Digestive, etc.).
Hypothesis: lift HepG2 + SK-N-SH on eval_01 without losing K562 → net
eval_01 improvement.

## Result — 011 vs 008 by cell type (eval_01)

|         | 008    | 011    | Δ      |
|---------|-------:|-------:|-------:|
| K562    | 0.6227 | 0.6198 | −0.003 |
| HepG2   | 0.5350 | 0.5333 | −0.002 |
| SK-N-SH | 0.5436 | 0.5534 | +0.010 |
| **mean**| 0.5671 | 0.5688 | +0.002 |

eval_07 mean: 0.5833 → 0.5989 (+0.016, real)
eval_13 mean: 0.5603 → 0.5769 (+0.017, real)
eval_08:      0.2111 → 0.1813 (−0.030, smaller)
overall mean: 0.554  → 0.555  (statistical tie)

## Interpretation
Cell-type-aware DHS sampling moved **SK-N-SH** up by 0.010 (the
specific component oversampling worked, slightly). HepG2 was unmoved —
either Cancer/epithelial / Digestive component labels don't capture
HepG2-specific accessibility, or HepG2's chronic disadvantage is not
about library composition.

**Critical finding**: K562 stayed at 0.62 even though I deliberately
under-weighted Lymphoid + Myeloid components. The K562 advantage
**does not come from K562-relevant elements being over-represented in
the DHS pool**. K562 must be intrinsically easier to predict (cleaner
ENCODE measurements, higher signal-to-noise, or trained-model bias
within the harness). Library design cannot close this gap.

eval_07 and eval_13 lifted nicely (+0.016, +0.017). These two evals
seem to *reward* cell-type-aware sampling — interesting because they
correlate (instructions noted eval_07/eval_13 pair). Cell-type
diversity matters more for them than for eval_01.

## Theory update
- **K562 advantage is intrinsic, not compositional.** Stop trying to
  close the K562/HepG2/SK-N-SH gap. The path to higher eval_01 is via
  general library improvements, not cell-type rebalancing.
- **eval_07 + eval_13 are cell-type-diversity-rewarding.** Component-
  balanced sampling helps them more than uniform/random does.
- eval_01 is plateauing around 0.57. The gains from DHS→cCRE→mix→
  cell-type-aware are all in the +0.001 to +0.004 range. Need a
  qualitatively new direction to break through.

## Numbers
mean_r: 0.555
eval_01: 0.5688
eval_07: 0.5989
eval_13: 0.5769
eval_08: 0.1813
time_s: 29

## Next
Three ideas:
1. **DHS:cCRE ratio sweep** — try cCRE-heavy (15 k DHS + 35 k cCRE)
   to see if more cCRE-balanced content helps further.
2. **Positional diversity** — same elements, different windows
   (summit-centered + summit-offset).
3. **Mutation augmentation** — k-mer perturbation to densify local
   sequence space.

Starting with (1) — it's the simplest test of "is 50/50 actually
optimal" and directly extends the 008 finding.
