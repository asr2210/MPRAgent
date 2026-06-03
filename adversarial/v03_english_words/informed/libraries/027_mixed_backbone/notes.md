# 027 — Mixed backbone library (25k random + 25k dinuc)

## Plan
Hedge backbone composition at library level: 25k random-uniform
backbone + 25k mammalian dinuc backbone. Both with 8 motifs (overlap)
from 289 pool. Tests T18 via averaging.

## Result
**eval_01 = 0.4181.** K562 0.585, HepG2 0.614, **SK-N-SH 0.056**
(lost the bump). 14-eval avg: 0.4181.

## What this teaches
- Library-level backbone hedging FAILS — averages to ~midpoint of the
  two parent designs (020=0.4284, 026=0.4153 → 0.4181 is below midpoint).
- SK-N-SH bump from dinuc backbone (0.062) gets washed out (0.056), not
  preserved. The K562/HepG2 oracles' negative response to dinuc seqs
  dominates per-eval mean.
- This is the SAME failure mode as exp 024 (hybrid dense). Both
  experiments confirm: library multiplexing averages performance,
  doesn't max it per cell type.
- **Strong evidence that mean_r over all 3 cell types is a saturated
  objective**: gains in one cell type must come without losses in
  others, but every "lever" we find is anti-correlated across types.

## Theory T19
The 008-family plateau (~0.422-0.428) represents the optimum for an
oracle that must score well on a MEAN of 3 anti-correlated cell types.
Each individual cell type oracle can be pushed higher (HepG2 0.628 in
025, SK-N-SH 0.062 in 026, K562 likely ~0.60 in 020), but no single
library can do all 3 simultaneously because their preferences conflict:
- K562: high GC, motif-rich, uniform background
- HepG2: high GC, motif-rich, uniform background
- SK-N-SH: composition-rich, low-CpG mammalian background

The plateau is a CONSTRAINT from the metric definition, not a library-
design failure.

## Next
If the plateau is metric-imposed, two remaining angles:
1. Find sequences that score well on ALL THREE simultaneously rather
   than relying on per-cell-type optima. Look for "universal enhancer"
   features (SP1+CTCF+TBP+NFY backbones with conserved motif spacing).
2. Diversity-driven library: oversample 100k seqs across many designs,
   keep most diverse 50k. The eval might reward diversity per se
   (richer training distribution → better generalization to held-out
   cells).

Exp 028: dense motif library using only HIGH info-content PFMs
(filter to top 30% by entropy). Tests whether the 289 pool dilutes
signal with noisy/weak PFMs.
