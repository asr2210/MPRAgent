# 019 — Motif pair clusters (5-20bp spacing)

## Plan
Each sequence has 2 PAIRS of motifs. Within-pair spacing 5-20bp. 4 motifs/seq.
Tests whether surrogate uses motif syntax / co-occurrence at biological
spacing scales.

## Result
**eval_01 = 0.4198.** Within noise of plateau.
K562: 0.591, HepG2: 0.624, SK-N-SH: 0.044.

## What this teaches
- Spatial syntax doesn't help. Surrogate appears to use TF PRESENCE as the
  feature, not pair spacing or co-occurrence patterns.
- Combined with consistent results across exps 010-019, the 008-family is
  a STRUCTURAL OPTIMUM for this pipeline at ~0.42-0.43.
- Confirms T13: surrogate's learning is dominated by what TFs are present;
  spacing, density, balance, syntax are all near-noise interventions.

## Next
Strategically, since structural variations don't move the needle, try
RADICAL design changes:
- Exp 020: aggressive DENSITY with overlap allowed (7 motifs/seq, can
  overlap). Tests if packing more TF signal in per sequence helps.
- Exp 021: motif length stratification (long-only vs short-only).
- Exp 022: HYBRID library (50% pure random + 50% motif-enriched). Tests
  library-level diversification.

If exp 020-022 all plateau too, ceiling is fundamental and document T14
as "no library-level lever exists to clear 0.43 in this pipeline."
