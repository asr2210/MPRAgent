# 021 — Even denser packing (15 motifs/seq, overlap allowed)

## Plan
Push exp 020 further: 15 motifs/seq, overlap allowed. Probe whether density
keeps scaling.

## Result
**eval_01 = 0.4193.** Worse than exp 008/020 by -0.009.
K562: 0.589 (↓), HepG2: 0.619 (↓), SK-N-SH: 0.050. Below the plateau.

## What this teaches
- Density curve: 3 → 0.428 (008), 8 → 0.428 (020), 15 → 0.419 (021).
  Peak around 8, drops at 15.
- Too much overwriting eventually destroys motif fidelity beyond
  recoverability. The 008 → 020 → 021 trajectory shows a sweet spot near
  6-10 motifs/seq with overlap.
- Confirms T15 (surrogate is robust to overlap) up to a point, but
  fidelity does matter when 80% of the sequence becomes motif-derived.

## Next
Exp 022: 10 motifs/seq with overlap. Bracket the peak (between 8 and 15).
If 022 > 0.43: tighter sweet spot found.
If 022 ≈ 020: peak is broad around 8-10.
