# 024 — HYBRID library (25k @ 8 motifs + 25k @ 10 motifs)

## Plan
Library-level multiplex: 25k seqs from exp 020 design (K562/HepG2 best),
25k from exp 022 (SK-N-SH best). Same 289 pool, both with overlap.

## Result
**eval_01 = 0.4253.** Between parents (020=0.4284, 022=0.4247).
K562: 0.595, HepG2: 0.624, SK-N-SH: 0.057. Essentially averages 020+022.

## What this teaches
- Library-level multiplexing → surrogate trains on a mixed distribution
  and learns features from both, but the FINAL cell-type performance is
  averaged, not max'd.
- No emergent boost from combining designs. The surrogate doesn't sub-
  divide its learning by sequence type.

## Theory T17 (updated)
- Lever stacking fails: combined SK-N-SH levers anti-stack (exp 023).
- Library multiplexing also fails: hybrid averages instead of taking the
  per-cell-type max (exp 024).
- 008-family plateau at ~0.4253 (14-eval avg of 008) is robust to library-
  level composition tricks too.

## Next
Different angle: structured per-cell-type DENSE libraries.
- Exp 025: 15k seqs with 8 K562 motifs/seq, 15k with 8 HepG2 motifs/seq,
  15k with 8 SK-N-SH motifs/seq, 5k random. Tests T16 (per-cell-type
  density preference) with proper cell-type clustering at high density.
  If SK-N-SH benefits from 8 of its OWN motifs per seq, mean_r might
  improve.
