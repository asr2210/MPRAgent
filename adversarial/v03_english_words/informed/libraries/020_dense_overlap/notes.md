# 020 — High-density motif packing with OVERLAP allowed

## Plan
8 motifs/seq, each can overwrite earlier ones. Same 289-pool. Tests whether
saturating sequence with motif content (at the cost of motif fidelity)
helps. Qualitatively different from 008-family (3 strict-non-overlap motifs).

## Result
**eval_01 = 0.4284.** TIES exp 008 (0.4283). First library to match the
plateau best, suggesting dense packing is at least as good as careful
3-motif placement.

K562: 0.599 (↑ 0.003), HepG2: 0.629 (≈), SK-N-SH: 0.057 (~).

eval_07: 0.4336 (higher than 008's 0.4314).
eval_13: 0.4321 (also high).

Average across all 14 evals: 0.4255 (vs 008's 0.4253) — marginally higher.

## What this teaches
- The strict-no-overlap constraint of 008-family wasn't a critical
  feature. Overlap-tolerant dense packing is COMPETITIVE.
- Saturating the sequence with motif content (8 attempted inserts) gives
  similar K562/HepG2/SK-N-SH signal as 3 careful motifs.
- The surrogate appears to extract TF presence features from messy /
  overlapped motif content — supports T11 (presence matters more than
  fidelity).

## Theory T15
The surrogate is robust to MOTIF FIDELITY. As long as recognizable TF
patterns appear in the sequence, it can extract them. Packing more motifs
(even with overlap-corrupted instances) doesn't hurt because:
1. With overlap rate ~30%, ~5-6 motifs survive intact per sequence.
2. The corrupted "fragments" may also contribute to TF presence detection.

## Next
Exp 021: push density further. 15 motifs/seq, overlap OK. Tests whether
density scales further. If 021 > 0.435 → density is a continuing lever.
If 021 ≈ 0.42 → 8 was the sweet spot.
