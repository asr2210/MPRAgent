# E019 — Gosai dinucleotide-shuffled (Euler-walk per sequence)

Each Gosai sequence shuffled with dinucleotide-frequency preservation
(motifs destroyed but 2-mer stats kept).

## Result
eval_01 = 0.3204. Per-cell K562=0.146, HepG2=0.201, SKNSH=0.614.

## Interpretation — major finding
Dinucleotide-shuffled (0.3204) vs intact Gosai (0.3233): difference is
only **+0.003**. The full information from "real sequences" beyond
dinucleotide structure is essentially zero (within noise).

## Refined decomposition of pipeline ceiling
Starting from random uniform = 0.244:
- + GC composition matched to Gosai: +0.061 → 0.306 (E017/E016)
- + dinucleotide structure          : +0.014 → 0.320 (E019)
- + higher-order motifs/grammar     : +0.003 → 0.323 (E004 intact Gosai)
- + Gosai chr 7/9/13/21/X quirk     : +0.013 → 0.336 (E008)

## Striking implication
The pipeline's model is effectively a **dinucleotide-frequency regressor**.
Triplet motifs, transcription factor binding sites, and longer-range
grammar contribute negligibly to its predictions. This explains why:
1. K562/HepG2 ceilings are immune to library design
2. Motif injection (E013) didn't help
3. Source mixing (E009) didn't help
4. Quality filtering (E005, E011) had marginal effect

The 0.34 ceiling is genuinely bounded by the model's representational
capacity, not by library content.
