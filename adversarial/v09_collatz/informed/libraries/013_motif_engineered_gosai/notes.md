# E013 — Motif-engineered Gosai backbones

Plant 2 random cell-type-specific consensus motifs (GATA1/KLF1/HNF4A/HNF1A/
FOXA1/ASCL1/NEUROD1/REST) per Gosai backbone, random position + orientation.

## Result
eval_01 = 0.3189 (vs E004 random Gosai = 0.3233). K562=0.148, HepG2=0.200,
SKNSH=0.609. Per-cell pattern unchanged.

## Interpretation
Injecting strong consensus motifs slightly HURTS (-0.004). The pipeline's
K562/HepG2 ceilings are not raised by extra motif content; the bottleneck
is downstream (model capacity, label noise on those cells), not motif
representation in input. Confirms structural ceiling.

Slight regression likely because injection scrambles the natural Gosai
co-occurrence structure that the model would have learned cleanly.
