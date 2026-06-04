# 018 — Per-sequence cell-type-balanced motifs

## Plan
Every sequence has EXACTLY 1 K562 + 1 HepG2 + 1 SK-N-SH motif. 3 motifs/seq.
Hypothesis: guaranteed per-cell-type coverage in every training sample
should help each oracle. Tests T13.

## Result
**eval_01 = 0.4223.** Essentially equal to 016 (0.4222 = noise reference).
K562: 0.594, HepG2: 0.626, SK-N-SH: 0.046.

## What this teaches
- Per-sequence cell-type balance does NOT break the plateau.
- The 008-family's uniform sampling already gives each cell type sufficient
  exposure across the 50k library; concentration per sequence isn't the lever.
- T13's "trade-off between cell types" persists even with forced balance.
  K562/HepG2 stayed similar to plateau, SK-N-SH stayed at floor.
- The structural plateau is real; cell-type balance doesn't unlock it.

## Theory T13 refined
The plateau at ~0.422-0.428 is robust to MANY librarian choices:
density, pool size, per-seq composition, balance, weighting. To break it
likely requires a structurally different SEQUENCE design (not motif-pool
or counting change). Possible levers:
- Motif length distribution
- Motif spacing / syntax (enhancer-like clusters)
- Backbone composition (not uniform random)
- Sub-library mixtures with very different members
- Iterative refinement (filter or score)

## Next
Exp 019: motif PAIR CLUSTERS at biological spacing (5-20bp between paired
motifs). Mimics real enhancer architecture where TFs cluster within ~50bp
windows. 4 motifs/seq arranged as 2 pairs.

This tests whether motif CO-OCCURRENCE within a small window is a feature
the surrogate uses (it should be, if the surrogate's mechanism is
biology-aware). If yes → structural change opens a new ceiling. If no →
spatial structure doesn't matter at this length scale.
