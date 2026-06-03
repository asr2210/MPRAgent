# 015 — Variable motif density (1-6 motifs/seq)

## Plan
Same 289-pool as exp 008. Per-seq motif count drawn from uniform[1, 6].
Mean density ~3.5, range covers exp 008 (3), 009 (5), and beyond.

## Result
**eval_01 = 0.4220.** Worse than exp 008 (0.4283) by -0.006.
K562: 0.592, HepG2: 0.625, SK-N-SH: 0.049.

## What this teaches
- Density mixing doesn't help.
- Confirms exp 008/009 finding: fixed-3 motifs/seq is near-optimal.
- A distribution of densities doesn't improve coverage; the surrogate
  benefits from CONSISTENT moderate density.
- Possibly: sparse sequences (1-2 motifs) are too close to noise; dense
  ones (5-6) overlap and waste pixels. Fixed-3 is a sweet spot.

## Theory T11 unchanged
TF identity diversity per sequence is the lever, not density variance.

## Next
6 consecutive variations of exp 008 (010-015) have all lost. Need to
estimate library NOISE: re-run exp 008 with SEED=42 (exp 016) to see how
much mean_r moves when only library randomization differs. If it varies
by ±0.005, my "regressions" may be real but small; if ±0.01, all my
recent comparisons are noise-dominated.
