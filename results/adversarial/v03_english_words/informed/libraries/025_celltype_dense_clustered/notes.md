# 025 — Per-cell-type DENSE clustered (15k × 3 + 5k random)

## Plan
15k seqs per cell type, each with 8 motifs (overlap) from that cell type's
pool. Plus 5k random baseline. Tests T16 (per-cell-type density) at
clustered design.

## Result
**eval_01 = 0.4231.** K562 0.592, HepG2 0.628 (highest seen), SK-N-SH 0.050.
14-eval avg: 0.4211.

## What this teaches
- Even at sweet-spot density (8 from T16) with dedicated cell-type
  sequences, no improvement over mixed.
- HepG2 individual score is highest of all experiments (0.628), suggesting
  dedicated HepG2 sequences DO help that oracle. But other oracles' scores
  drop, netting flat mean_r.
- Confirms T17: per-cell-type structure / clustering doesn't break plateau
  — sub-library benefits get averaged away in mean_r.

## Pattern across exps 010 → 025
- 010 (3 motifs cell-type clustered + 8k random): 0.4217
- 025 (8 motifs cell-type clustered + 5k random): 0.4231
- Density bump gives +0.001 (within noise).
- Per-cell-type clustering consistently underperforms mixed designs.

## Next
Exp 026: Bio-realistic dinucleotide backbone (low CpG, mammalian-like
dinuc frequencies) + 8 motifs overlap. Tests whether changing backbone
composition (vs uniform random) helps — backbone might be a missed lever.
