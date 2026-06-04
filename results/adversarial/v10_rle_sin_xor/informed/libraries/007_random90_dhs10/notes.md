# Experiment 007 — 90% random + 10% DHS mix

## Result
eval_01 = 0.4996. K562=0.93 (vs 0.99 random!), HepG2=0.567, SK-N-SH=0.

## Verdict
Even 10% DHS contamination caused a substantial K562 drop (0.99→0.93).
Tiny HepG2 gain (~0.003) doesn't compensate. **Mixed library underperforms pure random.**

## Insight
The K562 oracle is hypersensitive to non-random content. Any biological content
reduces K562 r disproportionately to its mass fraction. Pure random is K562-optimal.
