# Experiment 009 — Gosai filtered for extreme |SKNSH_log2FC|

## Result
eval_01 = 0.4972 (worse than uniform Gosai 0.5031).
SK-N-SH = 0.013 (worse than uniform Gosai's 0.026).

## Verdict
**Filtering for extremes HURTS generalization.** Extreme SK-N-SH sequences may
be a biased subset (specific motifs/regions) that don't generalize to typical
test sequences.
