# Experiment 007 — multi-source diverse library

## Design
- 10k cCRE uniform
- 4k cCRE CTCF-only boost
- 4k cCRE DNase-H3K4me3 boost
- 5k Table_S2 (UKBB+GTEx, non-eval-chrom)
- 4k DHS-topic
- 3k random synthetic
- 20k paired flanks (from cCRE+DHS positives)
- Total = 50,000

## Result — eval_01 = 0.1592, mean14 = 0.1532 (BOTH new personal bests)
- K562_r = +0.011 (eval_01), +0.089 (eval_06/11)
- SKNSH = 0.456 (eval_01), 0.507 (eval_06)

### Per-eval vs my prior best
- WIN: eval_01 +0.011, eval_02 +0.012, eval_03 +0.018, eval_06 +0.069,
  eval_11 +0.069, eval_12 +0.018, eval_14 +0.011
- LOSS: eval_04 -0.014, eval_07 -0.062, eval_09 -0.014, eval_10 -0.065,
  eval_13 -0.053

## Interpretation
Combining sources covers the UKBB/GTEx eval space (02/03/05/06/11/12) at
cost of chr-held-out evals (07/10/13). Net mean_r still rose because
UKBB/GTEx wins were larger.

K562 unlock on enhancer evals (06/11) reached +0.089 — biggest yet. This
is the CTCF/DNH3-boost effect, amplified by Table_S2's variant context.

## Theory v5
Multi-source library "covers eval distribution union." Specific composition
ratios matter: too much DHS dilutes ground-truth signal, too little
Table_S2 forgoes UKBB/GTEx wins.

## Next
Exp 008: triple Table_S2 share. Test if eval_07 recovers without losing
eval_06/11.
