# Experiment 015 — max-SE quality filter (stricter)

## Method
max(K562_se, HepG2_se, SKNSH_se) < 0.3. All cells must have clean
measurements. Random sample from 465K passing sequences.

## Result
**eval_01 = 0.0094** — WORSE than mean-SE filter (014 = 0.0168).

## Interpretation
Stricter per-cell quality HURTS. Probably because max-SE filter biases against
cell-specific regulatory sequences (which by definition have one cell with weak
or noisy signal). The library becomes constitutively-active-biased — bad for
learning cell-specific regulation.

## Theory consolidated (post exp 014/015)
The "good library" needs:
- Mean-SE quality filter (~0.5) to remove garbage labels
- Natural distribution of cell-specificity (don't filter out cell-specific seqs)
- Natural distribution of activity (stratification optional, not essential)

Bottleneck: per-sequence information content. To improve we need MORE
INFORMATION PER SLOT than measured Gosai variants provide. Next direction:
data augmentation (revcomp, mutation) that recycles measured sequences with
the oracle's predictions.
