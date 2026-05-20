# 001 — Random baseline

## Goal
Establish the floor. Predicted essentially-zero model performance because random
200bp sequences should contain almost no functional regulatory elements.

## Method
- 150,000 uniformly random sequences, 200bp, ACGT, seed=0 (50% GC, every base i.i.d.).
- `prepare.py libraries/001_random_baseline/sequences.txt`.

## Result
Mean across 14 evals ≈ **0.82** (range 0.739 – 0.908). Runtime 57m14s (3434s).

| eval | mean_r | k562   | hepg2  | sknsh  |
|------|--------|--------|--------|--------|
| 01   | 0.7760 | 0.7650 | 0.7806 | 0.7824 |
| 02   | 0.8718 | 0.8620 | 0.8679 | 0.8854 |
| 03   | 0.8562 | 0.8474 | 0.8510 | 0.8701 |
| 04   | 0.8178 | 0.8099 | 0.8239 | 0.8196 |
| 05   | 0.7756 | 0.7637 | 0.7807 | 0.7823 |
| 06   | 0.8722 | 0.8622 | 0.8685 | 0.8858 |
| 07   | 0.7910 | 0.7821 | 0.7912 | 0.7996 |
| 08   | 0.9080 | 0.9056 | 0.9053 | 0.9131 |
| 09   | 0.8890 | 0.8804 | 0.8937 | 0.8928 |
| 10   | 0.8670 | 0.8652 | 0.8614 | 0.8745 |
| 11   | 0.7623 | 0.7521 | 0.7684 | 0.7664 |
| 12   | 0.7394 | 0.7334 | 0.7409 | 0.7438 |
| 13   | 0.7762 | 0.7610 | 0.7712 | 0.7964 |
| 14   | 0.8722 | 0.8625 | 0.8684 | 0.8858 |

## What this means
This is a **major** surprise. A model trained on pure random sequences with their
MPRA measurements predicts activity on held-out real sequences with r ≈ 0.82.
The MPRA must produce a meaningful spread of activity even on random DNA —
random 200bp regions contain accidental weak TF binding sites, varied GC,
CpG-like patterns, etc., and 150k examples is apparently enough to learn the
"soft" grammar.

The floor is much higher than expected. The bar to beat is r ≈ 0.82.

## Eval clustering (potential)
Several evals report near-identical numbers, suggesting they share data or
scoring:
- {02, 06, 14}: 0.8718 / 0.8722 / 0.8722 — essentially identical
- {01, 05}:    0.7760 / 0.7756
- {03, 10}:    0.8562 / 0.8670 (similar but not identical)

Worth tracking whether these always move together — if so, they're effectively
one signal each.

## Per-cell-line
SKNSH ≥ HepG2 ≈ K562 in most evals. The model isn't systematically worse on
any one cell line. Generalisation across cell types from random seems OK.

## Cost
~57 min wall-time per experiment on this GPU. 30 experiments = ~28h compute.
Need to be selective.

## Next
Experiment 002: real human genomic sequences (tiled 200bp windows from
hg38 chr22). Tests whether evolution-shaped sequence space teaches the model
more than uniform random. If yes → fold in more genome. If no → random
sequence space is surprisingly hard to beat for general grammar.
