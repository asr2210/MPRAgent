# Experiment 022 — push CpG slice to 10k (20%)

## Result
eval_01 = 0.5710 (−0.0036 vs 021's 0.5746). Too much CpG hurts.

eval_08 keeps rising: 022 = 0.2167 (+0.033 vs 021). Linear in CpG slice.
eval_04/09: 0.5776 (best ever for 04/09, same as 012).
eval_13: 0.5669 (loss).

## Interpretation
The CpG slice has roughly linear effect on eval_08 but a concave-down
effect on eval_01: 0 % (015)=0.5736, 12 % (021)=0.5746, 20 % (022)=0.5710.

**021's 12 % CpG is the sweet spot.** More CpG steals from the orth-
DHS diversity that lifts non-K562 cells.

## Numbers
mean_r: 0.556
eval_01: 0.5710
eval_04: 0.5776 (best for 04/09)
eval_08: 0.2167
time_s: 24

## Next
Stop pushing CpG. Try cCRE class rebalance — drop dELS (huge, generic
distal enhancer pool) and boost rare classes (TF, CA-TF, PLS).
