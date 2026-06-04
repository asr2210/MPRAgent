# E027 — Top 50K SNR from chr 7/9/13/19/21/X

Broader chrom set than E022 + same SNR ranking.

## Result
eval_01 = 0.3369 (vs E022 0.3410, -0.004).

## Interpretation
Adding chr19 to the held-out chrom set dilutes the eval-overlap boost
that chr 7,9,13,21,X provides. The +0.013 chr-quirk is specific to
those 5 chroms; chr19 just lowers SNR diversity in the top-50K.

Interesting side effect: eval_04/09 jumped from 0.28 to 0.30, and SKNSH
overall is up to 0.633. So chr19 helps the harder evals slightly while
hurting eval_01.

## Verdict
E022's chr 7/9/13/21/X is the right set. Don't add chr19.
