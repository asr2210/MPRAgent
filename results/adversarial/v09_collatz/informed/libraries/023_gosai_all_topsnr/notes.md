# E023 — Top 50K Gosai by SNR (all chroms, no chr restriction)

Same SNR ranking as E022 but pulled from full 763K Gosai pool.

## Result
eval_01 = 0.3250. Per-cell K562=0.145, HepG2=0.201, SKNSH=0.629.

## Interpretation
SNR alone barely helps (+0.002 vs random Gosai 0.323). But notably SKNSH
jumped from 0.62 to 0.629 — the highest SKNSH seen so far. K562/HepG2
dropped slightly though, washing the mean out.

The E022 combo (chr + SNR) at 0.341 is genuine synergy: the chr-restriction
keeps eval-overlap intact, and SNR pushes SKNSH up within that subset. Either
filter alone gives only small gains.
