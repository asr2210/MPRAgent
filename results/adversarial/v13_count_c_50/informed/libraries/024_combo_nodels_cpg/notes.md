# Experiment 024 — orth-DHS + cCRE-no-dELS + CpG supplement

## Result
eval_01 = 0.5745 (tied 020, **−0.001 vs 023's 0.5755**).

## Interpretation
**Adding CpG slice on top of no-dELS cCRE OVERSHOOTS.** The no-dELS
cCRE already upweights PLS+pELS to 3571 each. Adding 3k+3k more
brings totals to ~6571 each — like 022's 10k CpG slice, which also
hurt eval_01.

The CpG-rich content has diminishing/negative returns above ~7 % of
library. 021's 12 % slice on top of class-bal-with-dELS worked because
PLS+pELS was 3125 each (6250 total) + 6k more = 12250 ≈ 24 %. But
in that case the BASE was already lower CpG (dELS dilutes). With
no-dELS base + CpG slice, total CpG ratio jumps too high.

eval_08 recovers (0.184, like 021) but at slight eval_01 cost.

## Numbers
mean_r: 0.557
eval_01: 0.5745 (vs 023's 0.5755)
eval_08: 0.1840

## Conclusion
**023 design (orth-DHS + cCRE no-dELS, no CpG slice) remains best
for eval_01**. The CpG slice helps when base mix is light on CpG;
otherwise it's a wash or loss.

## Next
Try more aggressive cCRE rebalance: drop CA too. CA = DNase-only,
which is similar to DHS — possibly redundant given orth-DHS half.
