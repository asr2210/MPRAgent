# Experiment 029 — 025 design + small CpG slice (3k)

## Result
eval_01 = 0.5740 (**−0.0022 vs 025's 0.5762**).
eval_08 = 0.1658 (+0.007 vs 025) — small CpG lift, as expected.

The 3k CpG slice steals from cCRE class balance (22k → 22k cCRE means
each of 6 classes drops from 4166 → 3666). The loss of regulatory
class diversity outweighs the eval_08 gain.

## Numbers
| eval | 025    | **029** | Δ      |
|------|-------:|--------:|-------:|
| 01   | 0.5762 | 0.5740  | −0.0022|
| 04/09| 0.5620 | 0.5611  | −0.001 |
| 07   | 0.6131 | 0.6109  | −0.002 |
| 08   | 0.1591 | 0.1658  | +0.007 |
| 13   | 0.5925 | 0.5904  | −0.002 |

CpG slice mean: 30.2 CpG/200bp (vs orth-DHS=1.9, cCRE=4.5).

## Conclusion
**Adding CpG on top of 025 is net-negative for eval_01**, same conclusion
as 021/022/024 found on earlier bases. Budget for CpG always comes at
the cost of class diversity. eval_08 is a small tax most designs pay.

The CpG slice does what it says (lifts eval_08), but eval_01 dominates
the overall metric.

## Final answer for eval_01
**025 remains the best design**: 25k orth-DHS (nsamp≤5) + 25k cCRE
balanced across 6 classes (no dELS, no CA).
