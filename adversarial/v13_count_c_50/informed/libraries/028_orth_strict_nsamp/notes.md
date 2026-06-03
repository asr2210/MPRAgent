# Experiment 028 — orth-DHS strict numsamples (≤2) + 025's cCRE

## Result
eval_01 = 0.5755 (−0.0007 vs 025's 0.5762).

Tightening numsamples from ≤5 to ≤2 hurt slightly. The ≤5 cap was
already capturing the cell-type-specific tail; pushing tighter strips
out informative low-mid-frequency elements and shrinks the eligible
pool.

## Numbers
| eval | 020 (≤5, no-prune) | 025 (≤5, no-dELS-no-CA) | **028 (≤2, no-dELS-no-CA)** |
|------|-------------------:|------------------------:|----------------------------:|
| 01   | 0.5745             | 0.5762                  | 0.5755                      |
| 04/09| 0.5509             | 0.5620                  | 0.5572                      |
| 07   | 0.6168             | 0.6131                  | 0.6165                      |
| 08   | 0.1404             | 0.1591                  | 0.1452                      |
| 13   | 0.5968             | 0.5925                  | 0.5956                      |

## Conclusion
**≤5 is the sweet spot for numsamples.** Two effects offset:
- Stricter nsamp lifts eval_07/13 (regulatory grammar tasks) by ~+0.003
- But hurts eval_04/09 (cell-type-specific) by −0.005

Net eval_01: −0.001. The medium-specificity DHS (numsamples 3–5)
carries enough cell-type signal to matter; cutting them was over-
selecting.

## Pool sizes
nsamp≤2 pool: ~340k orth-DHS (vs ~610k at nsamp≤5).
Still 13x N_DHS=25k, so coverage isn't the issue — quality is.

## Next
Specificity lever exhausted at the numsamples knob. Try complementary
levers: 029 = small CpG slice on top of 025; 030 = final consolidated
best design.
