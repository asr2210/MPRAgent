# 008_multisource_s2_heavy — Table_S2 tripled to 15k

**Composition (50k):**
- 5k cCRE uniform
- 3k cCRE CTCF-only
- 3k cCRE DNase-H3K4me3
- 15k Table_S2 UKBB+GTEx (3x exp 007)
- 4k DHS-topic
- 20k paired flanks
- (synthetic dropped vs exp 007)

**Result:** mean_r(14) = 0.1503  eval_01 = 0.1569

**Vs exp 007 (mean_r=0.1532, eval_01=0.1592):**
| eval | 007 | 008 | Δ |
|---|---|---|---|
| eval_07 (SEI) | 0.131 | 0.166 | **+0.035** |
| eval_06 (K562) | 0.228 | 0.187 | -0.041 |
| eval_11 (K562) | 0.228 | 0.187 | -0.041 |
| eval_01 (chr-GT) | 0.159 | 0.157 | -0.002 |
| eval_13 (genomic) | 0.127 | 0.148 | +0.021 |
| eval_08 (synth) | 0.046 | 0.041 | -0.005 |
| eval_10 (DHS) | 0.121 | 0.132 | +0.011 |

**Reading:**
- Table_S2 dose-response: more Table_S2 → eval_07 (SEI chr7/13) rises.
- BUT it directly competes with the cCRE volume that drives K562 unlock on eval_06/11.
- Net mean_r drops slightly (-0.003) because eval_06/11 are equally weighted.
- Synth removal: trivial cost (-0.005 on eval_08, which is already noise-floor 0.04).
- eval_13 improves: genomic chr7/13 benefits from more genomic-distributed Table_S2 seqs.

**Takeaway:** The 5k-S2 dose in exp 007 is near-optimal for the eval mix. Doubling S2 only helps if SEI weight matters more than K562. For mean_r maximization, exp 007's recipe wins.

**Implication for next:** Don't keep scaling Table_S2. Instead, ATTACK the K562 unlock mechanism directly — what about exp 007 (which had 4k+4k cCRE CTCF/DNH3) drove eval_06/11 to 0.228? Try doubling CTCF specifically.
