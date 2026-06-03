# 030 — gene-density² × repeat-rich filter (stacking)

## Result
**eval_01 = 0.5064** (3-seed). Slightly below E28 (0.5086).

## Design
Combine the two confirmed positive signals:
  (a) gene-density² weighting at 500kb (E28 = 0.5086)
  (b) repeat-rich content filter (E18 = 0.4937)

Acceptance rate ~75% (50k from ~66k attempts).

## Interpretation
Stacking is **mildly negative**. Removing the 25% high-complexity tiles
from gene-dense regions slightly hurts.

Likely explanation: gene-dense regions (especially promoter-proximal,
CpG islands) naturally include high-complexity sequences (TF binding
sites with diverse motifs). The repeat-rich filter that helped on a
chr1/17/19/22 random pool (E18 vs E2) actually removes useful signal
when the pool is already gene-density-biased.

The two signals overlap rather than stack additively.

## Final result of the campaign
**E28 (gene-density² at 500kb, 3-seed) = 0.5086 is the campaign best.**

## Bin/exponent/EPS surface (squared, full genome, 3-seed)
| Config                       | eval_01 |
|------------------------------|---------|
| 1Mb,   EXP=2.0, EPS=0.5      | 0.5071  |
| 500kb, EXP=2.0, EPS=0.5      | **0.5086** ← best |
| 500kb, EXP=2.0, EPS=0.1      | 0.5083  |
| 250kb, EXP=2.0, EPS=0.5      | 0.5084  |
| 100kb, EXP=2.0, EPS=0.5      | 0.5063  |
| 1Mb,   EXP=1.0, EPS=0.5      | 0.5023  |
| 1Mb,   EXP=3.0, EPS=0.5      | 0.5030  |

Saturated at ~0.508 across the wide plateau (250-500kb, EXP~2, EPS<=0.5).
