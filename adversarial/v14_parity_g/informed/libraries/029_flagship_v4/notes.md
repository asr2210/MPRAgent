# 029_flagship_v4

**Final flagship principled library.**

Composition (50k total, every component RC-augmented):
| component | parents | + RC | total |
| --- | --- | --- | --- |
| cCRE-PLS               | 3,000 | 3,000 | 6,000  |
| cCRE-dELS              | 3,000 | 3,000 | 6,000  |
| cCRE-pELS              | 1,500 | 1,500 | 3,000  |
| cCRE-CTCF-only         | 1,250 | 1,250 | 2,500  |
| cCRE-DNase-H3K4me3     | 1,250 | 1,250 | 2,500  |
| Sharpr top by activity | 5,000 | 5,000 | 10,000 |
| Sharpr bot by activity | 5,000 | 5,000 | 10,000 |
| synthetic motif-packed (50-vocab) | 2,500 | 2,500 | 5,000 |
| motif-in-cCRE-backbone (50-vocab) | 2,500 | 2,500 | 5,000 |
| **TOTAL** | 25,000 | 25,000 | **50,000** |

Design rationale (after 28 prior experiments):
* The v14 metric is uninformative — all scores fall in [-0.005, +0.005]
  regardless of library. Seed-to-seed variance of the *same recipe* is
  ~0.001 (019 vs 022). Therefore optimize for **biological coverage**, not
  the metric.
* Real regulatory-element class breadth (5 cCRE classes) covers promoters,
  proximal/distal enhancers, CTCF insulators, broad-H3K4me3 marks — the
  major regulatory contexts a generalizing MPRA model needs to learn.
* Real measured MPRA activity poles (Sharpr top/bottom) provide anchored
  high/low activity exemplars from a published dataset.
* Synthetic motif-packed sequences (50-TF consensus vocabulary) give
  controlled combinatorial coverage of TF grammar above natural rates.
* Motif-in-cCRE hybrid combines authentic genomic context with controlled
  motif insertion.
* Reverse-complement augmentation of every component encourages strand-
  invariant motif representations.

Result: eval_01 = +0.0001. Noise, as expected.
30 (flagship v4 with SEED=1) = +0.0007 — seed-to-seed reproducibility = 0.0006.

This library is the deliverable: the most defensibly designed 50k 200bp
DNA library we could ship, given the constraint that the optimization
signal (v14 metric) is essentially uninformative.
