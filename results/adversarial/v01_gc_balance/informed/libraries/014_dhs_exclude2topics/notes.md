# 014_dhs_exclude2topics — notes

**Design**: Exclude only the 2 most label-specific NMF topics (Myeloid/erythroid, Neural). 14 topics retained. Uniform sample 50k.

**Result**: eval_01 = 0.6691 (vs exp 005 (4 topics excluded) = 0.6752, exp 003 (0 excluded) = 0.6604). Mean = 0.6253.

**Pattern (topic-exclusion dose-response)**:
| topics removed | eval_01 |
|---|---|
| 0 (all 16) | 0.6604 |
| 2 (Myel/Neu) | 0.6691 |
| **4 (4 most label-aligned)** | **0.6752** |
| 12 (only 4 retained, topic-balanced) | 0.6631 — note: also balanced not uniform |

Removing 4 is better than 2. Suggests Digestive + Cancer/epithelial were also confounding (HepG2 is HCC, hits both). 4-topic exclusion is a local optimum on this axis.
