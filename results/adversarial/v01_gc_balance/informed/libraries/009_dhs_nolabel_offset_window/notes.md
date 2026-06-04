# 009_dhs_nolabel_offset_window — notes

**Design**: Same pool as exp 005 (12 non-label-aligned topics) but each window has a random offset of -50..+50 from summit (regulatory element can sit anywhere from position 50 to 150 within the 200bp window).

**Result**: eval_01 = 0.6724, essentially flat vs exp 005 (0.6752, -0.003). All evals nearly identical.

**Conclusion**: positional rigidity is NOT what's limiting performance. The model is already position-invariant (consistent with a convolutional architecture that's translation-equivariant by design). Positional augmentation is a non-issue.

**Implication**: stop tweaking the DHS sampling recipe — every variant lands ±0.01 of 0.675. The 0.05 gap to baseline `dhs_topic` (0.7232) is likely due to multi-seed averaging (baseline uses 5 seeds) and some protocol detail I haven't reproduced.

Pivot to a different axis: **motif/TF-content selection** of DHSs (orthogonal to topic / signal / breadth).
