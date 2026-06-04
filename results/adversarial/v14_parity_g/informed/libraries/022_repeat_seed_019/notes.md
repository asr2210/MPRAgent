# 022_repeat_seed_019

Identical recipe to 019 (50-vocab motif + RC) but SEED=1. Tests whether
019's apparent +0.0025 was reproducible or a noise draw.

Result: eval_01 = +0.0013 (vs 019's +0.0025).

**Conclusion**: same recipe, different seed → score moves by 0.0012.
Confirms that 019's "best" score was sampling variation, not signal.
The metric's seed-to-seed noise alone is ~0.001, comparable to the
inter-strategy variation across all 30 experiments. v14 is uninformative.
