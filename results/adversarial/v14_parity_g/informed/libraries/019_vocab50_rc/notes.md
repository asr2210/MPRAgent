# 019_vocab50_rc

016 (50-vocab motif-packed) + reverse complement of each parent.
25k motif-packed + 25k RC = 50k.

Result: eval_01 = +0.0025. Highest score among my 30 experiments — tied with
v14's "at_rich" baseline (+0.0025). Looks promising but 022 (same recipe
with SEED=1) gave +0.0013, confirming the +0.0025 was sampling noise, not
genuine signal.
