# 012_dhs_nolabel_rc_aug — notes

**Design**: 25k DHSs uniformly from non-label-aligned pool, each emitted as forward + reverse complement = 50k sequences from 25k biological sources.

**Result**: eval_01 = 0.6760 (vs exp 005 = 0.6752, +0.001). Mean = 0.6296 (vs 0.6282, +0.001). **Effectively flat — within seed noise.**

**Conclusion**: model already learns RC-equivariance internally (consistent with conv-net architecture). RC augmentation provides no measurable lift. Same flat outcome as positional jitter (exp 009).

**Pattern**: any augmentation that respects an invariance the model already has → flat. Augmentation only helps when it breaks a model bias.

**Implications**:
- Strand and small-position invariances are free.
- Need to find augmentations the model DOESN'T have implicitly.
- The 0.05 gap to baseline `dhs_topic` (0.7232) is not strand or position; likely multi-seed averaging.

**Next**: try genuinely new biological variation — DHS centered on `core_midpoint` instead of `summit` (paired window axis), or length-stratified DHS sampling (favor short/sharp peaks vs broad). Or counter-test: motif-SPARSE selection (exp 013 candidate) to confirm H5.
