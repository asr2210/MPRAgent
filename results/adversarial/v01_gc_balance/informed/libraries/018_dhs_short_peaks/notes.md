# 018_dhs_short_peaks — notes

**Design**: filter DHSs to those with peak width ≤300bp (sharp peaks only). 2.33M of 2.61M no-label DHSs retained. Uniform sample 50k.

**Result**: eval_01 = 0.6685 (vs exp 005 = 0.6752, **-0.007**). Slight hurt.

**Conclusion**: broad peaks (>300bp) carry useful training signal — filtering them out removes biology the model uses. Not a productive filter axis.

**Bigger picture**: 11 filters/selections/augmentations have now been tested on top of exp 005 (uniform from 12-topic-excluded DHS pool). None beats the baseline. The 0.6752 plateau seems robust to:
- Annotation source switch (cCRE)
- Window anchor (summit vs core)
- Window offset / RC augmentation
- Multi-window paired views
- Weighting (forward, inverse numsamples)
- Selection (motif-density, GC-extreme)
- Filtering (peak-width)

**The 0.05 gap to the `dhs_topic` baseline (0.7232)** is most likely:
1. Multi-seed averaging (baseline = 5 seeds, mine = 1)
2. Per-DHS NMF loading matrix (not in my Index file; would need separate Zenodo download)
