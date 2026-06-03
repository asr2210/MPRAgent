# 017_dhs_multiwindow — notes

**Design**: 25k unique DHSs from non-label-aligned pool, each emitted as two windows at offsets (-50, +50) from summit. 50k total from 25k biological sources.

**Result**: eval_01 = 0.6715 (vs exp 005 (50k unique) = 0.6752, **-0.004**). Mean = 0.6266 vs 0.6282. Slight hurt.

**Conclusion**: source count beats per-source augmentation. 50k unique DHSs > 25k DHSs × 2 paired views. Cutting unique sources in half — even with augmentation back to 50k — loses biological diversity that the model uses.

**H8**: Biological source count is the bottleneck. Augmentation can't substitute for unique DHSs. With the library budget fixed at 50k, maximize unique biological sources within a thoughtfully filtered pool.

Combined with prior negative aug results:
- exp 009 (single random offset): flat
- exp 012 (RC aug × 2 strands, 25k sources): flat
- exp 017 (paired offsets, 25k sources): slight hurt

All augmentation strategies are at best neutral when they preserve source count, at worst hurtful when they cut sources. Conclusion: augmentation is the wrong lever in this regime.
