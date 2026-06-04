# 015_dhs_celltype_plus_cancer — no lift

**Composition:** 7k Neural + 5k Myeloid + 5k Digestive + 5k Cancer/epi + 4k cCRE-uni + 4k S2 + 20k flanks.

**Result:** mean_r(14) = **0.1467**, eval_01 = 0.1504. *Slight regression vs recipe baseline.*

**Per-eval (selected):**
- eval_06/11 mean = 0.170 (vs baseline 0.190 multi-source seed avg)
- eval_07 = 0.188 (decent, S2 helped)
- eval_10 = 0.118 (regressed; less topic-DHS)
- K562 r near zero everywhere (no lift from cancer-epithelial)

**Conclusion:** Cancer / epithelial component doesn't differentially help K562/HepG2. Probably because that component is a chromatin-program signature, not a cell-line-specific one. The cell line being "cancer-derived" doesn't equate to "cancer/epithelial chromatin program" — K562 is myeloid leukemia, HepG2 is hepatocellular, both with very different chromatin landscapes from epithelial cancers.

**Next:** Stop tweaking DHS components. Switch to maximizing the dominant SK-N-SH r since it has the most headroom (current ~0.45, theoretical ceiling probably ~0.55-0.60).
