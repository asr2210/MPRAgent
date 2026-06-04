# 010_ccre_dELS_uniform — notes

**Design**: Orthogonal annotation source. 50k uniform sample from 1.47M ENCODE cCRE Registry-V4 distal-enhancer-like-signature (dELS) regions, centered on midpoint, ±100bp.

**Result**: eval_01 = 0.6671 (vs exp 005 DHS-no-label 0.6752, exp 003 DHS-uniform 0.6604). Mean across 14 evals = 0.6275 (vs exp 005 = 0.6282).

**Conclusion**: cCRE-dELS and DHS-no-label-topic train essentially equivalent models. cCREs are curated for enhancer-like signature (DNase + H3K27ac + H3K4me1) but the additional curation doesn't translate to a measurable lift in this pipeline. Both annotation sources are sampling the same underlying biology: DNase-accessible regulatory elements.

**Implication**: switching annotation source (DHS → cCRE) is not the right axis. The bottleneck is not "which catalog you pull from" — it's something more fundamental about either (a) the eval-set ceiling or (b) some selection axis not yet tried.

**Theory update** (still H4 — diversity wins):
> cCREs and DHSs are largely overlapping in biological content. cCRE classification adds functional labels (enhancer/promoter/CTCF) but doesn't change the underlying sequence pool meaningfully when sampling 50k from millions.

**Next ideas**:
1. **Mix cCRE classes** (dELS + pELS + PLS + CTCF + CA) — adds promoter/CTCF context. Tests if class diversity within cCREs matters.
2. **Motif-rich DHS selection**: PWM-scan DHSs for JASPAR motif content, sample top-scored. Tests "regulatory grammar density" hypothesis directly.
3. **Cross-source mix**: 25k DHS + 25k cCRE — but exp 006 (mixing hurt) suggests this might fail.
4. **GC-content stratification**: match eval-set GC distribution.
5. **Length variation**: try shorter/longer effective windows or random crop+pad.
