# 015_chr19_22_only

## Design
50k random 200bp tiles from chr19 + chr22 (most gene-dense human
chromosomes: chr19 = 23 genes/Mb, chr22 = 14.5 genes/Mb). Length-
weighted; 53% chr19 / 47% chr22. seed=0.

## Result
- eval_01 = **0.4902** (vs E2 0.4992; Δ -0.009)
- mean over 14 ≈ 0.489

Per-eval: eval_04 and eval_06 slightly UP (+0.003 to +0.006);
eval_07 and eval_13 DOWN (-0.02). Aggregate slightly down.

## Interpretation
Maximizing gene-density doesn't beat broad gene-rich. The
narrowness penalty (less sequence diversity, more tile overlap at
10Mb sampling from 111Mb pool) offsets the gene-density gain.

This confirms: within gene-rich human DNA, the gene-density axis
is SATURATED at chr1/17/19/22. Further narrowing buys nothing.

## Theory update v6.1
The 0.50 ceiling for 50k libraries is robust. Tested negative for:
- Different DHS strategies (E3 stratified, E10 uniform, E13 top-signal)
- Cross-species mixing (E12)
- Tightening to top-density chroms (E15)
- Adding ectopic motifs (E7, E8)

The "narrowness penalty" axis is well-characterized. The "diversity
gain" axis has been weakly positive (E12: +0.01 from mouse, hidden
under -0.025 dilution).

## What's left to try

I'm running out of easy "add more curation" experiments. The path
forward requires a qualitatively different approach. Options:

1. **Element-class mix from cCRE**: explicit balance across promoter,
   enhancer, CTCF, intergenic classes. Different cut than DHS.
2. **Small-fraction admixture**: 10% mouse, 10% synth, etc. to see
   if dilution penalty scales sub-linearly with fraction.
3. **Quality-filtered random**: reject low-complexity tiles.
4. **Conservation-stratified**: use mouse-syntenic regions as proxy
   for vertebrate conservation. Requires alignment data I don't have.
5. **Activity-stratified**: train a quick predictor on DHS mean_signal
   then sample for activity diversity.

## Plan for E16
Try the cCRE element-class balanced mix. Requires downloading ENCODE
cCRE BED (~5MB) and sampling 12.5k from each of {PLS, pELS, dELS,
CTCF-only}. This is a fundamentally different stratification than
NMF or mean_signal — it's based on ChIP-seq histone marks at
elements.

If this lands at 0.49-0.50: element-class balance is also saturated.
Time to consider radical designs (synth + human hybrids, conservation,
or accepting the ceiling).
