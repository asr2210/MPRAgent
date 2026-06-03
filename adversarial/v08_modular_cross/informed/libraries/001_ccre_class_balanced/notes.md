# Experiment 001 — cCRE class-balanced library

## Hypothesis
Real human cCRE regulatory sequences, stratified across 9 functional classes (dELS,
pELS, PLS, CTCF, etc.), should match or exceed the dhs_topic baseline (0.7232) from
instructions.md Table 1, because cCREs are the ENCODE V3 union across 1518
biosamples — broader coverage than the Meuleman 2020 DHS catalog (733 biosamples).

## Method
- Downloaded GRCh38-cCREs.bed (1.06M elements) and hg38.fa
- Stratified-sampled 50,000 cCREs (~5,555 per class)
- Extracted 200bp centered on each cCRE midpoint

## Result
**eval_01 = 0.0002** — essentially zero. All 14 eval sets score in the range
[-0.0075, +0.0051]. Time: 67s, single seed.

## Interpretation
This is **catastrophic and unexpected**. Real, biologically-curated regulatory
sequences should achieve >0.5 if the instructions.md baselines transferred to this
environment. Instead, they perform at chance level — same as random_uniform (0.0004)
from strategies.md.

Three possible explanations:
1. **strategies.md is the real baseline for v08.** All listed baselines score ~0 there.
   The instructions.md Table 1 baselines may be reference numbers from a different
   (earlier, non-adversarial) eval and don't apply here.
2. **My cCRE methodology subtly differs from dhs_topic.** dhs_topic uses NMF topic
   weighting on the Meuleman 2020 DHS catalog; I used class-balanced cCRE sampling.
   The detail of weighting may matter.
3. **OOD failure.** The eval sequences may be from a distribution very different from
   cCREs (e.g., designed/synthetic), so a model trained on cCREs can't generalize.

## Notes on eval structure
The eval scores show high redundancy:
- eval_01 ≈ eval_14 ≈ eval_06 ≈ eval_11
- eval_02 ≈ eval_05 ≈ exact match
- eval_03 ≈ eval_12 ≈ exact match
- eval_04 ≈ eval_09 ≈ exact match
- eval_07, eval_08, eval_10, eval_13 are distinct

Same redundancy pattern visible in strategies.md baseline rows. So the 14 sets
collapse to ~8-9 distinct evaluations.

## Updated theory
The v08 environment is genuinely hard. The straightforward "real biological sequences"
approach failed. Next: probe what kind of library escapes the chance-level plateau.
