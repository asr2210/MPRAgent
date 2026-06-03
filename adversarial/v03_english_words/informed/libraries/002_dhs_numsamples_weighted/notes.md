# 002 — DHS weighted by numsamples

## Plan
Test T1 hypothesis: by upweighting elements active in many biosamples
(numsamples-weighted sampling), shift the 160k synthseqs pool toward
constitutive elements that should be active in K562/HepG2/SK-N-SH and thereby
provide measurable signal for training.

Sampled with weight ∝ numsamples (linear). Resulting library: mean numsamples
= 33.3 (vs full-pool mean 15), max 401. Mean NMF proportion still 0.825.

## Result
**eval_01 = 0.3901.** Essentially identical to exp 001 (0.3883).
SK-N-SH still ~0.07 across all evals.

## What this teaches
T1 hypothesis is REFUTED. Even shifting the synthseqs pool toward
constitutive elements does not recover performance. The synthseqs subset is
fundamentally limited:
- Mean proportion is still 0.825 even after numsamples reweighting — the pool
  was curated for cell-type specificity, so even its "constitutive-ish" elements
  are still relatively specific.
- The pool likely excludes most ubiquitously-open promoters/enhancers, which
  are the strongest signal-producing elements in any individual cell line.

**Updated finding**: the synthseqs `classifier_light` files are not a usable
foundation for MPRA library design, regardless of reweighting. We need a
broader source.

## Theory update
T2 (revised): library value comes from "elements whose activity is
*measurable* in our cell types AND that span the regulatory grammar." The
synthseqs pool fails on the first axis no matter how we reweight. Compare:
even uniform random sequences (synth_oracle baseline) achieve 0.684 eval_01 —
clearly better than my biased-toward-specificity natural sequences.

The instruction-given baselines all sample from the *full* 3.6M DHS index, not
this 160k curated subset. The full pool naturally contains many constitutive
elements that the subset omits.

## Next
Two paths:
A) **Sanity-check pipeline**: run a random-uniform library. Should hit ~0.68 if
   pipeline matches the baseline table. This proves whether my failures are
   real (pool problem) or a methodological issue.
B) **Switch sources**: download the full DHS Index TSV + hg38 FASTA, properly
   sample with topic-loading weights.

I'll do (A) first as exp 003 — it's instant to generate and proves whether my
prepare.py invocation is correct.
