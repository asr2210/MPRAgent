# 003_encode_ccre

## Design
50k 200bp windows centered on midpoints of randomly-sampled ENCODE v4
cCREs (~2.35M total, classed as dELS/pELS/PLS/CA/CA-CTCF/CA-H3K4me3/TF/CA-TF).

## Hypothesis
cCREs are a different annotation than DHS — typed by function (enhancer/
promoter/CTCF). If the eval was built using cCRE-style elements (as Agarwal
et al. 2025 lentiMPRA does), this should beat DHS-random (0.041). If not,
all standard genomic libraries hit the same ~0.04 ceiling.

## Result
**eval_01 = 0.0466** — slightly above cCRE > DHS > random by ~0.005 each.
Still near the universal ~0.04 ceiling.

eval_08 = 0.055, vs DHS 0.045 and random_uniform 0.121 in strategies.md.
cCRE is BETTER than DHS on eval_08 but much worse than random.

## Interpretation
- Tiny ordering: cCRE > DHS > random by ~0.005. Effect is small but
  consistent across most evals.
- The "ceiling" of ~0.04-0.05 holds for all annotation choices tested so
  far. The eval doesn't reward any specific genomic regulatory annotation.
- eval_08 specifically favors unstructured/random sequences. Bio sequences
  cost ~0.07 r on eval_08 alone.

## What to try next
- Mix random + cCRE (test if hybrid helps)
- Saturation mutagenesis-style libraries (variation within small set of seeds)
- Pure "max compositional diversity" libraries
- Eventually: explore whether eval rewards a very specific design (e.g.,
  positional motif placement, cell-type-specific TF clusters)
