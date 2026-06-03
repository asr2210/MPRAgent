# 018 — Random hg38 200bp regions

## What I built
Sampled 50k random 200bp windows from hg38 autosomes + chrX (~3Gbp total). Rejected
windows containing N. Sampled 50k from 51,962 attempts (3.8% N-rejection rate).

## Result
- eval_01 = 0.0735 (vs DHS NMF-strat 0.0739, vs Dirichlet(0.3) 0.0786). **≈ DHS.**
- mean ≈ 0.0903 (vs DHS 0.0901). Essentially tied with DHS.
- eval_08 = 0.0640 (vs DHS 0.0693). Slightly worse on noisiest eval.

## Interpretation
Random hg38 ≈ DHS. Generic genomic sequences and regulatory-enriched DHS perform
identically. This means:
- The regulatory enrichment in DHS is NOT what limits biology performance.
- Both genomic sources have the same compositional concentration (~50% GC) and the
  same overall composition distribution at scale.
- The biology ceiling at ~0.074 is set by **natural genomic composition concentration**,
  not by what subset of the genome we sample.

To break the 0.078 ceiling, we'd need a fundamentally different lever — not biology.

## Hypothesis killed
H18: random hg38 > DHS → REJECTED (tied, very slightly worse on eval_08).

## What to try next
Test if more granular 4D composition stratification (each base count axis) helps over
natural Dirichlet(0.3) draws. If it doesn't, the Dirichlet(0.3) ceiling is definitive.
