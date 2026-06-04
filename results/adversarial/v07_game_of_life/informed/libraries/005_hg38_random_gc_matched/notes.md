# Experiment 005 — Random hg38 genomic windows at GC ∈ [0.45, 0.55]

## Design
50k 200bp windows uniformly sampled from autosomes (chr1-22 + X + Y),
weighted by chromosome length, then filtered for all-ACGT and GC in
[0.45, 0.55]. Mostly intergenic/intronic — biologically "unselected" but
composition-matched to gc_50 baseline.

## Hypothesis
Tests whether the CURATION of natural sequence matters when composition is
controlled. If 005 ≈ 0.392 (= cCRE GC-matched), the curation status doesn't
matter; any natural human DNA at GC=50% gives the same result. If 005 >
0.397 (> gc_50), natural genomic k-mer structure helps even outside
regulatory regions.

## Result
- eval_01 = **0.3900** (vs cCRE GC-matched 0.3921, gc_50 0.397, cCRE bal 0.392)
- Mean across 14 evals: **0.3805**
- Per-cell-type on eval_01: K562 0.6013, HepG2 0.4232, SK-N-SH 0.1454
- Runtime: 1209s

## Interpretation
Random hg38 windows at matched GC score **the same** as cCRE-derived sequences
at matched GC. **The curation status of natural sequence is irrelevant** at this
50k scale — only composition (and maybe k-mer structure shared by all natural
human DNA) matters.

This is a fifth datapoint cementing theory v3:

| library              | mean GC | eval_01 | source                  |
|----------------------|---------|---------|-------------------------|
| gc_50 (baseline)     | 0.50    | 0.397   | i.i.d. random           |
| 001 cCRE balanced    | 0.51    | 0.392   | cCREs, all categories   |
| 002 motif+gc50       | 0.50    | 0.394   | random + JASPAR implants|
| 003 PLS only         | 0.62    | 0.375   | cCRE PLS only (bad GC)  |
| 004 cCRE GC-matched  | 0.50    | 0.392   | cCREs filtered to 0.5GC |
| 005 hg38 random      | 0.49    | 0.390   | random genomic windows  |

## What this rules out
- Natural-genome dinucleotide bias does not help over i.i.d. uniform random
- Annotation/curation status does not matter
- Biological provenance (regulatory vs intergenic) is invisible at this scale

## What's left to try
- **Activity-range coverage**: use a pretrained oracle (e.g., Malinois) to
  pre-select sequences that span the full predicted activity distribution.
  Hypothesis: more dynamic-range gives the small model more discriminative
  signal per sequence.
- **Cell-type-specific motif clustering**: not just single motifs but motif
  pairs / triplets with appropriate spacing (real enhancer grammar).
- **Adversarial / max-information sequences**: use the trained model on a
  prior library to find sequences that maximize prediction uncertainty.

## Key takeaway
After 5 experiments, I have a strongly-supported single-variable model: the
expected score on this evaluator is almost entirely predicted by mean GC
content. To move beyond ~0.40, I need an entirely different mechanism — most
likely **activity-range selection** via an oracle.
