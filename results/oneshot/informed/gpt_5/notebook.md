# MPRA library design notebook

## Initial interpretation

The objective is a 50,000 sequence, 200 bp MPRA library for training a model of
general regulatory grammar, not a model specialized to K562, HepG2, or SK-N-SH.
The most informative prior result is that DHS-topic sampling was best at 50k on
the first evaluation set and remained very strong across nearly all anonymous
sets. This says that real accessible regulatory DNA is the safest high-density
training substrate. Fully random sequence was surprisingly competitive on some
sets, especially the hardest sequence-space-diversity-looking set, but it gave up
too much on likely genomic regulatory evaluations.

My working theory is that a strong one-shot library should be mostly real
regulatory sequence sampled broadly across cCRE classes, with a deliberate
minority of synthetic perturbation/control sequences. The real component should
cover promoters, enhancer-like elements, CTCF/architectural elements, accessible
TF-bound regions, and accessible-only regions. The synthetic component should
not be pure iid noise only; it should include motif-density, GC, spacing, and
orientation regimes that are hard to see uniformly in a finite genomic sample.

## Data source decision

The instructions say I do not have access to the baseline DHS/SEI pools. For a
public substitute, I chose ENCODE SCREEN Registry v4 human cCRE class BED files.
This is not exactly the DHS pool used in the baselines, but it captures a broad
union of candidate regulatory elements across many cell and tissue contexts. The
available class files include promoter-like, enhancer-like, CTCF-bound and other
accessible/TF classes, which is well aligned with the goal of general regulatory
grammar.

I will fetch centered 200 bp windows from hg38 for sampled cCREs. For elements
shorter than 200 bp, centering around the annotated midpoint should retain the
core signal plus flanking context. For longer elements, the midpoint is a simple
class-agnostic anchor; adding small random jitter around the midpoint avoids
training only on one exact annotation convention.

## Mixture plan before implementation

I will not try to optimize for one named cell line. The planned mixture is:

- 38,000 real cCRE-centered genomic windows, class-balanced enough to avoid
  being dominated by distal enhancers while still respecting that enhancers are
  the largest source of diverse grammar.
- 4,000 genomic negative/background windows sampled from non-N bases away from
  explicit cCRE windows, to help the model learn what inactive-ish sequence looks
  like in a natural genomic distribution.
- 5,000 motif-synthetic sequences with combinations of common regulatory motifs
  placed in varied GC backgrounds, spacings, orientations, and copy numbers.
- 3,000 dinucleotide-shuffled cCRE-derived controls, preserving local composition
  while disrupting motif syntax.

This is intentionally more conservative than a 1/3 synthetic mixture because
the prior 50k table shows real DHS-heavy designs dominate at this size. The
synthetic/control part is large enough to increase grammar coverage but small
enough not to replace the strongest signal source.

## Implementation notes

I downloaded the ENCODE SCREEN Registry v4 GRCh38 cCRE BED file and the UCSC
hg38 2bit genome. The cCRE class counts are very uneven: dELS is much larger
than the other classes, while CA-TF is small. I decided against sampling
proportional to raw class frequency because that would make the library mostly
distal enhancers and would underrepresent promoter and architectural syntax.

The implemented generator uses:

- 42,000 cCRE-centered genomic windows:
  - 13,000 dELS
  - 7,500 pELS
  - 6,000 PLS
  - 4,500 CA-CTCF
  - 4,000 CA
  - 3,000 CA-H3K4me3
  - 2,500 TF
  - 1,500 CA-TF
- 3,000 natural genomic background windows avoiding annotated cCRE intervals.
- 3,000 dinucleotide-shuffled cCRE-derived controls.
- 2,000 synthetic motif-combination sequences.

The cCRE windows use the midpoint plus up to 40 bp of deterministic random
jitter and a random reverse-complement orientation. The jitter is meant to avoid
overfitting to exact annotation centering while still preserving the cCRE core.
The reverse-complement step increases orientation invariance without changing
the underlying grammar.

For synthetic sequences, I used small motif programs rather than arbitrary
consensus piles: AP-1/ETS/GATA, SP1/NRF1/E-box, CTCF pairs, NFY/SP1/ETS,
FOX/SOX/RUNX, CREB/AP-1/IRF, p53/E-box/SP1, and MEF2/GATA/AP-1. Background GC
is drawn from a broad beta distribution and motifs are lightly mutated. This is
not meant to be biologically exhaustive; it is a compact way to expose the
oracle/model to spacing, copy number, GC, and orientation regimes that may be
sparse in 42k real windows.

The generator includes a pure-Python 2bit reader. I initially tried py2bit, but
it could not build in this environment because Python headers are unavailable.
Avoiding that dependency makes the generator reproducible under the system
Python.

## Pre-submission validation

I ran `python3 library/generate.py` and independently checked the output before
calling the evaluator. `library/sequences.txt` has exactly 50,000 lines, every
line is exactly 200 bases, every character is in A/C/G/T, and there are no
duplicate sequences. The aggregate GC fraction is 0.4778, with a wide per-sequence
range from 0.055 to 0.925 because the library intentionally mixes real genomic
windows, genomic background, shuffled controls, and synthetic motif-rich designs.

At this point the design is fixed. The next command is the one allowed
`prepare.py` evaluation.

## Evaluation result

I ran `python3 prepare.py library/sequences.txt` once. It evaluated the 50,000
sequences against the 14 anonymous test sets in 892.1 seconds and wrote
`library/result.json`.

The eval mean_r values were:

- eval_01: 0.7199
- eval_02: 0.8147
- eval_03: 0.7981
- eval_04: 0.7834
- eval_05: 0.7199
- eval_06: 0.8146
- eval_07: 0.7561
- eval_08: 0.6665
- eval_09: 0.8531
- eval_10: 0.7756
- eval_11: 0.7075
- eval_12: 0.6848
- eval_13: 0.7547
- eval_14: 0.8149

The average of the 14 reported mean_r values is 0.7617. This is close to, but
slightly below, the prior DHS-topic result on eval_01 (0.7232 vs. this library's
0.7199). It appears that ENCODE cCREs were a good substitute for DHSs, but the
extra background/shuffle/synthetic portion likely cost performance on eval_08
relative to the random-heavy baselines and did not fully recover the DHS-topic
advantage elsewhere. If I had another shot, I would likely reduce shuffled
controls, increase real cCRE/DHS-like active windows, and add a stronger
accessibility-signal proxy rather than class balance alone.
