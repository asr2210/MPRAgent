# 005_markov_k4

## Design
50k 200bp sequences from a 4th-order Markov chain trained on
hg38 chr1+chr17+chr19+chr22 (same chromosomes as 002 real genome).
Laplace smoothing on 256 contexts × 4 outputs. seed=0.

Predicted to isolate the value of 4-mer composition. If
matching k=4 reproduces 002's score (0.499), real-DNA prior =
k-mer composition. If much worse, higher-order structure matters.

## Result — strong negative
- eval_01 = **0.2681** (vs synth random 001: 0.3068; vs real 002: 0.4992)
- mean over 14 ≈ 0.267
- Markov k=4 is WORSE than uniform random.
- eval_07 = 0.3902, eval_13 = 0.3599 — distinct evals also worse than synth.

## What this means
4-mer composition matching is NOT sufficient to capture the real-DNA
prior. Stronger result: it's actually counter-productive. Possible
mechanisms:
1. Markov k=4 occasionally produces low-complexity stretches (long AT
   or CG runs that match training corpus stats). The model trained on
   these confuses normal genomic patterns with extreme ones.
2. The Markov model averages over chromosome-specific biases (chr19
   GC-rich, chr22 CG-island-rich) producing a chimeric composition
   that no real region matches.
3. Sequences that "look like real DNA compositionally" but lack motif
   structure may produce specific incorrect activity predictions — the
   model overfits to spurious local patterns. By contrast, uniform
   random is obviously noise; the model defaults to near-flat output
   which correlates trivially with the eval mean.

Whatever the mechanism, the **empirical conclusion is firm**: matching
4-mer composition does not produce the real-DNA effect. The +0.19
r-point gap from 001→002 lives in higher-order structure: motifs,
repeats, gene/intergenic structure, or specific functional sequences.

## Theory update
v0 said: real DNA's value is its "functional grammar density". I had
no decomposition. Now I do:
- ~0% of the real-DNA effect is captured by 4-mer composition.
- The effect must live in {motif content, long-range structure,
  specific functional sequences, repeat structure}.

## Implications for design
- Synthetic-with-k-mer-matching is NOT a viable substitute for real
  DNA. Drop this design direction.
- Either (a) use real DNA as the backbone, with refinements, or (b)
  build synthetic sequences with EXPLICIT motif content.
- A motif-implanted-scaffold design becomes much more attractive:
  random sequence + planted TF motifs at known positions and
  frequencies should approximate the "motif content" component of
  real DNA's value.
