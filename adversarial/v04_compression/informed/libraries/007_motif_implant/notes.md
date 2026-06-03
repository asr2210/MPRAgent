# 007_motif_implant

## Design
50k 200bp uniform random scaffolds, each with K=3 non-overlapping
TF consensus motifs planted at random positions (random strand).
Motif library: 53 curated strong vertebrate TFs (4–15bp consensus).
seed=0.

## Result
- eval_01 = **0.3468** (vs synth 001: 0.3068; vs real 002: 0.4992)
- mean over 14 ≈ 0.348
- eval_07 = 0.4448 (vs synth 0.4024 — biggest gain on this eval)
- eval_13 = 0.4272 (vs synth 0.3809)
- eval_08 = 0.0879 (worse than synth 0.1098)

## Interpretation
Planting 3 motifs in random scaffolds gives a +0.04 r-point bump on
eval_01 over uniform random — small. The motifs explain at most ~20%
of the real-DNA prior (0.04 of the 0.19 gap to real genome).

Combined with the Markov k=4 result (0.27, worse than synth), this
means real DNA's value is dominated by structures that are NEITHER
4-mer composition NOR isolated TF motifs:
- Likely repeat content (LINE/SINE/ALU, simple repeats)
- Long-range context (gene bodies, promoters, intronic vs intergenic)
- Motif co-occurrence and spacing (regulatory grammar)
- Composite elements (multi-TF CREs)

The model wants to see *natural sequence context*, not isolated
"signal-bearing" motifs.

## Useful signal: eval_07 and eval_13
These two evals were most responsive to motif planting. They might
be testing TF-motif-style sequences and reward motif content
directly. Watch them on future motif-augmented designs.

eval_08 keeps drifting around 0.09. It's not responsive to any
sequence-design choice tested so far. May test something fundamentally
different (designed/synthetic sequences, single-mutant scans, etc.).

## Generalization argument (revised)
Original argument: TF motifs are universal across cell types, so
planting them should give cross-cell-type generalization. Result:
they DO transfer some signal, but it's small.

The real DNA prior must include cell-type-spanning info beyond just
TF motifs. Probably: the natural co-occurrence statistics of motifs
in their gene-regulatory context. A library that captures this needs
to either (a) use real DNA directly, or (b) engineer composite
multi-motif elements in realistic context.

## Plan for next experiment
Test: hybrid design where each sequence is real-genome BASED but with
1-2 extra strong motifs PLANTED. Tests whether real-DNA prior PLUS
explicit motif signal beats pure real DNA.

Or alternatively: test cCRE-curated regulatory elements as the source.
cCRE includes promoters + enhancers + CTCF, providing functional
diversity that DHS alone may have missed.

I'll go for cCRE: it tests a well-curated functional element catalog
that the published baseline didn't compare. If cCRE matches or beats
real genome, functional curation does add value.
