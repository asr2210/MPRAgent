# E6: dinuc_markov_natural

## Design
Symmetric 1st-order Markov chain w/ CpG depletion (P(C→G)=P(G→C)=0.10
vs uniform 0.25) and mild AT depletion (P(A→T)=P(T→A)=0.20). Symmetric
→ uniform 50% GC stationary. 50,000 seqs from this chain.

## Result (mean across 14 evals ≈ 0.777)
- eval_01: 0.7966 (vs E2 0.8565, **-0.060**)
- K562: ~unchanged (~0.85)
- HepG2: slightly higher (~0.87)
- SK-N-SH: crashed (~0.65, eval_07 = 0.43)

## Interpretation
Even MILD natural-genome-like dinucleotide structure (CpG depletion +
AT depletion) hurts substantially. SK-N-SH again the most sensitive.

**Strong confirmation of theory v4**: the eval distribution is essentially
IID uniform random at 50% GC. ANY change to per-position statistics
(motifs, per-base balance, Markov dependencies) hurts SK-N-SH severely.

## What this tells us about generalization
SK-N-SH may be served by training sequences whose per-position statistics
match IID uniform PRECISELY. This is consistent with SK-N-SH's neural
biology: neural regulatory grammar may be more dependent on rare
combinations of features that only appear at the expected rate under
pure random. Any structural bias (motif planting, dinuc constraints,
balance) suppresses these rare combinations.

K562/HepG2 are more robust — perhaps because their regulatory grammar
has higher built-in redundancy or because the model can leverage even
slightly-biased data for them.

## Theory v4 (stable, confirmed)
The eval distribution ≈ IID uniform random 50% GC. The optimum library
is random_uniform. The ceiling is ~0.857 eval_01 / ~0.841 mean.

Beating this would require either:
1. Active selection within IID per-position distribution (e.g., based on
   model-predicted activity diversity) — limited by no oracle access.
2. Multi-seed aggregation — not allowed (1 file).
3. Something subtle that improves within-IID diversity.

## Pivot
E7 will attempt true_gc_50 (IID with rebalance to exact GC=100) to
properly replicate the published baseline. If it matches 0.859, I have
my new strongest baseline. Then E8 onward, I'll try diversity-based
selection experiments.
