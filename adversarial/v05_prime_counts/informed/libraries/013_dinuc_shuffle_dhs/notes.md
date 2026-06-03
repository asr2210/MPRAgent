# 013_dinuc_shuffle_dhs — notes

## Design
Take dhs_specific (the best library, 0.049) and apply Altschul-Erickson
dinucleotide-preserving shuffle to each sequence. Preserves per-sequence
GC + dinucleotide frequencies. Destroys motifs/repeats/conservation.

## Result
eval_01 = 0.0376  (vs dhs_specific 0.049) — DROP of 0.012
eval_08 = 0.0570  (vs dhs_specific 0.055) — slight gain
mean_r = 0.0367

## Interpretation
Motif content explains ~0.012 of the eval_01 score for dhs_specific. The
remaining ~0.038 comes from dinucleotide composition / GC.

The eval IS sensitive to motif content (orientation, position), not just
composition. This means motif-augmenting the library should be productive.

eval_08's slight INCREASE on shuffle confirms eval_08 prefers more random-
looking content over real motifs.
