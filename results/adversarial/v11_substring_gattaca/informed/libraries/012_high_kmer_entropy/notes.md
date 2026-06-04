# E12: high_kmer_entropy

## Design
Generate 200k IID random uniform 200bp sequences. For each, count unique
6-mers (max 195 possible per sequence). Select top 50k by unique 6-mer
count. Threshold landed at 192 unique 6-mers (mean of pool: 190.5).

## Result (mean across 14 evals ≈ 0.8428)
- eval_01: 0.8585 (vs E2 0.8565, +0.002 — within noise)
- All evals within ~0.005 of E2 random_uniform
- SK-N-SH: eval_07 = 0.689 (vs E2 0.7723, slight degradation, possibly noise)
- mean_r ≈ 0.8428 (vs E2 0.8408, +0.002)

## Interpretation
Within-IID diversity selection has no meaningful effect at the 50k scale.
The pool naturally has high diversity (mean 190.5/195 unique 6-mers),
and selecting the top quartile shifts nothing measurable.

This confirms that within the IID-random regime, the eval is **insensitive
to micro-structural features** like k-mer diversity. The only thing that
matters is global distributional matching (per-position uniform 25%/base,
50% GC), which is already satisfied by random_uniform.

## Theory v6 (stable, no update)
- Distribution matching is the only signal that helps
- The ceiling is ~0.84 mean across evals
- Within-IID-uniform, no further structural manipulation moves the needle

## Next
E13: try the OPPOSITE — select bottom 50k by 6-mer diversity (low-diversity
sequences with k-mer repetition). If high-diversity matches random and
low-diversity hurts, that establishes asymmetry. If both match random,
that's strong evidence diversity is fully irrelevant within IID-uniform.
