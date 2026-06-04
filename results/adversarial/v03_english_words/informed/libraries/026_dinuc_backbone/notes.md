# 026 — Mammalian dinucleotide backbone + 8 motifs overlap

## Plan
Backbone generated via Markov chain with hg38-like dinucleotide
frequencies (low CpG, ~41% GC). Insert 8 motifs (overlap) from 289 pool.
Tests whether bio-realistic backbone composition unlocks a new lever.

## Result
**eval_01 = 0.4153.** K562 0.578 (↓ from 0.599 in 020), HepG2 0.606
(↓ from 0.629), **SK-N-SH 0.062 (highest mean across all exps)**.

eval_07 **SK-N-SH = 0.0704** — breaks the 0.07 ceiling for the first time.

## What this teaches
- **Backbone composition is a real lever for SK-N-SH.** Going from
  uniform random → mammalian dinuc backbone shifts SK-N-SH from 0.057 →
  0.062 (and one eval to 0.070). The SK-N-SH oracle is responsive to
  background sequence statistics, not just motif content.
- K562/HepG2 PREFER higher-GC / more uniform backgrounds. With low-GC
  mammalian backbone they drop ~0.02 each. So backbone is also a lever
  for them, but in the opposite direction.
- The 008-family's strong K562/HepG2 performance was partly driven by
  random uniform backbone being well-matched to their oracle preferences.

## Theory T18
SK-N-SH oracle uses sequence COMPOSITION statistics (dinucleotide
frequencies, low-complexity patterns) more than discrete motif features.
K562/HepG2 oracles use motif features primarily. The plateau is partly
because uniform random optimizes K562/HepG2 backbone preference, not
SK-N-SH.

## Generalization implication
For an unknown cell type, library design should hedge across BACKBONE
composition too: some sequences random uniform (good for motif-focused
oracles), some sequences with bio-realistic dinuc backbone (good for
composition-sensitive oracles). This is a previously missed dimension.

## Next
Exp 027: test the OPPOSITE — HIGH GC backbone (60%) + 8 motifs. If K562
jumps further, GC tuning is a K562 lever. If SK-N-SH crashes, confirms
T18. This completes the backbone-composition picture.

Or: exp 027 = MIXED BACKBONE library (50% random uniform + 50% mammalian
dinuc), each with 8 motifs. Tests T18 via averaging — gain SK-N-SH +0.003
while losing only +0.01 on K562/HepG2 = net mean_r improvement.
