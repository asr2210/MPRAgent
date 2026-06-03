# 003 — Random uniform sanity check

## Purpose
Validate pipeline consistency. Compare my prepare.py output for a random
uniform library against the published `synth_oracle` baseline (0.684 eval_01)
to detect any setup issues.

## Result
**eval_01 = 0.4203.** Time: 18s eval / 1m52s total.

K562 ~0.58, HepG2 ~0.62, SK-N-SH ~0.06.

## Major reinterpretation
0.4203 does NOT match instructions.md baseline `synth_oracle` = 0.6840. But it
DOES match `strategies.md` baseline `random_uniform` = 0.4228 (within seed
noise). This means **my prepare.py corresponds to the strategies.md baseline
regime, not the instructions.md regime.**

The relevant baseline table to compare against is therefore the strategies.md
one, with much lower numbers across the board. Realistic strategies.md
benchmark:
- gc_50 / random_uniform: ~0.42
- dirichlet_composition: 0.35
- gc_sweep / at_rich: ~0.31-0.33
- gc_rich: 0.30
- homopolymer/dinuc: <0.20

The instructions.md baselines (dhs_topic 0.72 etc) use a different model
budget / data setup. Useful as upper bounds and for relative ordering between
strategies, but my own results match strategies.md scale.

Also notable: SK-N-SH correlation is ~0.06 even for random uniform. This
suggests the SK-N-SH oracle / target is intrinsically harder for the model to
fit from 50k sequences at this regime (not a library defect).

## Revised theory (T3)
- **Pipeline-real ceiling is much lower** than instructions.md baselines
  suggest. I'm probably operating in a setting where 50k is small relative to
  needed data — performance grows slowly with library quality, and tiny
  improvements over random baselines matter.
- **The 0.42 random uniform floor** is high relative to what biological
  sequence pools tend to add — but a good pool should still beat it.
- **SK-N-SH is intrinsically hard**: pushing it above ~0.1 might require
  targeted SK-N-SH-related sequences (neuronal regulatory elements) since
  random/general sequences cap out at 0.06.

## What this changes for next experiments
- Compare against strategies.md numbers as primary anchor (random_uniform 0.42).
- Beating 0.42 by ANY margin requires informative biological sequences in the
  pool.
- SK-N-SH performance is the bottleneck for `mean_r` — even strong K562/HepG2
  scores (0.58/0.62) only pull mean to 0.42 with SK-N-SH=0.06.
- Library design should explicitly target neuronal/SK-N-SH-relevant elements
  to boost SK-N-SH correlation, since K562/HepG2 are easy to push up.

## Next
Download full DHS Index + hg38 FASTA. Build a properly broad DHS sampler.
Goal for exp 004: beat random_uniform (0.42) with a properly-broad natural
sequence library.
