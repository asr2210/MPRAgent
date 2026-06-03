# 002_genome_random

## Design
50,000 random 200bp windows from hg38 chr1, chr17, chr19, chr22.
Sampled with probability proportional to chromosome length. N-containing
windows discarded (reject rate 6.9%). seed=0.

This isolates "real human DNA distribution" from "open chromatin
annotation" — there's no functional enrichment, just raw genomic
substrate.

## Result
- eval_01 = **0.4992** (vs 001 synth: 0.3068; Δ +0.193)
- eval_07 = 0.5985 (vs 0.4024; Δ +0.196)
- eval_13 = 0.6020 (vs 0.3809; Δ +0.221) — biggest gain
- eval_10 = 0.5529 (vs 0.3566; Δ +0.196)
- **eval_08 = 0.0916** (vs 0.1098; Δ -0.018) — went DOWN

mean over 14 evals ≈ 0.503 (vs 0.308; Δ +0.195).

## What this means
- Real DNA priors are worth ~0.2 r-points across most evals, even with
  no annotation enrichment whatsoever. This is the single largest
  effect I expect to see from any design choice.
- eval_08 *anti-correlates* with real DNA — it gets slightly worse.
  Plausible interpretations: eval_08 is dominated by synthetic /
  designed / motif-implanted sequences whose distribution is far from
  the genomic prior, and the model's genomic prior actively confuses
  it. Or eval_08 has very low total variance and is being driven by
  noise.

## Justification for generalization
Real genomic 200-mers carry k-mer composition, dinucleotide biases,
CpG distribution, repeat content, and motif occurrences that any
mammalian cell-type sequence-to-activity model has to handle. None of
that is K562/HepG2/SKNSH-specific. So training on real DNA exposes
the model to the substrate any future cell-type evaluation will see.
The big jump from synth → real-DNA confirms this prior matters a lot.

## What this tells me for next experiment
- Real DNA matters a lot; DHS annotation might add only a little more
  on top. I want to measure that gap explicitly in E3 (DHS-balanced
  sequences, same overall structure).
- Watch eval_08 — it behaves opposite to the others. Any design that
  improves eval_08 without crashing eval_01/02/05/14 is doing
  something orthogonal that I should study.
