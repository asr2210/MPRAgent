# 018_repeat_rich_only

## Design
50k random tiles from chr1/17/19/22, ONLY those FAILING the E17
filter (i.e., repeat-rich / low-complexity by any of:
max single-nuc > 0.45, 2-mer entropy < 3.0, homopolymer ≥ 5).
Acceptance: 70% (50k of 71k attempts).

## Result — informative tie
- eval_01 = **0.4937** (vs E2 0.4992; Δ -0.006) — at noise floor
- mean over 14 ≈ 0.490

Repeat-rich-only is effectively TIED with full random.

## Interpretation — major insight refinement
The 70% repeat-rich half of random carries virtually ALL the signal.
The 30% high-complexity half (E17 = 0.407) is much less informative
per slot.

Per-slot informativeness:
- Repeat-rich tile (70% of random): contributes most of E2's score
- High-complexity tile (30% of random): contributes much less

This INVERTS the typical "high complexity = more information"
intuition for this pipeline. The model learns the human prior
DOMINANTLY from repeat-rich regions (LINEs, SINEs, ALUs, simple
repeats, mononucleotide tracts), with high-complexity coding-like
regions adding marginal value.

## Theory v8

The "+0.19 from real DNA" prior breaks down as:
- ~0.17 from the model seeing realistic REPEAT GRAMMAR (mostly ALU,
  LINE, microsatellites, mononucleotide stretches).
- ~0.02 from high-complexity gene-body / coding-like regions.

This is the OPPOSITE of what I'd guess from "regulatory elements
matter most". The eval set is dominated by something that REPEAT
content captures and curated regulatory subsets miss.

Possible explanation: many MPRA libraries have substantial repeat
content because random fragments often hit them. If the eval set
itself is random/DHS-style sampling, repeat content is well-
represented in eval. Training on repeat-rich tiles matches the
eval distribution best.

## Implications

- **Stop curating away from natural sampling.** Any filter that
  shifts the complexity profile (high OR low) hurts.
- **Repeat-rich grammar is the primary learnable signal.** Tools
  that explicitly model repeat structure might help (out of scope).
- **The 0.50 ceiling is robust because no design EXCEEDS the natural
  complexity profile in informativeness.** I cannot beat broad random.

## Plan for E19

Try **gene-density weighted full-genome** sampling: weight 1Mb bins
by gene count, sample length-weighted across ALL chromosomes
(including chrX, chrY, chr2-22). This is the most "natural-flavored"
broadening of E2 (chr1/17/19/22).

Predicted: ~0.50. Probably saturated at the natural-prior ceiling.
If <0.49: gene-density continuous weighting hurts. If >0.50:
continuous gene weighting beats hard chromosome filter (small upside).
