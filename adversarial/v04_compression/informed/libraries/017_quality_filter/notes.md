# 017_quality_filter

## Design
50k random tiles from chr1/17/19/22, but ONLY tiles passing all filters:
- max single-nucleotide frequency ≤ 0.45
- 2-mer Shannon entropy ≥ 3.0 bits
- max homopolymer run ≤ 4

Filter rejected 70% of random tiles (151k of 216k attempts). seed=0.

## Result — STRONG NEGATIVE
- eval_01 = **0.4069** (vs E2 0.4992; Δ -0.092)
- mean over 14 ≈ 0.421
- All evals UNIFORMLY DOWN by 0.05-0.10.

This is roughly the same magnitude as cCRE-balanced (E16 = 0.402).

## Interpretation — BIG INSIGHT
**Repeat-rich / low-complexity DNA is NOT wasted slots — it's
INFORMATIVE.** Removing these tiles cripples the library.

The model needs exposure to:
- Simple repeats and homopolymers
- Tandem repeats and microsatellites
- Mononucleotide tracts

These regions are common in real DNA and presumably common in the
eval distribution. Filtering them shifts the training distribution
away from the eval distribution.

This OVERTURNS my prior intuition that "median random has wasted
slots". The opposite is true: the natural complexity range, including
the seemingly-low-info repetitive content, is what the model needs.

## Theory update v7

**The real-DNA prior IS the broad complexity distribution.** Any
deviation that shifts the complexity profile (filtering up OR down)
hurts. Two strong confirmations:
- E17 (high-complexity only): -0.09
- E5 (Markov k=4, smoothed composition): -0.04 vs uniform synth

The real-DNA win (+0.19 over synth) is NOT about choosing "good"
sequences — it's about exposure to the FULL complexity landscape:
- Repeats (LINE, SINE, ALU)
- Tandem repeats
- Mononucleotide / dinucleotide tracts
- Plus all the "regular" gene-body / intergenic content

Filtering ANY of these out costs ~0.09. Match the natural complexity
distribution or lose.

## What this means for design

1. **Stop filtering**: don't reject "low quality" tiles, ever.
2. **Match natural complexity**: any design choice that shifts the
   complexity distribution will hurt.
3. **The 0.50 ceiling is the natural information content of 50k
   200bp tiles from broad gene-rich human DNA.** It's not a soft
   ceiling — it's where the training distribution best matches
   whatever the eval distribution is.

## Plan for E18
Test the converse: **inverse filter** — keep ONLY tiles that FAIL
the E17 filter (repeat-rich, low-complexity).

Predicted:
- If repeat-rich alone is too narrow: ~0.40 (similar to E9 / E16 /
  any narrow-content design).
- If repeat-rich is INFORMATIVE in its own right but lacking
  broad-complexity diversity: ~0.45.
- If repeat-rich is the dominant useful signal: >0.50 (would be a
  huge finding).

Likely outcome: ~0.42-0.45. Repeat-rich is informative WITH the rest
but not on its own.
