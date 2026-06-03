# Experiment 004 — cCRE GC-matched (composition control)

## Design
50k 200bp windows from ENCODE cCREs, but filtered so only sequences with
GC content in [0.45, 0.55] are kept. Mean GC 0.496, std 0.030 — tightly
matched to the gc_50 baseline composition.

## Hypothesis
If composition is the dominant lever (from exp 003), composition-matched
biological cCREs should score ≈ gc_50 (0.397). If biology adds value when
composition is controlled, > 0.40.

## Result
- eval_01 = **0.3921** (vs cCRE balanced 0.392, gc_50 0.397, PLS 0.375)
- Mean across 14 evals: **0.3811**
- Per-cell-type: K562 0.606, HepG2 0.425, SK-N-SH 0.145
- Runtime: 1209s

## Interpretation
Composition matching DID lift PLS-style cCRE selection back to the cCRE
baseline level, but did NOT exceed gc_50 random. **Biological cCRE content
provides zero measurable lift over random sequences at matched composition.**

Composition + biology comparison:
| library         | mean GC | eval_01 |
|-----------------|---------|---------|
| gc_50 (baseline)| 0.50    | 0.397   |
| 001 cCRE bal    | 0.51    | 0.392   |
| 003 PLS only    | 0.62    | 0.375   |
| 004 cCRE-gc-mat | 0.50    | 0.392   |
| 002 motif+gc50  | 0.50    | 0.394   |

ALL composition-matched libraries fall in 0.390-0.397, regardless of content.
The composition explanation is the FULL explanation for these data.

## What this rules out
- Biological element CURATION (cCREs) doesn't help
- Motif IMPLANTATION doesn't help
- Specific element types (PLS only) don't help (and hurt via GC bias)

## What's left to test
- Is GC=0.50 truly optimal, or is there a slightly different sweet spot?
- Does the SOURCE of GC=0.50 sequences matter? (Random natural genomic
  windows vs random i.i.d. vs cCREs vs synthetic)
- Does cell-type-specific motif content help raise SK-N-SH score (lagging)?
- Does combining multiple GC-50 strategies (hybrid) lift above 0.40?

## Key uncertainty
I'm at ~0.39 on 4 different strategies, the gc_50 baseline is 0.397, and
strategies.md max is 0.397. Suggests a hard ceiling around 0.40 for ANY
50k-sequence library on this evaluator. May not be possible to break through
substantially.
