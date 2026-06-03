# Experiment 012 — ccre_heavy_mix

## Design
15 k DHS uniform + 35 k cCRE class-balanced (4375/class × 8).
Tests cCRE-heavy ratio vs 008's 25/25 split.

## Result
| eval | 008 (25/25) | 011 (CT-tgt) | **012 (15/35)** |
|------|------------:|-------------:|----------------:|
| 01   | 0.5671 | 0.5688 | 0.5678 |
| 04/09| 0.5723 | 0.5579 | **0.5772** |
| 07   | 0.5833 | 0.5989 | 0.5828 |
| 08   | 0.2111 | 0.1813 | 0.2169 |
| 13   | 0.5603 | 0.5769 | 0.5589 |
| mean | 0.554  | 0.555  | 0.554  |

eval_01 = 0.5678 — within 0.002 of 008. cCRE-heavy doesn't help eval_01.

## Interpretation
**eval_01 plateau confirmed at 0.567–0.569.** Three different mix
designs (008, 011, 012) all land in this band — diminishing returns
from rebalancing DHS:cCRE within an element-level sampling regime.

**eval_04/09 (paired) jumps to 0.5772 — best so far.** Looking at cell
type splits: K562=0.520, HepG2=0.584, SK-N-SH=0.628. eval_04/09 is
clearly NOT K562-biased like other evals (where K562 leads by 0.05+).
It's a non-K562 eval. cCRE class-balancing — which includes lots of
dELS (cell-type-specific) and rare classes — helps non-K562 cells.

**eval_07 loss from 011's 0.5989 → 012's 0.5828**: 011 was specifically
component-targeted, which helped eval_07. 012 went back to uniform
DHS, which gave up that gain.

## Theory update
- **eval_01 needs a different lever**, not ratio tuning. Three designs
  (008/011/012) clustered at 0.567–0.569 → ratio is not the bottleneck.
- **Different evals reward different priors**:
  - eval_07/eval_13: cell-type-component-balanced DHS (011 wins)
  - eval_04/09: cCRE-rich (012 wins)
  - eval_01: ??? (no design dominates)
- Could try a *mosaic* library: pick the design that wins each eval
  and blend in proportion, but that risks "fairy-mean" with no peak.

## Numbers
mean_r: 0.554
eval_01: 0.5678
eval_04: 0.5772
eval_07: 0.5828
eval_08: 0.2169
time_s: 49

## Next
Element-level sampling has plateaued. Try a *qualitatively different*
augmentation strategy — exp 013: take 25 k unique elements (008 design
halved) and add 25 k mutated copies (5 random SNPs each). Tests local
sequence-space densification, which RC augmentation (exp 010) couldn't
test because RC preserves motifs exactly.
