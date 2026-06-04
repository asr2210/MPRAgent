# Experiment 010 — dhs_ccre_mix_rc_aug

## Design
25 000 unique elements (12.5 K DHS uniform + 12.5 K cCRE class-balanced)
+ reverse-complement of each = 50 000 sequences total.

## Hypothesis
Explicit RC augmentation might help the model learn strand-symmetric
regulatory grammar. If yes, beats 008 (0.5671). If no, halving unique
elements was a wash and the model learns strand symmetry implicitly.

## Result
| eval | 008 (50 k unique) | **010 (25 k × 2 strands)** |
|------|------------------:|---------------------------:|
| 01   | 0.5671 | 0.5659 (−0.001) |
| 07   | 0.5833 | 0.5852 |
| 13   | 0.5603 | 0.5621 |
| 04/09| 0.5723 | 0.5674 |
| 08   | 0.2111 | 0.2083 |
| mean | 0.554  | 0.553 |

**Statistical tie.** RC augmentation neither helps nor hurts.

## Interpretation
Two readings, both interesting:
1. **Halving element diversity didn't hurt eval_01.** At 25k unique
   elements with RC-paired duplicates, the model performs essentially as
   well as 50k unique elements. Element diversity is saturating —
   somewhere between 25k and 50k, additional elements stop adding new
   information.
2. **RC augmentation didn't help either.** The model probably learns
   strand symmetry implicitly from the genomic distribution; an
   architectural prior would help but adding RC examples doesn't.

## Theory update
- **Element diversity saturates around 25k.** This has a counterintuitive
  consequence: spending the 50k budget on *fewer unique, higher-quality
  elements* with augmentation might be as good as 50k unique. Could be
  useful if I find a small high-quality source.
- **Strand augmentation is a wash.** Skip this from now on.

## Numbers
mean_r: 0.553
eval_01: 0.5659
eval_08: 0.2083
time_s: 54

## Next
Now that I know element diversity saturates and RC augmentation doesn't
help, the path to higher eval_01 is probably either (a) finding a new
data source with higher per-element info density on the eval distribution,
or (b) cell-type-balanced design that addresses the K562-bias I see in
all results so far (K562 score is always 0.05–0.10 above HepG2/SK-N-SH).
