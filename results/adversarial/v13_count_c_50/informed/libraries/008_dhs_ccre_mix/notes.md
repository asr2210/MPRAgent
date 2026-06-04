# Experiment 008 — dhs_ccre_mix

## Design
25 000 uniform random DHS + 25 000 class-balanced cCRE (3 125 / class × 8).

## Hypothesis
DHS contributes accessibility/cell-type-program coverage (strongest on
eval_07/13). cCRE class balance contributes functional-class coverage
(strongest on eval_04/08/09). Complementary → mix beats either alone.

## Result — **new best**
| eval | 002 DHS | 006 cCRE CB | **008 mix** |
|------|--------:|------------:|------------:|
| 01   | 0.5627 | 0.5637 | **0.5671** |
| 04/09| 0.5405 | 0.5841 | 0.5723 |
| 07   | 0.5974 | 0.5707 | 0.5833 |
| 08   | 0.1663 | 0.2333 | 0.2111 |
| 13   | 0.5771 | 0.5473 | 0.5603 |
| mean | 0.516  | 0.547  | **0.554**  |

eval_01 +0.003 over 006 (small but the best across all four candidates).
The mean improvement is real (+0.007 over 006, +0.038 over uniform DHS).
Per-eval scores mostly land between the two source-alone scores — i.e.
mixing averages each eval's strength weighted by what fraction of the
library carries that strength.

## Interpretation
1. **Sources are complementary, not redundant.** Different evals favor
   different sources, and the mix banks each strength partially.
2. **No source wins on every eval.** eval_07/13 want DHS-style
   accessibility info; eval_04/09 want cCRE-style functional class
   balance. There's no single best regulatory pool.
3. The eval_01 gain is small but the mean gain is meaningful. **Mean
   across 14 evals is increasingly tracking generalization potential
   better than eval_01 alone** — eval_01 has flattened near 0.56 across
   most reasonable libraries, while the mean keeps moving.

## Theory update
- Best library design will likely combine 3+ complementary regulatory
  sources, each contributing a different axis of diversity:
  (i) DHS for cell-type / accessibility-program coverage,
  (ii) cCRE class balance for functional regulatory grammar,
  (iii) ?? — possibly motif-explicit content, ChIP-seq, promoter
       regions, conservation-prioritized regions, or real MPRA seqs.
- A fourth axis that's notably absent: **explicit TF motif content**.
  All my libraries today carry motifs implicitly via genomic context.
  Inserting known motifs into varied backgrounds is the canonical way
  to boost a TF→activity model's ability to decompose.

## Numbers
mean_r averaged across 14 evals: 0.554
eval_01: 0.5671 (new best)
eval_08: 0.2111
time_s: 21.1

## Next
Three directions worth testing in parallel:
- 009: cCRE with PLS+pELS oversampled (active regulatory tilt).
- 010: Add a third source — promoters / 5'UTR / conserved-element
  region pool.
- 011: Augmentation — reverse-complement each element (half × 2).
