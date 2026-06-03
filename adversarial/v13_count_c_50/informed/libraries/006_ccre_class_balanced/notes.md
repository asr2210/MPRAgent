# Experiment 006 — ccre_class_balanced

## Design
Sample 6 250 from each of the 8 SCREEN cCRE classes (PLS, pELS, dELS,
CA, CA-CTCF, CA-H3K4me3, CA-TF, TF). 200 bp centered on midpoint.

## Hypothesis
cCRE classes encode functionally distinct regulatory grammars (promoter
vs distal enhancer vs CTCF-only vs TF-only), unlike DHS NMF topics
which are co-accessibility patterns. Forcing balance across these
functional categories should improve generalization vs uniform DHS
(0.5627).

## Result — **new best**
| eval | uniform DHS (002) | cCRE class-balanced (006) | Δ |
|------|------------------:|--------------------------:|------:|
| 01   | 0.5627 | **0.5637** | +0.001 |
| 04/09| 0.5405 | **0.5841** | +0.044 |
| 07   | 0.5974 | 0.5707 | −0.027 |
| 08   | 0.1663 | 0.2333 | +0.067 |
| 13   | 0.5771 | 0.5473 | −0.030 |
| mean | 0.516  | **0.547**  | **+0.031** |

eval_01 barely moves but the mean across all 14 evals jumps by 0.031.
Big winners are 04/09 (+0.044) and 08 (+0.067 — without adding any
synthetic). Losers are 07 and 13.

## Interpretation
1. **cCRE > DHS** as a substrate, at least at 50k. Likely because the
   chromatin-mark filter selects elements with stronger evidence of
   actual regulatory function.
2. **Class balance is doing real work.** In a uniform cCRE pool dELS
   would be 62 % of draws and PLS only 2 %. Balancing forces ~6 % each.
   PLS (active promoters) carries dense, learnable TF-grammar that
   the model can extract.
3. **eval_08 boost without synthetic** is interesting. cCRE midpoint
   windows have slightly different base composition than DHS summit
   windows — could be that. Or PLS / CA classes happen to look more
   like the eval_08 distribution.
4. eval_07 and eval_13 dropped — those favored DHS programs (per
   strategies.md analysis). The cCRE library doesn't contain all the
   information DHS does. Suggests a mix should win on both fronts.

## Theory update
- Functional-class balance > co-accessibility-topic balance. Diversity
  matters but the axis matters too — divide the space by *what the
  element does*, not by *which cells it's open in*.
- Mixed-source libraries (DHS + cCRE) are now the most promising
  direction. DHS contributes accessibility/cell-type diversity; cCRE
  contributes regulatory-class diversity. They're complementary.

## Next
Two diagnostic experiments to disentangle the source-vs-class-balance
effects:
- **007** cCRE uniform random (no class balance). If ≈ 002 = source
  isn't the win, class balance is. If ≈ 006 = source is the win,
  balance is gravy.
- **008** 50/50 DHS + cCRE mix, both uniform. Tests whether sources are
  complementary — should beat both 002 (0.5627) and 006 (0.5637).

## Numbers
mean_r averaged across 14 evals: 0.547
eval_01: 0.5637 (new best)
eval_08: 0.2333
time_s: 22.8
