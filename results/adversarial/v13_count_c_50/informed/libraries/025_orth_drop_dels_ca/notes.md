# Experiment 025 — orth-DHS + cCRE-no-dELS-no-CA

## Design
25 k orth-DHS (numsamples ≤ 5) + 25 k cCRE balanced across 6 classes:
CA-CTCF, CA-H3K4me3, CA-TF, PLS, TF, pELS (4166 each).

Dropping CA (DNase-only) on top of dropping dELS (generic distal).
Rationale: DHS half already covers DNase signal; CA is redundant.

## Result — NEW BEST eval_01
| eval | 020    | 023    | **025** | Δ vs 023 |
|------|-------:|-------:|--------:|---------:|
| 01   | 0.5745 | 0.5755 | **0.5762** | +0.0007 |
| 04/09| 0.5509 | 0.5528 | 0.5620 | +0.009 |
| 07   | 0.6168 | 0.6183 | 0.6131 | −0.005 |
| 08   | 0.1404 | 0.1402 | 0.1591 | +0.019 |
| 13   | 0.5968 | 0.5981 | 0.5925 | −0.006 |

### Cell-type lift on eval_01
|         | 023    | **025**    | Δ      |
|---------|-------:|-----------:|-------:|
| K562    | 0.6119 | 0.6155 | +0.004 |
| HepG2   | 0.5532 | 0.5529 | −0.000 |
| SK-N-SH | 0.5615 | 0.5602 | −0.001 |

## Interpretation
Dropping CA recovered K562 score (+0.004) while keeping HepG2/SK-N-SH
stable. Net eval_01 +0.0007.

**Why K562 recovers**: Removing 4166 CA elements (DNase-only) and
giving that budget to the 6 functionally-marked classes means more
high-signal regulatory elements per class. The model trains on a
denser regulatory grammar signal, which K562 (the cleanest cell line
in ENCODE) benefits from most.

eval_08 also lifts (+0.019) because removing CA freed budget for
PLS+pELS (both increased to 4166 from 3571 in 023, +17 %), adding
more CpG-rich content.

eval_07/13 slightly drop because the dropped classes (CA, dELS) had
diverse content — losing them tightens the regulatory grammar pool
but reduces breadth.

## Theory update
- **Aggressive cCRE class pruning continues to help eval_01**. Three
  rebalances now confirmed:
  - Drop dELS (023): +0.0010 vs 020
  - Drop dELS + CA (025): +0.0017 vs 020
- **The 6 retained classes (CA-CTCF, CA-H3K4me3, CA-TF, PLS, TF, pELS)
  are the "functionally specific" cCRE subset** — each requires
  evidence beyond DNase alone. They carry more regulatory information
  per element.

## Numbers
mean_r: 0.555
eval_01: 0.5762 (new best, +0.0007 over 023)
eval_04: 0.5620
eval_08: 0.1591
time_s: 26

## Next
Test if even more aggressive pruning helps:
- 026: drop CA-CTCF too (CTCF = insulator, not enhancer activity)
- Or: differently weight the 6 classes by some quality metric.
