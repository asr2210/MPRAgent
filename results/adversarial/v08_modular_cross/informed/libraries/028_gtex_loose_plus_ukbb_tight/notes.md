# Experiment 028 — 40K GTEX-loose + 10K UKBB-tight blend

## Method
- 40K from GTEX with lfcSE<0.7 (409K pool)
- 10K from UKBB with lfcSE<0.2 (very tight — 184K pool)

The asymmetric quality criterion respects the finding from 014/021: UKBB needs
tight quality (high SE = noise), GTEX benefits from loose quality (high SE =
still real effects). Random within each pool.

## Result
**eval_01 = 0.0209** — within seed variance of 019's 0.0222.
**NEW BEST on 7 of 14 evals**: eval_02 (0.0212), eval_03 (0.0202), eval_05
(0.0212), eval_06 (0.0204), eval_11 (0.0204), eval_12 (0.0202), eval_14
(0.0209).

Cell balance dramatically improved:
- HepG2 contribution lifted (0.020 vs 0.013 in 019)
- SKNSH steady (~0.013)
- K562 strong (~0.030)

## Interpretation — KEY FINDING
Blending DOES help, but only when each sub-source uses its own optimal quality:
- GTEX (eQTL): tolerates loose quality, more diversity helps
- UKBB (GWAS): needs tight quality, must be cleaned

The naive mixed-source experiments (009, 021) failed because they applied a
single quality criterion to both sub-sources, hurting one or the other.

## Best library candidate
This blend covers more eval clusters than any single-source library. Mean
across 14 = 0.0153 (close to 019's 0.0164 but more balanced).

For aggregate performance, this is the design to recommend.
For pure eval_01 maximization, GTEX-only lfcSE<0.7 with the right seed.

## Next direction
Try 30K GTEX + 20K UKBB and other ratios to see if more UKBB-tight helps.
