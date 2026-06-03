# Experiment 021 — Mixed Gosai lfcSE<0.7 random

## Method
Full Gosai (CRE+GTEX+UKBB) with lfcSE<0.7 (726K). Random 50K.

## Result
**eval_01 = 0.0137** — worse than both reference points:
- Mixed Gosai lfcSE<0.5 random (014): 0.0168
- GTEX-only lfcSE<0.7 random (019):    0.0222

## Interpretation
Loose quality + mixed sub-source is BAD. The UKBB sequences with looser quality
contribute noisy labels that hurt training. For UKBB, only the cleanest
sequences are useful; for GTEX, even noisier ones contribute.

Why? GTEX (eQTL) sequences were pre-screened for functional effect on
expression. Even at higher SE, they have non-trivial true effects. UKBB
(GWAS) sequences are mostly null variants — at high SE, the labels are
basically noise.

## Theory consolidated — GTEX advantage compounds
1. GTEX has higher base rate of TRUE non-zero effects (eQTL screening)
2. So GTEX tolerates more measurement noise per sequence
3. Looser quality on GTEX → more sequence diversity at acceptable label noise
4. Looser quality on UKBB → adds null variants with noise = pure dilution

The BEST library = GTEX-only with the loosest quality filter that still
excludes pure noise.

## Next direction
Fine-grained search around lfcSE=0.7 (try 0.6, 0.8) and seed variance check.
