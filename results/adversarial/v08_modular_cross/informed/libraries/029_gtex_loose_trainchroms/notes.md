# Experiment 029 — 30K GTEX-loose + 20K UKBB-tight (different blend ratio)

## Method
Same as 028 but shifted ratio: 30K GTEX (lfcSE<0.7) + 20K UKBB (lfcSE<0.2).
Tests whether MORE UKBB-tight content improves over 028's 40+10.

## Result
**eval_01 = 0.0188.** Below 028's 0.0209. Mean across 14: 0.0136 (vs 028's
0.0153).

## Interpretation
More UKBB-tight HURTS. The 40+10 ratio is closer to optimal — UKBB needs to be
a minority addition. Likely because:
- GTEX is the primary signal (driving eval_04/09 + eval_01)
- UKBB-tight only adds value for HepG2/SKNSH balance (eval_06/11 etc.)
- Too much UKBB dilutes the GTEX-driven eval_01 signal

## Final synthesis
Best library configurations identified:
1. **GTEX-only loose-quality**: eval_01 peak ~0.022 (seed-lucky), mean ~0.016
2. **40K GTEX-loose + 10K UKBB-tight blend**: eval_01 ~0.021, mean ~0.015,
   BEST on 7 of 14 evals (most balanced)

For final library, try one variant with smaller UKBB share (e.g., 45+5) to see
if even less UKBB but still some adds value.
