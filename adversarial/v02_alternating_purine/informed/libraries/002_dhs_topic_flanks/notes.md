# Experiment 002 — DHS topic + paired genomic flanks

## Design
- 25k DHS-topic positives (NMF loadings sum, single seed)
- 25k paired flanks: per positive, offset ±[1500, 3000]bp from summit,
  200bp window, reject N or DHS-overlap
- Tests whether "flanks help" (blind's finding) transfers across pools

## Result — eval_01 = 0.1327 (vs 0.1293 in exp 001)
- K562_r = -0.024 (was -0.005 in 001) — DROPPED
- SKNSH_r = 0.446 (was 0.399) — IMPROVED +0.047
- Net change on mean_r is +0.003 (essentially noise)

## Per-eval differences vs 001
- Better: eval_07 (+0.039), eval_10 (+0.017), eval_13 (+0.036)
- Worse: eval_06 (-0.019), eval_11 (-0.019), eval_04 (+0.015)
- Eval_08 unchanged at 0.047 (universal floor)

## Comparison
- My DHS-topic + flanks = 0.133
- Blind's cCRE + flanks (010) = 0.166
- Gap = 0.033 → cCREs ARE meaningfully better positives than DHS-topic
  under the same flank recipe. Source pool matters more than I thought.

## Theory update
**Theory v2**: The cCRE pool has structural advantage over DHS for this
harness. Possible causes:
- cCREs encode discrete functional classes (PLS, pELS, dELS, CTCF, DNH3) —
  these may map cleanly to features the model exploits
- cCREs are curated / smaller (~1M vs 3.6M); uniform sampling captures more
  of the relevant diversity
- DHS topic-weighting may over-represent tissue-specific elements at the
  expense of canonical regulatory elements

Adding flanks helps SKNSH (the only working cell channel) but slightly hurts
K562. Net win comes only when positives are themselves "good."

## Next
Exp 003: Use cCRE as positives + paired flanks to confirm I can reproduce
blind's 0.166 baseline. Then build variations from there.
