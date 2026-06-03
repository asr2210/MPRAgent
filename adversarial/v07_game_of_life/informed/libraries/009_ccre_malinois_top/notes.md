# Experiment 009 — cCRE + Malinois top-magnitude

## Design
All 1.06M ENCODE cCREs → 200bp center windows (with jitter) → filter to GC
[0.45, 0.55] → 1.9M candidates → score with Malinois → top 50k by max log2FC.

Selection stats:
- Pool max-cell score: range [-1.07, 12.62], median 0.66
- Selected: min max-cell score 3.84, max 12.62
- Selected pred means: K562 4.10, HepG2 3.43, SK-N-SH 3.72
  (vs 007 random+oracle means: K562 2.60, HepG2 2.36, SK-N-SH 2.63)
- GC: mean 0.501, std 0.029

## Hypothesis
Combining real enhancer sequence statistics with oracle activity selection
should outperform either alone (random+oracle = 0.397, cCRE alone = 0.392).

## Result
- eval_01 = **0.3879** (vs 007 = 0.3969, surprisingly WORSE)
- Mean across 14 evals = **0.3774** (vs 007 = 0.3861)
- Per-cell-type on eval_01: K562 0.598 (−0.020), HepG2 0.422 (−0.014),
  SK-N-SH 0.144 (+0.007)
- Runtime: 1201s

## Interpretation
**Surprising negative result.** Despite 60% higher predicted activity in the
selected pool, cCRE+oracle UNDERPERFORMS random+oracle on K562 and HepG2.

Likely mechanism: **Malinois OVERPREDICTS on cCRE-like sequences** (because
its training data was enriched in cCRE-derived sequences). Selecting "top
Malinois activity" from cCREs picks sequences that look very enhancer-like
to Malinois but whose REAL MPRA activity is normal. Training the small
model on these biased-high training labels then overshoots on similar
test sequences and hurts correlation.

Equivalent statement: when the oracle and target distributions differ,
selecting the oracle's "best" can bias training away from the target
distribution. Random GC-50 candidates are unbiased; Malinois may be more
calibrated on them than on cCREs.

The lone bright spot: SK-N-SH lifted slightly (+0.007). This might mean
neural-lineage TF motifs in cCREs are providing real signal — worth
following up.

## What this rules out
- "Biology + oracle" is not strictly additive
- Higher Malinois prediction does NOT translate to higher eval score
- The 007 ceiling (~0.397) holds against this attempt

## Theory v6 update
The Malinois oracle is calibrated on a random-like background. Using it to
select from a biased pool (cCREs) creates a distribution mismatch that
hurts. The best Malinois-based library is still 007 (random+top), capped
at 0.397.

To break the ceiling we likely need either:
1. A different oracle (or ensemble) less biased
2. A different lever entirely (motif syntax targeted at SK-N-SH)
3. Hyperactive generation in random sequence space (sequences not
   accessible by random sampling)
