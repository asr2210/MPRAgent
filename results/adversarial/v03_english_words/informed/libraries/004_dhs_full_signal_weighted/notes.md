# 004 — Full DHS Index, mean_signal weighted

## Plan
Switch to the full 3.6M-element Meuleman DHS Index. Sample with weight ∝
mean_signal, summit-centered 200bp extraction from hg38. Expected this to
recover ≥0.45 eval_01 by drawing from a representative regulatory pool.

## Result
**eval_01 = 0.3868.** Time 1m55. K562/HepG2 ~0.55, SK-N-SH ~0.06.

Essentially identical to exp 001/002 with the synthseqs subset. **Switching to
the full 3.6M-element pool did NOT help.**

## What this means
Three independent DHS-derived libraries (synthseqs-stratified, synthseqs-
numsamples, full-index-signal-weighted) all hit ~0.39 eval_01. Random uniform
hits 0.42. So **natural DHS sequences underperform random in MY pipeline by
~10%**, irrespective of source or weighting.

This is the OPPOSITE of what the instructions.md table shows (dhs_topic 0.72 >
synth_oracle 0.68). My pipeline doesn't reward DHS sequences.

## GC analysis
- Random uniform: 50% ± 3.5% GC (very tight)
- Full DHS (this exp): 46.6% ± 9.8% GC (wider, slightly AT-biased)
- Synthseqs: 44.1% ± 8.0% GC

Baseline ordering:
- gc_50 (50% ± 0): 0.4243
- random_uniform (50% ± 3.5): 0.4228
- gc_sweep (0% to 100%): 0.3334
- at_rich (20% GC): 0.3118

So GC variance hurts in this pipeline. My DHS libraries have wide GC variance.
GC variance is a candidate explanation for DHS underperformance.

## Next
Exp 005: filter DHS to 45-55% GC. If GC variance is the issue, this should
recover ≥0.42. If still ~0.39, the issue is deeper than GC.
