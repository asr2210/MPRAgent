# 009_ctcf_heavy — CTCF doubled, K562 went *negative*

**Composition:** 8k cCRE uniform + 8k CTCF + 4k DNH3 + 5k S2 + 3k DHS + 2k synth + 20k flanks

**Result:** mean_r(14) = **0.1354**, eval_01 = 0.1367.

**Vs exp 007:** -0.018 on mean14, -0.022 on eval_01. SIGNIFICANT REGRESSION.

**The shocker:** K562 r went *negative* on 12/14 evals (range -0.04 to +0.02). Compare exp 007: K562 r was +0.07 to +0.09 on eval_06/11.

| eval | exp 7 (4k CTCF) | exp 9 (8k CTCF) | Δ |
|---|---|---|---|
| eval_06/11 (UKBB/GTEx K562) | 0.228 | 0.144 | **-0.084** |
| eval_01 (chr-GT) | 0.159 | 0.137 | -0.022 |
| eval_07 (SEI) | 0.131 | 0.154 | +0.023 |
| eval_10 (DHS) | 0.121 | 0.158 | +0.037 |
| eval_13 (genomic) | 0.127 | 0.157 | +0.030 |

**Diagnosis:**
CTCF-only cCREs are *insulator* sites with low intrinsic enhancer/promoter activity. They have constrained, repetitive grammar (clustered CTCF motifs, often at TAD boundaries). Doubling them appears to:
1. **Bias the K562 head toward predicting "low/no activity"** — explaining negative correlation on K562-heavy evals (06/11)
2. **Crowd out enhancer-grammar examples** (uniform cCRE dropped from 10k→8k)
3. Marginally help genomic-distributed evals (07/10/13) because CTCF sites ARE genomic, but at huge cost elsewhere

**Theory v6 was wrong.** The composition-score relationship is NOT linear-additive per source. Sources have *interaction effects* — CTCF helps when balanced with enhancers; CTCF hurts when it dominates.

**Theory v7 — Interaction matters, not just counts:**
- Each source contributes *signal* (correlation lift) AND *noise* (anti-correlation drag from out-of-distribution examples for that head).
- The signal/noise ratio is non-linear in the source count: small additions are pure signal, but at high counts the source biases the model's output distribution.
- For CTCF: ~4k is the sweet spot. Going to 8k flips the K562 head's output mean toward "inactive."

**Implication:** Don't simply scale up the best-performing component. Instead, test *symmetric expansion* (more of everything) by scaling flanks down. Or test whether the K562 boost actually came from cCRE uniform + DHS (high-activity enhancers) NOT from CTCF.
