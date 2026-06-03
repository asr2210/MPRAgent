# Experiment 030 — FINAL: 45K GTEX-loose + 5K UKBB-strict

## Method
- 45K GTEX (data_project=="GTEX"), lfcSE<0.7 (409K pool)
- 5K UKBB (data_project=="UKBB"), lfcSE<0.15 (119K pool)
- Random within each pool, shuffle together

Synthesis of 29 prior experiments. Each sub-source uses its OWN optimal quality:
GTEX tolerates loose (more diversity wins), UKBB needs strict (cleanliness wins).
80/20 split came from 028 vs 029 ratio comparison.

## Result
**eval_01 = 0.0212** — best stable performance, top 2 result.
**eval_02 = 0.0216** — NEW BEST.
**eval_05 = 0.0216** — NEW BEST (tied with eval_02 by construction).
**eval_04/09 = 0.0224** — matches 019's seed-lucky best.
**eval_06/11 = 0.0207, eval_14 = 0.0212** — close to best.
**Mean across 14 evals = 0.0158** — highest aggregate of all 30 libraries.

## Best library identified

The recipe that breaks the original 0.018 plateau and achieves balanced
performance across most eval clusters:

```
GTEX (eQTL variants):  45K with mean(lfcSE) < 0.7   [allow noise, take many]
UKBB (GWAS variants):  5K with mean(lfcSE) < 0.15   [strict cleanliness]
Random sampling within each.
```

## Why this works (theory)
1. **GTEX variants are pre-screened for transcriptional effect.** Even at high
   measurement noise, their TRUE effects are non-trivial, so the noisy labels
   are still informative.
2. **UKBB GWAS variants are largely null on transcription.** Their high-SE
   measurements are pure noise, so only the cleanest cases contribute.
3. **80/20 GTEX/UKBB split** keeps the primary GTEX signal dominant while
   adding the complementary information UKBB provides on HepG2/SKNSH (which
   GTEX-only undersamples).

## What did NOT work
- Natural genomic sequences (cCRE, DHS): 0.0 — outside oracle distribution
- Random Gosai (no quality filter): 0.0 — labels too noisy
- Activity extremes / variance / SNR selection: hurts eval_01
- Tight quality (lfcSE<0.3): biases toward easy sequences, hurts diversity
- Per-cell or per-source stratification: distorts natural distribution
- Mixing extremes + middle: sub-additive (each dilutes the other)
- Single-source mix with shared quality threshold (021): worst of both worlds
- Augmentation via reverse complement: wastes slots on duplicates
- Top-K by |effect| (deterministic extremes): boosts eval_04/09 but tanks eval_01

## The three knobs that matter
1. **Sub-source** = composition of training data: GTEX > UKBB ≫ CRE
2. **Quality** = per-source U-curve: too loose ⇒ noise, too tight ⇒ bias.
   Peak depends on sub-source.
3. **Sampling shape** = natural random > any deterministic selection
