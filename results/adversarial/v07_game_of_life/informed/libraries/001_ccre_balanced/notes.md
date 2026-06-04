# Experiment 001 — cCRE category-balanced

## Design
50,000 sequences = 10,000 each from PLS, pELS, dELS, CTCF-only, DNase-H3K4me3
(ENCODE SCREEN cCREs v3). Each is a 200bp window centered on the cCRE midpoint,
extracted from hg38 (skip if any N). Random seed 1.

## Hypothesis
Balanced functional categories of regulatory elements should teach the model
universal regulatory primitives that generalize across cell types.

## Result
- **eval_01 = 0.3919** (mean_r across K562, HepG2, SK-N-SH)
- Mean across all 14 evals: **0.3812**
- Per-cell-type: K562 0.604, HepG2 0.428, SK-N-SH 0.143
- Runtime: 1237s (~20 min)

## What this tells us
**HUGE FINDING.** The cCRE library scored at the random baseline. Comparing to
strategies.md baselines (the file calibrated for this run's evaluator):
- gc_50 (random, 50% GC): 0.3972
- random_uniform: 0.3951
- my cCRE library: **0.3919**

So biologically curated regulatory elements give NO IMPROVEMENT over random
sequences on this evaluator. This contradicts the `instructions.md` baseline
table (which claims dhs_topic = 0.7232) — those numbers appear inflated and not
calibrated to my evaluator. The `strategies.md` (v07) numbers are real for me.

## What I'm changing in my mental model
1. **The instructions.md baseline table is not calibrated to my evaluator.**
   The strategies.md baselines (showing 0.18 - 0.40 for simple strategies) are
   the correct comparison. Top of leaderboard for me is ~0.40, not ~0.72.
2. **Biological elements do not automatically help.** Curation alone (selecting
   "real" regulatory elements) is insufficient. The model must have a reason
   to learn FROM the curated structure.
3. **Cell-type imbalance is severe.** K562 trains 4x better than SK-N-SH.
   Any improvement must lift HepG2 and SK-N-SH meaningfully, since K562 is
   already near the floor performance for any sane sequence set.
4. **GC composition matters.** gc_50 is best among simple baselines. Random
   sequences with controlled GC outperform AT-rich and GC-rich. This tells
   me the eval distribution favors balanced GC content.

## Eval set structure observation
My results show pairs of identical scores:
- eval_01 = eval_02 = eval_05 = eval_14 (0.3919-0.3924)
- eval_03 = eval_12 (0.3856)
- eval_04 = eval_09 (0.3948)
- eval_06 = eval_11 (0.3897)
- eval_08 is consistently hardest (0.2675)

So the 14 eval sets cluster into ~5-7 distinct test distributions, with
the others being duplicates. eval_08 is the outlier (much harder).
