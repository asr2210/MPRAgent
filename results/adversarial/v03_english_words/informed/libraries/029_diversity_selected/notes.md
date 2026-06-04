# 029 — Oversample 150k + greedy 6-mer diversity selection (BREAKTHROUGH)

## Plan
Pivot away from per-cell-type hand-tuned configs. Generate 150k
candidates from a MIX of three best-known designs (020 dense8 random
bg, 022 dense10 random bg, 026 dense8 dinuc bg). Greedily select 50k
maximizing 6-mer novelty (each pick scored by sum of 1/(1+covered[km])
over its 6-mers). Sample-200-per-round greedy for tractability.

## Result
**eval_01 = 0.4288 — FIRST TO BREAK PLATEAU.**
K562 0.594, HepG2 0.627, **SK-N-SH 0.0653 (new mean high)**.
eval_04 SK-N-SH = **0.0743**, eval_07 = 0.0730, eval_09 = 0.0743
(three evals above 0.07 vs zero before exp 026).
14-eval avg = **0.4265** (prior best 008/024 = 0.4253).

Pool composition picked: A (020-style) = 24832, B (022-style) = 17472,
C (026-style) = 7696. Selection naturally prefers random-bg over
dinuc-bg (dinuc seqs have lower 6-mer diversity due to dinuc
correlations), but keeps a non-trivial dinuc fraction for SK-N-SH lift.

## What this teaches
- **Library composition that LOOKS optimal per cell type isn't optimal
  per training signal.** The eval rewards diverse short-context
  patterns. A library can have great motif content but if the
  short-context distribution is redundant, the model can't generalize.
- **Diversity selection is a UNIFYING lever** that captures partial
  per-cell-type benefits without the lever-stacking problem:
  - K562/HepG2 get good motifs + GC content from pool A/B
  - SK-N-SH gets some composition diversity from pool C
  - Greedy 6-mer selection prevents redundant sequences from
    dominating, giving the surrogate richer training distribution
- **The plateau wasn't metric-imposed (T19) — it was DESIGN-imposed.**
  Single-design libraries inevitably oversample particular k-mer
  distributions; mixed-and-selected libraries don't.

## Theory T20 (PLATEAU BREAK MECHANISM)
The surrogate's generalization ability is limited by the EFFECTIVE
DIVERSITY of its training set, not by any per-sequence design choice.
A library of 50k sequences with high motif density but low k-mer
diversity is effectively "smaller" than the same nominal size with
high k-mer diversity. Greedy diversity selection from a candidate
pool of mixed designs produces a more informative training signal.

This explains:
- Why pool-restriction always hurt (less candidate diversity →
  selected pool has less variety either way)
- Why dense overlap helped (more motifs = more short-context variety
  per sequence)
- Why dinuc bg helped SK-N-SH (its oracle reads composition;
  mammalian dinuc adds new k-mer contexts not in random bg)
- Why straight hybrid (024) failed but selected hybrid (029) won
  (greedy selection picks complementary sequences, not redundant ones)

## Next
Exp 030: refine the diversity-selection approach to push further:
- Larger candidate pool (200k from 4 designs)
- Add 4th pool: 008-style (5 motifs no-overlap from celltype-targeted
  289) — adds another distribution shape
- Tune k? (k=5 might capture motif-scale, k=7 might be too sparse)
- Final library + multi-seed variance check
- Write comprehensive summary in notebook
