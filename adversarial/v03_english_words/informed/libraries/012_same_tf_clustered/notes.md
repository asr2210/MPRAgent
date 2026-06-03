# 012 — Same-TF clustered (3 instances of one TF per seq)

## Plan
Same 289-motif pool as exp 008 (best, 0.4283), 3 motifs per 200bp random
backbone, BUT all 3 motifs are stochastic samples of the SAME TF (picked
randomly per sequence). Tests whether per-TF density amplifies signal vs
3 different TFs.

## Result
**eval_01 = 0.4196.** Worse than exp 008 (0.4283) by -0.009.
K562: 0.590, HepG2: 0.620, SK-N-SH: 0.049.

## What this teaches
- Per-TF clustering doesn't help — within-sequence TF diversity beats
  TF repetition.
- The surrogate doesn't use motif "stamps" or per-TF density as a primary
  feature. It cares more about WHICH TFs are present (set membership) than
  HOW MANY copies of each.
- Combined with exp 010 (cell-type-clustered, 0.4217 < 0.4283) and exp 011
  (consensus, 0.4169 < 0.4283), three different clustering/specialization
  strategies have lost to the mixed/diverse strategy.

## Theory T10
The surrogate is essentially detecting a "TF presence" feature set.
Sequences are most informative when they probe MANY distinct TF presence
features simultaneously. Per-sequence diversity (different TFs in each seq)
covers more of the TF feature space per training sample than clustered
designs. Concretely: a sequence with {GATA1, FOXA1, NEUROD1} teaches the
model about all three TFs; a sequence with {GATA1, GATA1, GATA1} just
reinforces one feature it already saw from neighbors.

Generalization implication: a library designed for unknown cell types
should expose the model to as many distinct TF binding patterns as possible
per sequence — not double down on any one TF, even one that matters for the
labeling cell types.

## Next
Exp 013: focused 50-motif pool (drop universals and rare TFs, keep only
the 50 most-likely-active cell-type TFs). Tests whether 289 is too broad
and dilution from rarely-impactful PFMs hurts. If exp 013 > 0.4283: pool
quality > pool size. If <: 289 was already near-optimal.
