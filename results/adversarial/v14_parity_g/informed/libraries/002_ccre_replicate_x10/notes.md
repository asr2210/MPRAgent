# 002_ccre_replicate_x10

## Design
5,000 unique cCREs sampled randomly, each replicated 10× → 50,000 sequences.
All else identical to experiment 001 (same chromosome filter, same 200bp
centering, same seed for the unique pick).

## Hypothesis
Under v14's apparently brief training, replication should increase examples
per learnable pattern and average label noise within each unique seq. If
the bottleneck is examples-per-feature, this would push r above the
exp 001 baseline.

## Result
eval_01 mean_r = 0.0021. Indistinguishable from exp 001 (0.0016). Replication
does not help.

## Interpretation
Either (a) v14's evaluator is producing near-noise on any library that
isn't structured for a very specific purpose, or (b) library replication
isn't the bottleneck. Moving on to test motif-density (003) and bimodal
contrast next.
