# E7: true_gc_50 — IID + rebalance to exact GC=100

## Design
IID per-position sampling at P=0.25 each, then for each seq flip random
excess GC/AT positions until exactly 100 GC. Preserves per-position
IID more than E5 (position-typing) yet enforces exact GC=100 per seq.

## Result (mean across 14 evals ≈ 0.803)
- eval_01: 0.8218 (vs E2 0.8565, **-0.035**)
- SK-N-SH eval_07: 0.3757 (vs E2 0.7200, -0.34 — crashed)

## Conclusion
Even "true" gc_50 with maximally-preserved IID-per-position structure
LOSES big to random_uniform. The published gc_50 baseline at 0.859 is
either (a) inconsistent with my pipeline, (b) effectively equivalent to
random_uniform (the 0.003 difference being noise), or (c) implemented
differently than any interpretation I've tried.

**Most likely**: the published 5-seed baselines have ~0.005 SD per seed
and the gc_50 vs random_uniform difference is within noise. My
single-seed runs reproduce random_uniform (E2 = 0.857 vs published 0.857)
but show GC=100 constraint loses 0.034 on eval_01.

## Lock-in
random_uniform (IID 50% GC, no per-seq constraints) is the eval distribution
within available precision. The ceiling at ~0.857 eval_01 cannot be moved
by distributional tweaks. SK-N-SH is the canary for any deviation:
constraint → SK-N-SH crash → mean drop.

## Path forward
Move beyond distributional tweaks. Test:
- Active sampling from larger pool (E8 candidate: antithetic / RC pairs)
- Information-theoretic selection (k-mer entropy, but minimal expected gain)
- Run-to-run seed variance to calibrate noise floor
- Or: accept ~0.857 as near-ceiling and characterize WHICH random-like
  libraries reproduce it best
