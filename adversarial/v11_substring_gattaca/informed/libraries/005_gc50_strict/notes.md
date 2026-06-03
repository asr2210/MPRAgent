# E5: gc_50_strict (attempted gc_50 replication)

## Design
- Each seq: pick 100 GC positions + 100 AT positions, assign G/C with P=0.5
  and A/T with P=0.5 independently.
- Intent: replicate published gc_50 (5-seed 0.8591 eval_01).

## Result (mean across 14 evals ≈ 0.801)
- eval_01: 0.8197 (vs published gc_50 0.8591, **-0.039** — large miss)
- mean across 14 evals: ~0.801

## Cell type breakdown
- HepG2 high (~0.90 across most evals — gained vs random_uniform)
- K562 ~ unchanged
- SK-N-SH crashed (eval_07 SKNSH = 0.377, like E3)

## Interpretation
My implementation is NOT equivalent to the published gc_50. The
position-typing (pre-allocate which positions are GC vs AT, then sample
within type) reduces local per-position randomness compared to the
published approach. Effectively this gives sequences with reduced
positional variance — similar to E3 in spirit.

The published gc_50 is probably implemented as IID per-position sampling
followed by acceptance/rejection on total GC=100, which preserves
positional randomness while constraining only the global count.
Or — more likely — published "gc_50" simply uses P(C)=P(G)=0.25 IID
with no per-sequence constraint, which IS random_uniform. The 0.003
difference between published gc_50 and random_uniform is likely noise.

## Critical update
The path of "constrain composition more tightly" is closed. Both
E3 (per-base balance) and E5 (positional type constraint) hurt SK-N-SH.
The published gc_50 vs random_uniform difference of 0.003 is within
plausible 5-seed noise.

**Effective baseline: random_uniform at 0.857 eval_01 (0.841 mean across evals).**
Beating this requires a fundamentally different approach — not just
tightening the uniform distribution.

## Pivot
Stop probing within-distribution tweaks. Test whether natural-genome-like
DINUCLEOTIDE structure (esp. CpG depletion) helps despite mild
distribution shift — this is the most natural prior beyond uniform IID
and has been shown to matter in genomics. E6 will test this.
