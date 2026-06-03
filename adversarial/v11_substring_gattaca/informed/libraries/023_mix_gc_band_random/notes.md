# E23: mix_gc_band_random

## Design
25k IID random with GC ∈ [90, 110] + 25k unfiltered random IID.

## Result
- eval_01: 0.8672 (vs E15 0.8777, -0.011)
- mean_r: 0.8512 (vs E15 0.8591, -0.008)
- SK-N-SH eval_07: 0.7381 (vs E15 0.7183, +0.02 — small recovery)

## Interpretation
Mixing pure random in gives up too much overall to recover SKNSH eval_07.
Not a Pareto improvement. E15 remains the leader.

The SKNSH eval_07 tradeoff is intrinsic — there's no clean rescue.

## Next
E24: GC band + per-position EXACT balance via A↔T / C↔G swaps. Preserves
per-seq GC perfectly, eliminates per-position sampling noise (currently
~0.78% per position). Cleanly regularizes the marginal.
