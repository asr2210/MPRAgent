# E26: per_position_stratified

## Design
Each position p built as shuffled column of exactly 12500 of each base.
Per-position p exactly uniform. Per-seq GC = Binomial(200, 0.5) — same
as random_uniform.

## Result
- eval_01: 0.8570 (vs E2 0.8565, ~tied; vs E15 0.8777, -0.021)
- mean_r: ~0.8460 (vs E2 0.8408, +0.005; vs E15 0.8591, -0.013)
- SK-N-SH eval_07: 0.6783 (vs E2 0.7723, -0.09)

## Interpretation — attribution
Per-position exact balance contributes only +0.005 mean over random.
**GC filter alone (E15) contributes +0.018**, the much larger effect.
Combining (E24) hits the same E15 ceiling (~0.859), so no extra gain.

**Conclusion: GC filtering [90, 110] is the only meaningful improvement.**
Per-position balance is a noise-reduction afterthought; biology, motifs,
diversity, RC pairs all neutral or hurt.

## Next
E27: tighter GC band [91, 109] with per-pos balance — final attempt to
squeeze out a small win by combining knobs.
