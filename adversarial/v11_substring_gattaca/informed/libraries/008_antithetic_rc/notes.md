# E8: antithetic_rc

## Design
25,000 IID random uniform sequences + their 25,000 reverse complements.
Total 50k. Preserves per-position IID 50% GC distribution exactly.

## Result (mean across 14 evals ≈ 0.839)
- eval_01: 0.8546 (vs E2 0.8565, **-0.002**)
- Mean: 0.839 (vs E2 0.841, -0.002)
- SK-N-SH: ~0.84 (matches E2 — NO crash!)

## Interpretation
- Essentially matches E2 random_uniform.
- The model likely has internal RC augmentation (no lift from explicit RC pairs).
- Halving unique sequences (25k vs 50k unique) costs essentially nothing
  (-0.002, within seed noise).
- Per-position IID preservation = SK-N-SH stays healthy. Confirms the
  SK-N-SH sensitivity is to per-position distribution shifts, not to
  sequence-level transformations.

## Implication
- The model is robust to within-IID variations like halving unique seqs.
- Effective sample size at 50k random is near saturation.
- The "marginal value of an additional unique random sequence" is very
  low at this scale.
- Therefore I should NOT try to "improve diversity" via more unique
  sampling — that won't move the needle.

## Strong implications for generalization
- The model is data-efficient on random sequences. Even 25k unique IID
  random captures most of the learnable signal.
- For cross-cell-type generalization, this suggests the eval set is
  near-random (no cell-type-specific structure preferred). A model
  trained on random data should transfer well to any cell type whose
  eval is also random — but may NOT transfer to a cell type whose eval
  is biological/structured. So a random library is "good" for synthetic
  evals but may underperform on biology evals.

## Next
E9: random_uniform with a different seed (seed=1) to characterize
single-seed variance. If E9 mean ≈ 0.841 ± 0.005, baseline noise is ~0.005
and E2's 0.841 is reliable. If E9 differs more, single-seed comparisons
are noisy and I need bigger effect sizes to claim improvements.
