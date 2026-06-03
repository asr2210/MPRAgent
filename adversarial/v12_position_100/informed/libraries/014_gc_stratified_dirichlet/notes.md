# 014 — GC-stratified Dirichlet(0.3)

## What I built
Oversampled 200k Dirichlet(0.3) compositions, binned by GC into 10 quantile bins,
selected 5k from each bin → 50k Dirichlet(0.3) compositions with FORCED uniform
GC coverage. One 200bp sequence per composition.

## Result
- eval_01 = 0.0775 (vs exp 004 random Dirichlet(0.3) = 0.0786). **-0.0011 loss**.
- mean ≈ 0.0948 (vs 0.0954). Approximately neutral.
- eval_08: 0.0700 (vs 0.0716). Slight loss.

## Interpretation
Forcing uniform GC coverage of Dirichlet(0.3) compositions did not help — slightly hurt.

Possible reason: Dirichlet(0.3) natural draws already over-sample mid-GC compositions
(due to multinomial entropy peak at uniform), which MATCHES what the evals' target
sequence distribution likely contains. Forcing extra mass at extreme-GC bins
(GC ≈ 0.0 or 1.0) wastes capacity on sequences the eval doesn't care about.

This refines theory v4: **the model wants composition distribution to MATCH the eval
target distribution, not maximize spread.** Dirichlet(0.3) hits this naturally;
forcing flatter GC distribution moves away from the target.

## Hypothesis killed
H14 (GC stratification of Dirichlet improves coverage) → REJECTED.

## What to try next
Alpha refinement around 0.3. Test Dirichlet(0.5) to see if peak is slightly above 0.3
or right at 0.3.
