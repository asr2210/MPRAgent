# 022 — Dirichlet(0.3) + reverse complements

## What I built
25k Dirichlet(0.3) sequences + 25k of their RCs = 50k shuffled.

## Result
- eval_01 = 0.0775 (vs Dirichlet(0.3) 0.0786). -0.0011.
- mean ≈ 0.0951 (vs 0.0954).

## Interpretation
RC augmentation hurt slightly. RC of (pA, pC, pG, pT) gives compositions (pT, pG, pC, pA)
— different specific base allocations even if entropy is the same. If the harness model
distinguishes A from T at the count level (likely), RC pairs introduce noise: same
activity expected for compositionally-permuted sequences but the model has to learn
strand-symmetry from limited examples.

Net: RC augmentation = slight loss.

H22 → REJECTED.

## What to try next
Fine alpha mix near peak (16.6k each at alpha=0.2, 0.3, 0.4).
