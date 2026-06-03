# 016 — Dirichlet(0.3) + shuffled DHS combo

## What I built
25k Dirichlet(0.3) + 25k shuffled DHS (NMF-stratified, each seq internally shuffled).
Combines two best single strategies.

## Result
- eval_01 = 0.0776 (vs pure Dirichlet 0.0786, vs shuffled DHS 0.0754).
- mean ≈ 0.0946 (vs 0.0954 pure Dirichlet, 0.0919 shuffled DHS).

Better than pure shuffled DHS, worse than pure Dirichlet(0.3).

## Interpretation
**Combination DOES NOT win.** Even combining the two best single strategies dilutes the
better one (Dirichlet) without enough compensating gain from the other (bio).

Confirms strong pattern across all combo experiments:
- exp 002 (DHS + uniform Dirichlet): 0.0765
- exp 006 (top-DHS + Dirichlet(0.3)): 0.0780
- exp 016 (shuffled DHS + Dirichlet(0.3)): 0.0776

All combos cap around 0.077-0.078, all below pure Dirichlet(0.3) at 0.0786.

The composition signal Dirichlet(0.3) provides is **maximized when undiluted**. Any
biology added reduces the proportion of extreme-composition sequences and hurts eval_01.

## What to try next
The composition distribution shape matters. GC-stratification hurt (014 → 0.0775).
Maybe a SHARPLY TUNED bimodal Dirichlet matches the eval target better.

Exp 017: Dirichlet(0.3) with **rejection sampling**. Keep compositions whose GC is in
[0.3, 0.7]. Tests if MID-GC concentration helps (since GC-stratification with uniform
GC coverage hurt).
