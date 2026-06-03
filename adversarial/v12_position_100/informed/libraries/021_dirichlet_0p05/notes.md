# 021 — Dirichlet(0.05) extreme alpha

## What I built
50k sequences from Dirichlet(0.05). Very corner-concentrated; many near-homopolymers.
Only 27,284 unique sequences (54% dedup rate) because alpha=0.05 produces compositions
near corners → many sequences are near-pure (e.g., 198 A's + 2 other) which collide.

## Result
- eval_01 = 0.0747 (vs alpha=0.1 0.0752, vs alpha=0.3 0.0786). Worse than 0.1 as predicted.

Final alpha sweep:
| alpha | eval_01 |
|------|--------|
| 0.05 | 0.0747 |
| 0.1  | 0.0752 |
| 0.3  | 0.0786 (PEAK) |
| 0.5  | 0.0768 |
| 1.0  | 0.0765 |

Sharp single peak at 0.3. Confirmed.

## Interpretation
- Going more extreme (smaller alpha) collapses unique sequence count and hurts.
- 27k unique sequences ≈ what we'd get if the harness dedupes some, AND the composition
  distribution becomes too degenerate (mostly near-homopolymers).

This confirms there's NO further composition-axis improvement to find.

## What to try next
Last orthogonal angles: reverse-complement augmentation (exp 022), fine alpha mix near
peak (exp 023). Then finalize with best library (exp 004 = Dirichlet(0.3)).
