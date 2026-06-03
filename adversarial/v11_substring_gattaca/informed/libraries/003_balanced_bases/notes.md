# E3: balanced_bases

## Design
- Each of 50,000 sequences = random permutation of exactly 50A + 50C + 50G + 50T.
- Strict per-base balance. Single seed.

## Result (mean across 14 evals ≈ 0.797)
- eval_01: 0.8150 (vs E2 random_uniform 0.8565, **-0.042**)
- mean across 14 evals: 0.797 (vs E2 0.841, **-0.044**)

## Cell-type breakdown vs E2 (random_uniform)
| eval | E2 K562 / HepG2 / SKNSH | E3 K562 / HepG2 / SKNSH | Δ SKNSH |
|------|--------------------------|--------------------------|---------|
| 01 | 0.85 / 0.87 / 0.85 | 0.85 / 0.91 / 0.68 | -0.17 |
| 07 | 0.78 / 0.90 / 0.72 | 0.79 / 0.94 / 0.35 | **-0.37** |
| 13 | 0.77 / 0.89 / 0.82 | 0.78 / 0.93 / 0.50 | -0.32 |

K562 unchanged. HepG2 **gained** ~0.05 across evals. SK-N-SH **collapsed** by 0.15-0.37.

## Interpretation — strong update to theory
Constraining per-sequence base counts to exact equality:
- Helps HepG2 (maybe HepG2 eval is uniformly composed too).
- DESTROYS SK-N-SH performance — many SK-N-SH eval sequences must depend
  on compositional features that don't exist when bases are perfectly
  balanced (e.g., AT-rich neural enhancers, CpG islands, dinucleotide
  bias regions).

**Compositional diversity within sequences is critical for at least one cell
type's eval.** Over-constraining composition is a worse error than letting
it vary (per random_uniform). The slight published lift gc_50 vs
random_uniform (~0.003) does NOT extrapolate to "tighter is better."

## Theory v2 → v3
- v2 was "match the eval distribution as tightly as possible."
- v3: "match the eval CENTRAL distribution (50% GC IID) and PRESERVE
  compositional variance — different cell types appear to query different
  compositional regions of sequence space."
- The eval set is heterogeneous across cell types in compositional terms.
- A good library needs to cover SK-N-SH-relevant compositions (probably
  including some AT-rich or compositionally varied sequences).

## Implications for next experiment
- Adding compositional variance (controlled) might HELP SK-N-SH while
  not hurting K562/HepG2 much.
- But too much variance (dirichlet, gc_sweep) is destructive.
- Target: random_uniform + a small fraction of compositionally varied
  sequences in a moderate band (e.g., 35-65% GC).
