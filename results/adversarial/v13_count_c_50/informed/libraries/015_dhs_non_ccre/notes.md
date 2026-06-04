# Experiment 015 — DHS non-overlapping cCRE + cCRE class-balanced

## Design
25 k DHS uniform, restricted to summits NOT within 200 bp of any cCRE
midpoint + 25 k cCRE class-balanced. Forces the DHS half to contribute
orthogonal regulatory content.

## Key data finding
**73.5 % of DHS overlap cCRE within 200 bp.** All previous mixes
(008/011/012/014) were burning ~70 % of the DHS budget on positions
that were already represented in the cCRE half. Effective unique
genomic loci in 008 ≈ 31 k, not 50 k.

The DHS-non-cCRE pool is 950 k elements (26.5 % of DHS): open
chromatin sites without classical regulatory ChIP marks. Includes
insulators, weak/discovery-stage regulatory regions, cell-type-
specific accessibility not yet annotated.

## Result — BIG LIFT

| eval | 008 (best) | 011    | **015**          | Δ vs 008 |
|------|-----------:|-------:|-----------------:|---------:|
| 01   | 0.5671 | 0.5688 | **0.5736** | **+0.007** |
| 07   | 0.5833 | 0.5989 | **0.6131** | **+0.030** |
| 13   | 0.5603 | 0.5769 | **0.5924** | **+0.032** |
| 04/09| 0.5723 | 0.5579 | 0.5535 | −0.019 |
| 08   | 0.2111 | 0.1813 | 0.1492 | −0.062 |
| mean | 0.554  | 0.555  | **0.560**  | **+0.006** |

### Cell-type lift on eval_01
|         | 008    | **015**    | Δ      |
|---------|-------:|-----------:|-------:|
| K562    | 0.6227 | 0.6131 | −0.010 |
| HepG2   | 0.5350 | **0.5485** | **+0.014** |
| SK-N-SH | 0.5436 | **0.5592** | **+0.016** |
| mean    | 0.5671 | **0.5736** | **+0.007** |

## Interpretation
**The plateau was caused by sequence redundancy, not by a fundamental
ceiling.** Removing DHS-cCRE overlap broke through the 0.568 plateau
to 0.574.

**Specifically improved**:
- HepG2/SK-N-SH lifted significantly (+0.014/+0.016) — the new
  DHS-non-cCRE content is rich in non-K562 cell-type-specific
  accessibility (K562 is the most-studied biosample, so K562 elements
  are well-covered by cCRE; the leftover DHS skews toward other
  cell types).
- eval_07 and eval_13 jumped massively (+0.030, +0.032). These evals
  reward broad cell-type-spanning regulatory grammar — exactly what
  the orthogonal DHS pool provides.

**Trade-offs**:
- K562 on eval_01 dipped slightly (−0.010): the orthogonal DHS pool
  is K562-poor by construction.
- eval_04/09: −0.019. This eval seems to reward cCRE-heavy designs
  (012 still wins eval_04/09).
- eval_08: −0.062. Confirms eval_08 wants CpG-rich classical
  regulatory grammar that the orthogonal-DHS pool lacks.

## Theory update
- **Redundancy is the silent killer.** DHS and cCRE overlap ~75 %
  in 200 bp windows. Mixing them naively wastes budget on
  near-duplicates. Always deduplicate by coordinate.
- **Orthogonal-content sampling beats more-of-the-same sampling.**
  The 950 k DHS-non-cCRE pool is the highest-leverage source I've
  found in this dataset.
- The K562 ceiling and HepG2/SK-N-SH gap are *partly* compositional:
  K562 is over-represented in cCRE; orthogonal DHS lifts the
  underrepresented cell types.

## Numbers
mean_r: 0.560
eval_01: 0.5736 (new best, +0.007 over 008)
eval_07: 0.6131 (new best, +0.030)
eval_13: 0.5924 (new best, +0.032)
eval_08: 0.1492 (drop — confirms eval_08 wants CpG-rich)
time_s: 24

## Next
Push the orthogonal-content direction:
- 35 k DHS-non-cCRE + 15 k cCRE class-balanced (more orthogonal DHS,
  less redundant cCRE)
- Or, combine 015 with 011's component-targeting: DHS-non-cCRE
  component-weighted toward HepG2/SK-N-SH-relevant components.
