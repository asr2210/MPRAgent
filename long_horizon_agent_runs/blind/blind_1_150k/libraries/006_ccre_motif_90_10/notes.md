# 006 — 90/10 cCRE + motif-embedded synthetic

## Goal
Test whether a small (10%) diversity component preserves the eval_08
super-additivity while paying only minimal dilution cost on other evals.

## Method
- 135,000 cCRE-centered windows (same pool as exp 003), seed=5
- 15,000 motif-embedded synthetic (random + 2 JASPAR consensus motifs), seed=5
- Concatenated, shuffled, written.

## Result: NEW BEST. Mean = 0.8883 (vs 0.8862 pure cCRE, +0.002).
Runtime 4716s.

| eval | cCRE   | 90/10  | Δ |
|------|--------|--------|---|
| 01   | 0.8288 | 0.8295 | +0.001 |
| 02   | 0.9270 | 0.9280 | +0.001 |
| 03   | 0.9193 | 0.9207 | +0.001 |
| 04   | 0.8657 | 0.8661 | +0.000 |
| 05   | 0.8285 | 0.8293 | +0.001 |
| 06   | 0.9274 | 0.9284 | +0.001 |
| 07   | 0.9025 | 0.9062 | +0.004 |
| 08   | 0.8922 | 0.9088 | **+0.017** |
| 09   | 0.9465 | 0.9452 | −0.001 |
| 10   | 0.9272 | 0.9287 | +0.002 |
| 11   | 0.8145 | 0.8155 | +0.001 |
| 12   | 0.7959 | 0.7984 | +0.003 |
| 13   | 0.9038 | 0.9022 | −0.002 |
| 14   | 0.9276 | 0.9285 | +0.001 |

## Key observations
1. **Wins on 12 of 14 evals**, two essentially flat (−0.001, −0.002).
2. **Eval_08 jumps +0.017** — the super-additive diversity effect survives at
   10% dilution. Not as big as at 50% (+0.024), but still substantial.
3. **Diversity gradient on eval_08:** 0% (cCRE) 0.892, 10% 0.909, 50% 0.916,
   100% (random) 0.908. Eval_08 peaks somewhere between 10% and 50%
   diversity, then declines at 100%. Suggests the model needs *both* cCRE-
   like and non-cCRE sequences to do its best on eval_08.
4. **Hard evals (11, 12)** barely moved: 0.815, 0.798. Diversity strategy
   has hit a ceiling for these.
5. **No regression anywhere** worth worrying about — −0.001 / −0.002 are
   within noise.

## Theory update (v5 → v6)
- Pure cCRE is NOT the ceiling. A small diversity component (~10%) gives
  free improvement on mean (+0.002) by lifting eval_08 ~+0.017 while
  costing < 0.002 elsewhere.
- The 90/10 ratio looks close to optimal for the cCRE + diversity family.
  Going to 95/5 or 99/1 is worth one experiment to see if even less
  diversity preserves enough eval_08 gain.
- Hard evals (11, 12) have plateaued. They need a fundamentally different
  lever. Strongest candidates:
    - Class-balanced cCREs (over-represent rare classes: PLS, CTCF-only).
    - Cell-type-specific cCREs.
    - Conservation-filtered cCREs.
    - Activity-extreme sequences.

## Next
Experiment 007 will test **class-balanced cCREs** (no diversity component).
Forces equal representation of PLS, pELS, dELS, and CTCF-only/DNase-H3K4me3
classes. The natural cCRE distribution is 74% dELS; balancing exposes the
model to relatively underrepresented promoter (PLS) and insulator (CTCF-only)
grammars. If this lifts hard evals 11/12, internal class diversity is a real
lever and we should pair it with the 10% diversity component in a follow-up.
