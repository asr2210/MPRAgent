# Skill: MPRA library design findings (live notes)

Running summary of what works and what doesn't for designing 150k-sequence
200bp MPRA libraries when measured against the project's 14 anonymous eval
sets. Updated as experiments accumulate.

## Aggregate leaderboard (mean across 14 evals)
| exp | design | mean | notes |
|-----|--------|------|-------|
| 001 | uniform random ACGT (50% GC) | 0.820 | floor — random is shockingly informative |
| 002 | hg38 chr19-22 200bp tiles    | 0.873 | +0.053 over random; eval_08 regresses |
| 003 | ENCODE cCRE-centered 200bp   | **0.886** | current best; +0.013 over genome |
| 004 | 50/50 cCRE + uniform random  | 0.880 | mean drops, eval_08 super-additive |
| 005 | 50/50 cCRE + motif-embedded  | 0.883 | tiny improvement over 004 |
| 006 | 90/10 cCRE + motif-embedded  | ?    | ratio refinement, in progress |

## What we know works
1. **Use real regulatory regions.** Random → genome → cCRE consistently lifts
   mean_r. Function enrichment is a real axis.
2. **A small diversity component helps eval_08 specifically.** Mixing any
   non-cCRE sequence type (random or motif-embedded) into a cCRE library
   produces super-additive performance on eval_08 (≈ 0.916, higher than
   either pure source). The effect is about *source diversity*, not motif
   content.
3. **Drop heavily soft-masked windows** (>50% lowercase) to avoid pure
   repeat tracks dominating the library. Keeps ~70% of cCRE-centered
   windows.

## What we know doesn't help
1. **50/50 dilution of cCRE costs more than it gains** on mean. Eval_08 wins
   ≈ +0.023; every other eval loses 0.005–0.015; net mean −0.003 to −0.006.
2. **Motif content alone doesn't drive the diversity bonus.** Uniform random
   and motif-embedded synthetic give essentially identical eval_08 lift
   (0.9156 vs 0.9160). Motifs only help on the non-eval_08 evals by a
   small amount.

## Eval-set behaviour
- **{eval_02, eval_06, eval_14}** report essentially identical numbers
  across all 5 experiments — treat as one effective signal.
- **{eval_01, eval_05}** likewise.
- **Hard evals: eval_11 (~0.81), eval_12 (~0.79).** Have not budged with
  function enrichment past cCREs or any diversity strategy. Need a
  different lever.
- **Easy evals: eval_08, eval_09, eval_14 (~0.92–0.95).** Saturated; small
  deltas there don't tell us much.
- **Eval_08 is the diagnostic for diversity.** It alone benefits from mixing
  sources; everything else benefits from pure-cCRE focus.

## Per-cell-line
Across all five experiments, HepG2 ≥ K562 ≈ SKNSH. No cell-line is
systematically lagging. No need to optimise for one over the others.

## Open hypotheses (untested)
- Class-balanced cCREs (equal PLS/pELS/dELS/CTCF) might beat the natural
  distribution (74% dELS).
- Conservation-filtered cCREs (high phyloP) might give higher information
  per sequence.
- Cell-type-specific cCRE subsets might capture different sequence grammars.
- Activity-stratified selection (if we could infer activity post-hoc) might
  balance the activity distribution.
- Dinucleotide-shuffled cCREs as diversity source — tests whether
  composition (k-mer) or motif identity drives the diversity benefit.
- 3-way mix: cCRE + non-cCRE genome + random — orthogonal diversity sources.

## Caveats / things to remember
- `prepare.py` takes ~50–80 min wall-time per run. Plan experiments to
  test the largest hypothesis difference per minute of compute.
- The model trained inside `prepare.py` is stochastic. Repeating the same
  library design with a different sample seed would give some noise — we
  haven't measured it, but inter-experiment deltas of >0.005 on mean are
  probably outside the noise.
