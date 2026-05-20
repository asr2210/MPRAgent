# 010 — cCRE neighbourhoods (center + ±200bp flanks)

## Goal
Test whether sequences adjacent to cCREs carry additional generalisable
regulatory signal beyond cCRE centres themselves.

## Method
- Sample cCREs in random order; for each, try to extract 3 windows: center,
  upstream (mid−300:mid−100), downstream (mid+100:mid+300).
- Apply per-window filter (ACGT only, <50% softmasked).
- Stop when 150k windows collected. Used 72,925 cCREs; 68,774 windows
  failed the filter (mostly >50% softmasked, i.e. repeat-rich flanks).

## Result: REGRESSION. Mean = 0.8822 (vs 003 pure cCRE 0.8862, **−0.004**).

| eval | cCRE | neighborhoods | Δ |
|------|------|---------------|---|
| 01 | 0.8288 | 0.8249 | −0.004 |
| 02 | 0.9270 | 0.9237 | −0.003 |
| 03 | 0.9193 | 0.9149 | −0.004 |
| 04 | 0.8657 | 0.8607 | −0.005 |
| 05 | 0.8285 | 0.8248 | −0.004 |
| 06 | 0.9274 | 0.9240 | −0.003 |
| 07 | 0.9025 | 0.9018 | −0.001 |
| 08 | 0.8922 | 0.8844 | **−0.008** |
| 09 | 0.9465 | 0.9395 | −0.007 |
| 10 | 0.9272 | 0.9196 | −0.008 |
| 11 | 0.8145 | 0.8109 | −0.004 |
| 12 | 0.7959 | 0.7926 | −0.003 |
| 13 | 0.9038 | 0.9044 | +0.001 |
| 14 | 0.9276 | 0.9241 | −0.003 |

## Key observations
1. **Flanks dilute the signal.** Loses on 13 of 14 evals.
2. **Eval_08 loses the most non-trivial amount** (−0.008) — and no eval_08
   bonus since there's no synthetic component.
3. **Hard evals barely move** (−0.004, −0.003). Adding regulatory context
   doesn't unlock them either.
4. **cCRE centres are the right granularity.** Three regressions now confirm
   this: class balance (007 −0.004), non-cCRE genome (008 −0.003),
   neighbourhoods (010 −0.004). The natural cCRE catalogue at 200bp centres
   is genuinely the highest-information-per-sequence regulatory pool I can
   easily build.

## Theory update (v9 → v10)
- The cCRE-centre catalogue at 200bp is near the per-sequence information
  optimum among "biologically motivated" libraries.
- The remaining headroom (above 0.888) is unlikely to come from sampling
  *more* cCRE-related sequences. It will come from one of:
    - Higher information density per sequence (conservation filter,
      cell-type filter, "high-quality" cCRE subsets)
    - Data augmentation (RC, mutations, k-mer shuffles) — train more on
      the same sequences
    - Multi-source diversity (3-way mixes that hit different aspects of
      the test distribution)
    - Orthogonal data sources (FANTOM enhancers, DNase peaks, ChIP-seq)

## Next
Going with experiment 011 = **cCRE + mutated cCRE pairs** (75k unique
cCREs + 75k of the same cCREs with 10% random base substitutions). Tests
whether data augmentation via small perturbations gives the model
counter­factual training signal that beats raw cCRE diversity alone.
Compare directly to exp 003 (150k pure cCRE) to isolate the augmentation
effect.
