# 002 — hg38 chr19-22 tiled

## Goal
Test whether evolution-shaped sequence space teaches the model more general
regulatory grammar than uniform random.

## Method
- hg38 chr19, 20, 21, 22; pooled non-overlapping 200bp windows
- Drop windows with N or any non-ACGT char
- Drop windows with > 50% soft-masked (lowercase, i.e. RepeatMasker) bases
- Total pool after filtering: 447,739 windows
- Sample 150,000 uniformly (seed=1)

## Result vs random baseline
Mean across 14 evals: **0.8731** (vs 0.8200 for random, **+0.053**). Runtime 3783s.

| eval | random | genome | Δ      |
|------|--------|--------|--------|
| 01   | 0.7760 | 0.8182 | +0.042 |
| 02   | 0.8718 | 0.9162 | +0.044 |
| 03   | 0.8562 | 0.9080 | +0.052 |
| 04   | 0.8178 | 0.8567 | +0.039 |
| 05   | 0.7756 | 0.8181 | +0.042 |
| 06   | 0.8722 | 0.9166 | +0.044 |
| 07   | 0.7910 | 0.8882 | **+0.097** |
| 08   | 0.9080 | 0.8679 | **-0.040** |
| 09   | 0.8890 | 0.9364 | +0.047 |
| 10   | 0.8670 | 0.9062 | +0.039 |
| 11   | 0.7623 | 0.8044 | +0.042 |
| 12   | 0.7394 | 0.7857 | +0.046 |
| 13   | 0.7762 | 0.8947 | **+0.119** |
| 14   | 0.8722 | 0.9164 | +0.044 |

## Key observations
1. **Real genome wins on 13 of 14 evals**, by 0.04–0.12. Big effect.
2. **Eval 08 is the lone regression** (-0.040). Something about that eval set
   prefers random sequences. Hypotheses: eval 08 may contain random-like
   negative controls or may emphasise extreme high/low activity that a
   random library covers better.
3. **Eval 13 jumped most** (+0.119) — probably tests something genomic-context-
   dependent (chromatin / TFBS spacing / dinucleotide structure).
4. **Eval 07 also big** (+0.097) — same story likely.
5. **Eval clustering is now confirmed empirically:**
   - {02, 06, 14}: 0.9162 / 0.9166 / 0.9164 — three reads of the same test
   - {01, 05}:    0.8182 / 0.8181 — two reads of the same test
   So 14 evals = ~11 effective signals.
6. Per-cell-line: K562 0.873, HepG2 0.877, SKNSH 0.870. Roughly balanced, no
   cell-line bias from training on chr19-22.

## Theory update
- Real genomic sequence space carries large additional information vs random.
  Confirms the obvious-but-unverified hypothesis.
- But "more realistic ≠ universally better" — eval 08 says some part of the
  test distribution is *poorly* covered by genomic tiling. A library that
  fully replaces random with genome loses ground there.
- The right library is probably **not** "as genomic as possible" but
  "covers more axes of the activity distribution than any single source can".

## Next
Two strong candidates for experiment 003:
- **A. ENCODE cCREs** — does enrichment for putative regulatory elements
  give further gains beyond random tiling, or does it overfit to the
  cell-types ENCODE was profiled in?
- **B. 50/50 random + genome hybrid** — does a mixed library recover the
  eval-08 strength while keeping the genome wins?

Going with A (cCREs). It tests the "more enriched = better" hypothesis at a
finer grain than 001→002 already established. If cCREs underperform tiles,
it suggests random untargeted genome is doing more than 'just' regulatory
enrichment (e.g. context, repeats, background composition all matter).
B is the natural follow-up if A is unclear.
