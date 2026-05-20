# 019 — Homotypic motif cluster synthetic

## Goal
Test whether biologically-motivated synthetic architecture (3 copies of
same JASPAR motif at fixed spacing) beats random-placement motif synthetic.

## Method
- 67,500 cCRE-fwd + 67,500 cCRE-RC of different (135k unique, mixed strand)
- 15,000 synthetic: random ACGT background + 3 copies of the SAME JASPAR
  motif at positions ~50, 100, 150 (30bp anchor spacing).
- Motif sampled uniformly from JASPAR motifs ≤ 30bp (877 motifs).

## Result: tied / marginal regression. Mean = 0.8893 (vs 015 0.8905, −0.001).

| eval | 015 | 019 | Δ |
|------|------|------|---|
| 01 | 0.8330 | 0.8329 | 0.000 |
| 02 | 0.9303 | 0.9290 | −0.001 |
| 03 | 0.9230 | 0.9214 | −0.002 |
| 04 | 0.8657 | 0.8654 | 0.000 |
| 05 | 0.8328 | 0.8327 | 0.000 |
| 06 | 0.9307 | 0.9294 | −0.001 |
| 07 | 0.9063 | 0.9053 | −0.001 |
| 08 | 0.9115 | 0.9078 | −0.004 |
| 09 | 0.9467 | 0.9450 | −0.002 |
| 10 | 0.9313 | 0.9289 | −0.002 |
| 11 | 0.8187 | 0.8188 | 0.000 |
| 12 | 0.8011 | 0.7999 | −0.001 |
| 13 | 0.9049 | 0.9044 | −0.001 |
| 14 | 0.9308 | 0.9295 | −0.001 |

## Key observations
1. **Homotypic clusters ≈ random-placement** — essentially tied (-0.001
   mean, within noise).
2. **Eval_08 dropped the most** (-0.004). Structured synthetic less
   "random-like", weaker eval_08 bonus — consistent with earlier pattern
   (009 dense motif had the same direction).
3. **Architecture of the synthetic component doesn't matter much**.
   Whether 2 motifs at random positions, 5 motifs at random positions,
   3 copies of same motif at fixed positions — all give ~0.890.

## Theory update (v18 → v19)
- The synthetic component is fungible. What matters is "presence of
  random-like sequences in training", not the specific motif structure
  or count.
- The 0.890 mean ceiling is robust to ALL tested variations of the
  synthetic component:
    - 2 motifs random placement (006/015): 0.888/0.890
    - 5 motifs random placement (009): 0.889
    - Homotypic clusters (019): 0.889
    - Uniform random no motifs (004 50/50): different ratio, but synthetic-as-source proven equivalent
    - FANTOM5 partial substitute (016): 0.891

## Next
Time to test cCRE-pool quality variants. Specifically: **wide cCREs (≥300bp)**
might be higher-information per sequence (more complex regulatory regions).
Pool size: 488k cCREs ≥300bp wide (46% of total).

EXPERIMENT 020 = wide-cCRE bulk + 015 recipe. 67.5k wide-cCRE-fwd +
67.5k wide-cCRE-RC of different + 15k motif = 150k.
