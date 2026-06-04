# 008 — DHS GC-content-stratified

## What I built
Meuleman 160k DHS pool. Computed GC content per sequence. Made 10 GC-quantile bins. Sampled 5000 from each bin → 50k DHS sequences with explicit uniform GC spread.

## Result
- eval_01 = 0.0737 — basically identical to DHS-NMF-stratified (0.0739, exp 001). No gain.
- mean across evals: 0.0901 — slightly worse than exp 001 (0.0909).

## Interpretation
DHS GC content has limited natural spread (min=0.10, max=0.885, mean=0.44, std=0.085) — even the "extreme" GC bins are mostly 30-55% GC because that's where the DHS pool concentrates. So GC-stratified DHS doesn't actually achieve much more compositional spread than NMF-stratified DHS.

**Key insight:** Real biological sequences are compositionally CONCENTRATED around natural genomic GC (~40-50% with some spread). To get the FULL compositional range of Dirichlet(0.3), you cannot use real DHS.

Combined with prior experiments:
- DHS variants: all hover at 0.074
- Dirichlet variants: best at 0.079
- Biology is fundamentally CAPPED at ~0.074 in this harness because of natural compositional concentration. Composition variance (synthetic) is the only way above 0.077.

## What to try next
Test a totally different lever: noise averaging via REPLICATION.

EXP 009: 5000 unique Dirichlet(0.3) sequences × 10 copies each = 50k total. If the prepare.py MPRA measures each occurrence (per-line) and the model averages duplicates, this reduces effective noise by ~3× per unique sequence. Risky — if the harness dedupes, we lose 90% diversity.

This tests whether the bottleneck at 50k is (a) per-sequence MPRA noise (replication helps) or (b) library diversity (replication hurts).

If replication wins, future experiments use it. If it loses, we abandon and pursue other axes.
