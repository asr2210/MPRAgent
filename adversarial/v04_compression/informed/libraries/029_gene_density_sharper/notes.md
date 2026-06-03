# 029 — sharper gene-density (EPS=0.1)

## Result
**eval_01 = 0.5083** (3-seed). Essentially tied with E28 (0.5086).

## Design
Same as E28 (gene-density², 500kb bins) but EPS=0.1 instead of 0.5.
Top-50 bin mass: 0.190 (vs 0.183 at EPS=0.5).

## Interpretation
Lowering EPS sharpens the weighting (more mass on top bins, less on
gene-deserts) but moves the score essentially zero. The peak at
~0.508-0.509 is robust to EPS in [0.1, 0.5].

The gene-density signal has SATURATED — no axis (bin-size, exponent,
EPS) moves the score beyond ~0.509.

## Next
E30: final stacking attempt or fill exponent gap at new peak bin.
