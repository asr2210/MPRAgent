# 028 — gene-density² at 500kb bins

## Result
**eval_01 = 0.5086** (3-seed). New peak by 0.0002 vs E24 (within noise).

## Design
Same as E24 but BIN_SIZE=500_000. 6,164 bins genome-wide.
Top-50 bin mass: 0.183 (vs 0.140 at 250kb, 0.064 at 1Mb).

## Bin-size sweep (squared weighting, 3-seed)
| Bin   | eval_01 |
|-------|---------|
| 1Mb   | 0.5071  |
| 500kb | 0.5086 ← peak |
| 250kb | 0.5084 |
| 100kb | 0.5063 |

Saturation between 250-500kb. Below 250kb degrades (over-fits small clusters).

## Interpretation
Bin-size optimization has saturated. The gene-density signal at this resolution
captures gene-rich vs gene-poor megabase neighborhoods, which correlates well with
the regulatory-active fraction of the genome.

## Next
E29: sharpen weighting by lowering EPS at 500kb. Push more mass to highest-density bins.
