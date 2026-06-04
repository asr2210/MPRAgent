# Experiment 020 — GTEX no quality filter

## Method
GTEX-only, no lfcSE filter (429K with valid ACGT 200bp). Random 50K.

## Result
**eval_01 = 0.0172** — drops vs lfcSE<0.7 (0.0222). Confirms a U-curve.

## GTEX quality U-curve
- lfcSE<0.3: 0.0145
- lfcSE<0.5: 0.0190
- lfcSE<0.7: 0.0222  ← peak
- no filter: 0.0172

So lfcSE<0.7 sits at the optimal tradeoff: nearly all of GTEX (92% pass) but
filters out the worst measurement noise. Going stricter biases toward easy
sequences; going looser admits too much label noise.

## Next direction
- Test if looser quality helps mixed Gosai too (021): tells us if GTEX-only
  matters or if loose-quality is the dominant trick
- Fine-grained search around 0.7 (later)
