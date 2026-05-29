# 021 — 95/5 cCRE/synthetic ratio

## Goal
Test if reducing synthetic from 10 % → 5 % (and using slots for more
unique cCREs) helps. Maybe synthetic effect saturates at low %.

## Method
- 71,250 cCRE-fwd + 71,250 different cCRE-RC (142,500 unique, mixed strand)
- 7,500 motif-embedded synthetic
- Total 150k.

## Result: regression. Mean = 0.8878 (vs 015 0.8905, −0.003).

| eval | 015 | 021 | Δ |
|------|------|------|---|
| 08 | 0.9115 | 0.9038 | **−0.008** |
| 11 | 0.8187 | 0.8156 | −0.003 |
| 12 | 0.8011 | 0.7986 | −0.003 |

## Key observations
1. **Eval_08 dropped most (−0.008)** — confirms synthetic is needed at
   sufficient quantity to give the eval_08 bonus. 5 % is too little.
2. **Other evals each dropped ~−0.002** — the 7.5k extra unique cCREs
   don't make up for the loss of 7.5k synthetic.
3. **10 % synthetic is the optimum** on the ratio axis (with RC bulk):
   - 0 % (013): 0.8832 — too little synthetic, eval_08 = 0.888
   - 5 % (021): 0.8878 — slight loss vs 10 %
   - 10 % (015): 0.8905 — optimum
   - 50 % (005): 0.8829 — too much dilution

## Theory update (v20 → v21)
The 90/10 ratio is the optimum within the cCRE-RC-synthetic family.
Both ends (5 % and 50 % synthetic) underperform.

## Next
EXPERIMENT 022 = critical NOISE/REPLICATION test. Repeat 015 with a
different seed (seed=21 instead of seed=14). Tells me how much of the
small deltas I've been chasing is real signal vs training stochasticity.
This is essential context for interpreting the rest of the experiments
and for picking a final best library.
