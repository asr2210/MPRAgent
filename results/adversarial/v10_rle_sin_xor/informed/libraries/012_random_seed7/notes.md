# Experiment 012: random uniform seed 7

## Result
- eval_01: **0.5176** mean_r. K562=0.9945, HepG2=0.5626, SKNSH=-0.0043

## Comparison (variance across seeds)
- Seed 1:  0.5210
- Seed 7:  0.5176
- Seed 42: 0.5221
- Range: ~0.0045

## Interpretation
True variance across seeds is ~0.005. So 0.5176-0.5221 is the "random uniform"
band. To beat random_uniform clearly, an experiment needs eval_01 > 0.525.
