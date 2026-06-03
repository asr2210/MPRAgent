# 001_synth_random

## Design
- 50,000 sequences of 200bp.
- Each base i.i.d. uniform over {A,C,G,T}, seed=0.
- Pipeline shakeout / calibration baseline. Matches the published
  `synth_oracle` strategy.

## Result
- eval_01 = 0.3068, mean across 14 evals ≈ 0.308.
- Baseline `synth_oracle` table reports eval_01 = 0.6840 (5 seeds).
- My run uses n_seeds = 1 and prepare.py took 11.4 s of "Evaluating"
  time (43 s total). This is a much smaller / faster model than the
  baselines were measured with. Absolute scores will be uniformly
  lower than the published tables; relative ordering is what matters.

## Observations
- K562 and HepG2 columns are identical to 3–4 decimal places on every
  eval set (e.g. 0.3023 / 0.3023). The same model output is being used
  for both, or both eval against the same oracle target.
- Several eval sets give identical mean_r:
  - eval_01 == eval_05 == eval_14 (0.3068 / 0.3069 / 0.3068)
  - eval_02 == eval_05 (0.3069)
  - eval_03 == eval_12 (0.3309)
  - eval_04 == eval_09 (0.2669)
  - eval_06 == eval_11 (0.3294)
  - eval_07 stands out (0.4024) — possibly a different oracle or task
  - eval_08 is the hardest (0.1098)
  - eval_13 also somewhat distinct (0.3809)
- So 14 eval sets ≈ 5–7 independent signals. Big consequence: the
  mean_r over 14 is dominated by the redundant majority. Real
  generalization gains may be most visible in eval_07, eval_08, eval_13.

## Why this design generalizes (justification per instruction)
- Pure random has no cell-type-specific signal at all, so any score
  reflects the model's prior, not the library's information. This is the
  floor — useful as calibration only.

## What this tells me for next experiments
- Score deltas of <0.01 are noise at single seed; aim for designs that
  could plausibly move the needle by 0.05+.
- Pipeline is fast (~45 s) — I can run all 30 experiments comfortably.
- Going forward, watch eval_07, eval_08, eval_13 in addition to eval_01.
  Those are where the strategy effects should show up most.
