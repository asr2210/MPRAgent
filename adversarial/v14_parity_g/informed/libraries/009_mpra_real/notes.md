# 009_mpra_real

50,000 sequences from Sharpr-MPRA (Ernst et al. 2016), the gold standard
MPRA training dataset. Selected top 25k + bottom 25k by mean activity
(maximum measured activity contrast). Padded 145bp → 200bp with random ACGT
flanks (27bp left, 28bp right).

Result: eval_01 = -0.0002. Noise.

**Major implication**: Real, measured MPRA sequences with extreme activity
contrast produce no signal on the v14 evaluator. This is the most informative
null we have. v14 either:
  1. Uses a fixed model with no real correlation to MPRA activity, or
  2. Evaluates against a very narrow distribution that Sharpr does not cover.

Either way, optimizing for the metric is futile. From here we focus on
principled library design — making libraries that *would* train a strong
sequence-to-activity model in a real-world setting, regardless of what v14
reports.
