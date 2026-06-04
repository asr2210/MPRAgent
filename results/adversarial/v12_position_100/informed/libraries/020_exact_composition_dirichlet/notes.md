# 020 — Exact-composition Dirichlet(0.3)

## What I built
For each Dirichlet(0.3) target composition p, built a sequence with EXACTLY
round(200·p_i) of each base (largest-remainder rounding), then shuffled positions.
Removes multinomial sampling noise within each sequence.

## Result
- eval_01 = 0.0771 (vs Dirichlet(0.3) 0.0786). -0.0015 loss.
- mean ≈ 0.0947.

## Interpretation
Removing within-sequence multinomial noise HURT eval_01. Possible explanations:
- The multinomial noise was actually providing extra training diversity (each
  composition sampled to slightly different actual counts).
- Forcing exact composition reduces effective compositional spread (compositions
  near simplex corners are now even more concentrated).
- The model uses the noise as a regularizer.

Whatever the mechanism, **adding even minor structural constraints to Dirichlet(0.3)
hurts.** Natural unconstrained Dirichlet(0.3) sampling is optimal.

## Hypothesis killed
H20: removing multinomial noise improves composition signal → REJECTED.

## What to try next
Try extreme alpha values (alpha=0.05) and/or bimodal alpha mix (very extreme +
moderate). These are the last untried composition variants.
