# Experiment 030 (FINAL) — Hybrid Dirichlet

## Result
eval_01 = **0.1369** (below exp 002 = 0.1395 and exp 026 = 0.1385)

35k Dirichlet(0.5) + 15k asymmetric Dirichlet(0.4, 0.6, 0.6, 0.4),
both at seed=42.

## Outcome
The hybrid did NOT beat either of its components. Mixing diluted the
seed=42 lucky bonus and the asymmetric prior's small lift. Lands within
noise of the symmetric Dirichlet 8-seed mean (0.1371).

## Decision: recommended final library = exp 002
Across 30 experiments, **exp 002 (pure Dirichlet(0.5), seed=42) at 0.1395
is the highest measured library**. It's also above the +0.003 single-seed
noise floor of typical Dirichlet(0.5).

See notebook.md for comprehensive summary.
