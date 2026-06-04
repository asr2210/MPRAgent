# 022 — 10 motifs/seq with overlap (density sweep midpoint)

## Plan
Brackets exp 020 (8) and 021 (15). Same 289-pool, overlap allowed.

## Result
**eval_01 = 0.4247.** Slightly below plateau but **SK-N-SH = 0.061** —
highest mean SK-N-SH across all experiments. Eval_06 and eval_11 = 0.0619,
eval_07 = 0.0638.

K562: 0.591 (slightly down from 020's 0.599).
HepG2: 0.622 (slightly down from 020's 0.629).
SK-N-SH: 0.061 (up from 020's 0.057).

## Density sweep
| Density | Overlap | mean_r | K562  | HepG2 | SK-N-SH |
|---------|---------|--------|-------|-------|---------|
| 3       | No      | 0.4283 | 0.596 | 0.629 | 0.060   |
| 5       | No      | 0.4259 | 0.595 | 0.625 | 0.055   |
| 8       | Yes     | 0.4284 | 0.599 | 0.629 | 0.057   |
| 10      | Yes     | 0.4247 | 0.591 | 0.622 | 0.061   |
| 15      | Yes     | 0.4193 | 0.589 | 0.619 | 0.050   |

## Insight (T16)
Different cell types have different optimal motif densities:
- K562/HepG2 want moderate density (~3-8).
- SK-N-SH wants higher density (~10), possibly because more motif content
  → more "biological-like" sequence → activates whatever feature the
  SK-N-SH oracle uses for noisy/complex sequences.

## Next
Exp 023 = combine: 10 motifs/seq overlap + broader pool (60/40 cell-type/
full JASPAR). Hypothesis: gives SK-N-SH both density and variety it
likes, while keeping K562/HepG2 in their preferred cell-type pool.
If mean_r > 0.43: combined design wins.
