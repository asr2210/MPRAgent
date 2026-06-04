# 011 — Pure top-signal DHS

## What I built
Sorted 160k DHS by total_signal; took top 50,000 (signal range 11.4-3599, median 27.8).

## Result
- eval_01 = 0.0745 (vs exp 001 random DHS = 0.0739, vs exp 004 Dirichlet(0.3) = 0.0786)
- mean ≈ 0.0913 (vs DHS 0.0901, Dirichlet 0.0954)
- eval_08 DROPPED 0.0716 → 0.0653 (worse than baseline DHS)

Marginally above baseline DHS, far below composition ceiling.

## Interpretation
High-signal selection adds ~0.0006 over random DHS — tiny. Biology+signal-selection is
**NOT a stronger lever than composition variance**. Top-signal DHS are
**compositionally more concentrated** (real strong regulators cluster near genomic
GC ~50%), so even with stronger labels per construct, the diversity loss costs net info.

eval_08 specifically dropped: this is the noisiest eval, where you'd most expect strong
signals to help. Instead they hurt. Suggests eval_08 cares about something other than
high-signal regulatory sequences.

## Hypothesis killed
H11: biological signal-selection breaks the composition ceiling → **REJECTED**

## What to try next
Test position-dependent composition structure (block sequences). Untested axis.
Exp 012: each sequence = two 100bp blocks with different Dirichlet(0.3) compositions.
