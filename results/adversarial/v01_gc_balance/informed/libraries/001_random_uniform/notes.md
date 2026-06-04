# 001_random_uniform — notes

**Purpose**: smoke test. Validate pipeline, measure runtime, get a feel for `result.json` schema.

**Method**: 50k × 200bp i.i.d. uniform random over {A,C,G,T}, seed=0.

**Result**: eval_01 = 0.4704 (close to baseline `random_uniform` reported as 0.494). Wall clock 48s, eval portion 18.1s.

**Key observations from the result table**:
- Several eval sets are numerically identical:
  - eval_01 ≈ eval_02 ≈ eval_05 ≈ eval_14 (all ~0.470)
  - eval_03 ≈ eval_12 (0.4688)
  - eval_04 ≈ eval_09 (0.3947)
  - eval_06 ≈ eval_11 (0.4684)
  This implies ~8 distinct evaluation signals, not 14. eval_07, eval_08, eval_10, eval_13 are each unique. **eval_01 is the primary metric anyway**, so I'll mostly focus on that, but track distinct sets too.
- eval_08 is by far the hardest (0.158); even rich strategies only get to 0.70.
- SK-N-SH within-cell-type correlations are systematically lower than K562/HepG2 — SK-N-SH is the harder cell.

**Compute budget**: 48s × 30 experiments = ~24 min of pipeline time. Plenty of headroom for library generation work (downloads, motif scans, etc.).

**Theory update**: no surprises. Random is a floor. Real work begins exp 002.
