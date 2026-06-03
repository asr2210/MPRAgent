# 014 — CLEAN exact-match pool

## Plan
Replace substring matching with exact TF name matching. Curated list of
human cell-type TFs + universals. Tests whether cleaner pool beats exp 008.

## Result
**eval_01 = 0.4184.** Worse than exp 008 (0.4283) by -0.010.
K562: 0.593, HepG2: 0.628, SK-N-SH: 0.035.

176 PFMs in final pool (more than intended due to JASPAR dimer entries like
FOS::JUN, ETV5::FOXI1, FOXO1::ELK1 etc. — each ::-pair where one half
matched got included).

## What this teaches
- **Cleaner ≠ better.** Removing "noise" pollution (plant MYBs etc.) in
  exp 013 → polluted-but-broader pool (008) made K562/HepG2 similar but
  tanked SK-N-SH to 0.035 (worst SK-N-SH of all experiments).
- Sub-pattern: K562 went UP (0.593 vs 0.596 in 008), HepG2 went UP slightly
  (0.628 vs 0.629). SK-N-SH crashed (0.035 vs 0.060). The "noise" PFMs in
  exp 008 were apparently helping SK-N-SH somehow.
- Possibly: exp 008's broader pool included substring matches like
  ATMYB31, KLF10-17 that look RANDOM-LIKE — and SK-N-SH might benefit from
  random-like inserts mixed in (since the neural oracle may be sensitive
  to baseline composition rather than specific motifs).
- Or: exp 014's pool is *too* curated and over-represents shared
  patterns (lots of bZIP/ETS dimers) at the expense of variety.

## Theory T11
For mean_r, **pool VARIETY matters more than pool purity**. The surrogate
learns from a wider distribution of motif patterns when the pool is broad,
even if some patterns are biologically irrelevant. SK-N-SH is especially
sensitive — it crashed with the curated pool, suggesting its oracle
benefits from a diverse motif background, not just neural-specific motifs.

## Next
Exp 015: variable motif density. Per-sequence motif count drawn from
uniform[1, 6]. Same 289-pool from exp 008. Tests whether sequence-level
density diversity adds information beyond identity diversity.
