# Experiment 019 — Real cCREs, NO GC filter

## Design
Real ENCODE cCREs as 200bp centered windows with NO composition filter.

## Stats
- GC: mean=0.485, std=0.098, P10=0.37, P90=0.62
- Note: ~80% of sequences are in [0.37, 0.62] — natural-like range with
  a long tail. This is fundamentally DIFFERENT from 016's bimodal extremes
  (100% at 0.30 or 0.70).

## Result
- eval_01 = **0.3953** (Δ vs 007 = -0.002, within noise)
- mean_r = **0.3841** (Δ vs 007 = -0.002, within noise)
- Per-cell: K562=0.609, HepG2=0.431, SK-N-SH=0.146 (slight SKNSH bump)

## Interpretation
**Natural cCREs without GC filter match GC-filtered random.** The natural
composition variation in real enhancers is acceptable — most sequences
remain in the natural-like range, with only a small tail outside.

This RECONCILES with the 016 bimodal failure: 016 forced ALL sequences
to compositional extremes, whereas natural cCREs have a Gaussian-like
distribution centered at the natural range.

**Marginal SK-N-SH improvement.** SKNSH r=0.146 vs 007's 0.141. Real
cCREs include some neural enhancers, which may help marginally —
though still within noise.

## What this rules out
- That GC filtering is essential for cCRE-derived libraries
- That natural composition variance hurts (only extreme matters)

## What this reinforces
- The natural-like subspace is a BROAD plateau, not a narrow ridge —
  any distribution mostly within [0.40, 0.60] per-seq GC works similarly.
