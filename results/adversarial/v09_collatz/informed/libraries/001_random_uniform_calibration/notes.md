# E001 — Random uniform calibration

## Design
50,000 sequences × 200bp, each base i.i.d. uniform {A,C,G,T}, seed=1.

## Result
eval_01 = **0.2444** (baseline `synth_oracle` claims 0.6840 at 5-seed mean)
mean across 14 evals = 0.227
Time: 26s

## Per-cell-type breakdown (eval_01)
- K562 r = 0.1368
- HepG2 r = -0.0413 (anti-correlated!)
- SKNSH r = 0.6377 (surprisingly high)

eval_08 mean = 0.0850 — the lowest, opposite of baseline `synth_oracle` (0.7696). 
Strong evidence eval_08 is NOT random/synthetic-like — it scored well for `synth_oracle`
because that strategy used ORACLE labels (clean), while my pipeline uses REAL noisy MPRA labels.

## Key interpretation
**The baseline table in instructions.md uses ORACLE LABELS** (clean simulated activity)
for everything except `mpra_real`. My pipeline gives real (noisy) measurements via
prepare.py. For high-signal sequences (genomic regulatory regions), real-vs-oracle
gap is ~0.06 (mpra_real 0.6026 vs mpra_oracle 0.6643). For LOW-signal sequences
(random), the gap is huge: 0.6840 (oracle) → 0.2444 (real). Random sequences have
no biological activity so real measurements are mostly noise → unlearnable.

**Implication**: my real-label pipeline benchmark for random is ~0.24 floor, not 0.68.
The relevant comparison baselines for me are likely the `mpra_real` style numbers,
i.e. ~0.6 range for high-signal sequences.

## Next
Get real-labeled DHS-based baseline to calibrate the genome-derived strategies.
