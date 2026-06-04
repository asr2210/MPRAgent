# Experiment 008 — Gosai 2024 MPRA library (real biology)

## Result
eval_01 = 0.5031. K562=0.92, HepG2=0.56, **SK-N-SH=0.026** (first nontrivial signal!).

## Verdict
Real biological sequences from a large variant-centered MPRA dataset:
- Hurt K562 (0.99 → 0.92) like other biological libraries
- HepG2 unchanged (0.56)
- **SK-N-SH lifts off zero** for the first time (0.026)

Does NOT reach the Table 1 mpra_oracle baseline of 0.6643 — confirming Table 1
is not from this pipeline.

## Key finding
SK-N-SH r CAN be moved off zero, but only by sequences derived from real
variant-centered regulatory contexts. DHS regions (open chromatin) did NOT help.
The relevant property is likely "sequence is centered on a known functional
genetic variant" or "sequence has measured non-trivial activity dynamic range."
