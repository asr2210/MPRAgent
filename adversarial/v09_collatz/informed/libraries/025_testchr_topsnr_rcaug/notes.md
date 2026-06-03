# E025 — Top 25K SNR test-chr Gosai + their RCs (50K total)

## Result
eval_01 = 0.3321 (vs E022 = 0.3410, -0.009).

## Interpretation
RC augmentation HURTS when it replaces unique sequences. The library
lost 25K unique sequences for RCs of the same coordinates → less diversity.
The strand-invariance prior gives nothing.

Reconfirms: diversity > replication. The pipeline benefits from MORE
unique sequences, not from extra augmentation/training signal per seq.
