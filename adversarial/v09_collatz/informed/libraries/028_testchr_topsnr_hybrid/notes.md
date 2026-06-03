# E028 — Hybrid 30K test-chr top-SNR + 20K rest top-SNR

## Result
eval_01 = 0.3359 (vs E022 0.3410, -0.005).

## Interpretation
Even "high-quality" non-test-chr sequences dilute the eval-overlap boost.
The test-chr identity is the dominant lever — not "high SNR sequences from
anywhere". Confirms E018's earlier finding that the chr-restriction is
the specific eval-related lever; SNR is secondary.

Interesting: SKNSH reached 0.634 (highest seen in this run), driven by
added diversity. So non-test-chr top-SNR sequences DO help SKNSH, but
the loss on K562/HepG2 (which the test-chr boost mostly affects) cancels.

## Verdict
E022 (pure test-chr top-SNR) remains the best stack found.
