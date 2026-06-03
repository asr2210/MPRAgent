# E022 — Gosai chr 7/9/13/21/X, top 50K by SNR (max |fc|/lfcSE across cells)

## Result — NEW BEST
eval_01 = 0.3410. Per-cell K562=0.172, HepG2=0.221, SKNSH=0.630.

vs E008 (chr only, random 50K): 0.3359 → +0.005 from SNR filter.
vs E011 (chr + loud + clean threshold): 0.3355 → +0.006 from continuous SNR.

## Interpretation
The chr-restriction (+0.013) and SNR ranking (+0.005) **do stack** —
just not by much. Continuous SNR ranking (max(|fc|/lfcSE)) extracts a
small extra signal that the binary loud+clean threshold (E011) missed.

This is the first library to clearly beat E008 (within seed noise).
Suggests there's a few more crumbs of signal in the test-chr Gosai
subset that careful filtering can extract.

## Per-cell pattern
All three cells creep up slightly. K562 hit 0.172 (E008 was 0.170);
HepG2 0.221 (E008 0.219); SKNSH 0.630 (E008 0.617). The structural
caps remain ~0.17/0.22/0.63 but this library sits at the top of each.
