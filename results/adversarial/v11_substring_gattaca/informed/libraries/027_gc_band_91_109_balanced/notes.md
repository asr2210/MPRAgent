# E27: gc_band_91_109_balanced

## Design
GC band [91, 109] + per-position A↔T / C↔G balance. Slightly tighter
than E15's [90, 110].

## Result
- eval_01: 0.8765 (vs E15 0.8777, -0.001)
- mean_r: 0.8573 (vs E15 0.8591, -0.002)
- Within noise of E15

## Interpretation
[91, 109] is essentially the same as [90, 110]. Peak is a flat plateau.

## Next
E28: Try uniform-on-band GC. E15 uses rejection from Binomial(200,0.5),
which gives a truncated Binomial — more mass at GC=100. E28 forces
uniform on [90,110] (equal mass at each GC value). Tests if eval prefers
flat vs peaked GC on the band.
