# E17: gc_band_tighter

## Design
IID random rejection to per-seq GC ∈ [93, 107] (~70% acceptance,
SD=4.0). Narrower than E15.

## Result
- eval_01: 0.8683 (vs E15 0.8777, -0.009; vs E2 0.8565, +0.012)
- mean_r: 0.8474 (vs E15 0.8591, -0.012)
- SK-N-SH eval_07: 0.5832 (vs E15 0.7183, -0.135) — starting to crash
- SK-N-SH eval_13: 0.7203 (vs E15 0.8217, -0.10)

## Interpretation — peak found, descent confirmed
SK-N-SH starts crashing again as we tighten further. The GC band sweep:
- [99, 101]  → 0.7929 mean (E14, far too tight)
- [93, 107]  → 0.8474 mean (E17, tighter — SKNSH degrading)
- [90, 110]  → 0.8591 mean (E15, PEAK)
- [85, 115]  → 0.8522 mean (E16, wider — slightly past peak)
- no filter  → 0.8408 mean (E2 random)

**Peak around [90, 110]**. Could fine-tune to [88, 112] but the gain is
small. SK-N-SH is the sensitive component: too-tight = crashes; too-loose
= mild degradation.

## Theory v9 holds
The eval distribution is IID Uniform with soft GC filtering, peak
matching achieved at GC ∈ [90, 110]. SK-N-SH eval_07/13 are the
sensitive detectors — they pick up both too-tight and too-loose
filtering.

## Next
E18: bisect [88, 112] (between E15 and E16). If matches or slightly
beats E15, optimum is right there. Then move on to layered improvements.
