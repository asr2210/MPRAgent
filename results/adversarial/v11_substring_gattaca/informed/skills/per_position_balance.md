# Per-position A↔T / C↔G swap balance

## What it is
Stack on top of a GC-band library: for each of the 200 positions,
balance the count of A vs T (swap excess to underrepresented) and
the count of C vs G. Both swap types preserve per-sequence GC count
exactly, so they layer cleanly on any GC-filtered pool.

## When to use
Apply AFTER any GC-band filtering. Adds ~+0.001 mean (within seed
noise, but consistent across two test cases). Cheap, zero downside.

## Recipe
```python
import numpy as np

def balance_position(seqs, p, b1, b2, rng):
    """Swap excess of base b1 at position p to base b2 (or vice versa).
    seqs: (N, 200) int8 with A=0, C=1, G=2, T=3."""
    col = seqs[:, p]
    idx_b1 = np.where(col == b1)[0]
    idx_b2 = np.where(col == b2)[0]
    n1, n2 = len(idx_b1), len(idx_b2)
    target = (n1 + n2) // 2
    if n1 > target:
        swap_idx = rng.choice(idx_b1, size=n1 - target, replace=False)
        seqs[swap_idx, p] = b2
    elif n2 > target:
        swap_idx = rng.choice(idx_b2, size=n2 - target, replace=False)
        seqs[swap_idx, p] = b1

# Usage after building `seqs` from GC-band rejection:
for p in range(SEQ_LEN):
    balance_position(seqs, p, 0, 3, rng)  # A<->T
    balance_position(seqs, p, 1, 2, rng)  # C<->G
```

## Effect attribution
- Per-pos balance alone (E26 vs random): +0.005 mean
- GC band alone (E15 vs random): +0.018 mean
- GC band + per-pos balance (E30 vs E29): +0.0002 mean (within noise)
- GC band + per-pos balance vs random (E30 vs E2): +0.022 mean total

The two levers stack but are nearly additive only when applied to
different baselines. On top of GC-band, balance saturates fast — the
GC filter already brings per-position uniformity to within 0.8%.

## Why it preserves per-seq GC
- A↔T swap: changes A→T or T→A, neither is GC
- C↔G swap: changes C→G or G→C, both are GC
So total GC per row is invariant under both ops.

## Per-position output check
After balancing, per-position counts hit within 1 of perfect:
```
pos 0: A=12555, C=12445, G=12444, T=12556   # A↔T balanced, C↔G balanced
```
A+T may not equal C+G (random per-position GC variance is preserved).
