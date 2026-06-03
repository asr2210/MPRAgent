"""
E24: gc_band_perpos_balanced

Start with E15's GC-band [90, 110] random pool. For each position p,
balance the per-base counts to exactly 12500 each by swapping bases.
A↔T swaps and C↔G swaps preserve per-seq GC.

Tests if removing per-position sampling noise (~0.78% deviation per
position) helps.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
GC_LOW, GC_HIGH = 90, 110
TARGET_PER_BASE = N_SEQS // 4  # 12500

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])
# indices: A=0, C=1, G=2, T=3

def gen_gc_band_pool(n):
    BATCH = 100_000
    selected = []
    while sum(len(s) for s in selected) < n:
        batch = rng.integers(0, 4, size=(BATCH, SEQ_LEN), dtype=np.int8)
        is_gc = (batch == 1) | (batch == 2)
        gc = is_gc.sum(axis=1)
        mask = (gc >= GC_LOW) & (gc <= GC_HIGH)
        selected.append(batch[mask])
    return np.vstack(selected)[:n]

def balance_position(seqs, p, b1, b2):
    """Balance positions where base is b1 vs b2 by swapping. Preserves GC."""
    col = seqs[:, p]
    idx_b1 = np.where(col == b1)[0]
    idx_b2 = np.where(col == b2)[0]
    n1, n2 = len(idx_b1), len(idx_b2)
    target = (n1 + n2) // 2
    # If we want target of each, but total may be odd, the extra goes to b1 by default
    # We need to swap (n1 - target) elements from b1 to b2 if n1 > target
    if n1 > target:
        n_swap = n1 - target
        swap_idx = rng.choice(idx_b1, size=n_swap, replace=False)
        seqs[swap_idx, p] = b2
    elif n2 > target:
        n_swap = n2 - target
        swap_idx = rng.choice(idx_b2, size=n_swap, replace=False)
        seqs[swap_idx, p] = b1

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    print("Generating GC-band pool...")
    seqs = gen_gc_band_pool(N_SEQS)
    print(f"Pool: {seqs.shape}")
    is_gc = (seqs == 1) | (seqs == 2)
    gc_before = is_gc.sum(axis=1)
    print(f"GC before balance: min={gc_before.min()}, max={gc_before.max()}, mean={gc_before.mean():.2f}")
    print("Balancing per-position via A<->T and C<->G swaps...")
    for p in range(SEQ_LEN):
        balance_position(seqs, p, 0, 3)  # A vs T
        balance_position(seqs, p, 1, 2)  # C vs G
    # Verify
    is_gc = (seqs == 1) | (seqs == 2)
    gc_after = is_gc.sum(axis=1)
    # Per-position check
    print(f"GC after balance: min={gc_after.min()}, max={gc_after.max()}, mean={gc_after.mean():.2f}")
    print(f"Per-seq GC preserved exactly: {np.all(gc_before == gc_after)}")
    # Per-position base counts at first 5 positions
    for p in [0, 1, 2, 100, 199]:
        col = seqs[:, p]
        counts = [int((col == b).sum()) for b in range(4)]
        print(f"  pos {p}: A={counts[0]}, C={counts[1]}, G={counts[2]}, T={counts[3]}")
    rng.shuffle(seqs)
    with open(out_path, "w") as f:
        for row in seqs:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
