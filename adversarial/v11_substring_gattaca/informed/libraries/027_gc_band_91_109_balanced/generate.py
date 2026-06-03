"""
E27: gc_band_91_109_balanced

GC band [91, 109] (slightly tighter than E15) + per-position A↔T / C↔G
balance to remove residual per-position noise. Fine-grained exploration
of the [90, 110] / [89, 111] / [91, 109] neighborhood.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
GC_LOW, GC_HIGH = 91, 109

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])

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

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    seqs = gen_gc_band_pool(N_SEQS)
    print(f"Pool: {seqs.shape}")
    for p in range(SEQ_LEN):
        balance_position(seqs, p, 0, 3)
        balance_position(seqs, p, 1, 2)
    is_gc = (seqs == 1) | (seqs == 2)
    gc = is_gc.sum(axis=1)
    print(f"GC: min={gc.min()}, max={gc.max()}, mean={gc.mean():.2f}, sd={gc.std():.2f}")
    rng.shuffle(seqs)
    with open(out_path, "w") as f:
        for row in seqs:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
