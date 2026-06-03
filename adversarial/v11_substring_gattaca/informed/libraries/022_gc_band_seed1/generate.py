"""
E22: gc_band_seed1

Replicates E15 with seed=1. Confirms reproducibility of the +0.02
gain from GC band filtering.
"""
import numpy as np
import os

SEED = 1  # different seed
N_SEQS = 50_000
SEQ_LEN = 200
GC_LOW, GC_HIGH = 90, 110

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    BATCH = 100_000
    selected = []
    total_drawn = 0
    while sum(len(s) for s in selected) < N_SEQS:
        batch = rng.integers(0, 4, size=(BATCH, SEQ_LEN), dtype=np.int8)
        total_drawn += BATCH
        is_gc = (batch == 1) | (batch == 2)
        gc = is_gc.sum(axis=1)
        mask = (gc >= GC_LOW) & (gc <= GC_HIGH)
        selected.append(batch[mask])
    seqs = np.vstack(selected)[:N_SEQS]
    print(f"Got {len(seqs)} from {total_drawn} drawn.")
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
