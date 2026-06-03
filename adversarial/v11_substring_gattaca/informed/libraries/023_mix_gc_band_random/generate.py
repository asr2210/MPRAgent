"""
E23: mix_gc_band_random

25k IID random with GC ∈ [90, 110] + 25k pure IID random (no filter).
Tests if SK-N-SH eval_07 recovers while preserving most of E15's gain.
"""
import numpy as np
import os

SEED = 0
N_HALF = 25_000
SEQ_LEN = 200
GC_LOW, GC_HIGH = 90, 110

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    BATCH = 100_000
    selected = []
    while sum(len(s) for s in selected) < N_HALF:
        batch = rng.integers(0, 4, size=(BATCH, SEQ_LEN), dtype=np.int8)
        is_gc = (batch == 1) | (batch == 2)
        gc = is_gc.sum(axis=1)
        mask = (gc >= GC_LOW) & (gc <= GC_HIGH)
        selected.append(batch[mask])
    filtered = np.vstack(selected)[:N_HALF]
    random = rng.integers(0, 4, size=(N_HALF, SEQ_LEN), dtype=np.int8)
    all_seqs = np.vstack([filtered, random])
    print(f"Filtered: {filtered.shape}, Random: {random.shape}, total: {all_seqs.shape}")
    is_gc = (all_seqs == 1) | (all_seqs == 2)
    gc = is_gc.sum(axis=1)
    print(f"GC: min={gc.min()}, max={gc.max()}, mean={gc.mean():.2f}, sd={gc.std():.2f}")
    rng.shuffle(all_seqs)
    with open(out_path, "w") as f:
        for row in all_seqs:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {len(all_seqs)} sequences to {out_path}")

if __name__ == "__main__":
    main()
