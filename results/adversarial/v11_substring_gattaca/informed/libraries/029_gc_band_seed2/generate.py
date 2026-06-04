"""
E29: gc_band_seed2

Third reproducibility check on E15's GC band [90, 110] design.
E15 (seed 0): 0.8777 eval_01 / 0.8591 mean
E22 (seed 1): 0.8755 eval_01 / 0.8578 mean
E29 (seed 2): expect 0.875 ± 0.006.
"""
import numpy as np
import os

SEED = 2
N_SEQS = 50_000
SEQ_LEN = 200
GC_LOW, GC_HIGH = 90, 110

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    BATCH = 100_000
    selected = []
    while sum(len(s) for s in selected) < N_SEQS:
        batch = rng.integers(0, 4, size=(BATCH, SEQ_LEN), dtype=np.int8)
        is_gc = (batch == 1) | (batch == 2)
        gc = is_gc.sum(axis=1)
        mask = (gc >= GC_LOW) & (gc <= GC_HIGH)
        selected.append(batch[mask])
    seqs = np.vstack(selected)[:N_SEQS]
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
