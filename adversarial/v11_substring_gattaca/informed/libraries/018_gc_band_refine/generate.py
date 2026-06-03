"""
E18: gc_band_refine

IID random rejection to per-seq GC ∈ [88, 112] (between E15 [90,110]
and E16 [85,115]). Bisects to refine the optimum.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
GC_LOW, GC_HIGH = 88, 112

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
        gc_counts = is_gc.sum(axis=1)
        mask = (gc_counts >= GC_LOW) & (gc_counts <= GC_HIGH)
        selected.append(batch[mask])
    seqs = np.vstack(selected)[:N_SEQS]
    print(f"Got {len(seqs)} from {total_drawn} drawn.")
    is_gc = (seqs == 1) | (seqs == 2)
    gc = is_gc.sum(axis=1)
    print(f"GC count: min={gc.min()}, max={gc.max()}, mean={gc.mean():.3f}, sd={gc.std():.3f}")
    rng.shuffle(seqs)
    with open(out_path, "w") as f:
        for row in seqs:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
