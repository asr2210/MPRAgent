"""
E28: gc_uniform_on_band

Force per-seq GC distribution UNIFORM on [90, 110] — exactly the same
number of sequences at each integer GC value (2381 each across 21 GC
values * 21 = 50001, take 50000).

Vs E15 which gives truncated-Binomial distribution (peak at GC=100,
tapering toward 90/110).

Tests if eval prefers flat-on-band or peaked-on-band.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
GC_VALS = list(range(90, 111))  # 21 values
PER_BIN = N_SEQS // len(GC_VALS) + 1  # 2381

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    # For each GC value, generate PER_BIN sequences via rejection
    all_seqs = []
    BATCH = 50_000
    for target_gc in GC_VALS:
        collected = 0
        bin_seqs = []
        while collected < PER_BIN:
            batch = rng.integers(0, 4, size=(BATCH, SEQ_LEN), dtype=np.int8)
            is_gc = (batch == 1) | (batch == 2)
            gc = is_gc.sum(axis=1)
            mask = (gc == target_gc)
            bin_seqs.append(batch[mask])
            collected = sum(len(s) for s in bin_seqs)
        bin_arr = np.vstack(bin_seqs)[:PER_BIN]
        all_seqs.append(bin_arr)
        print(f"  GC={target_gc}: collected {len(bin_arr)} sequences")
    seqs = np.vstack(all_seqs)[:N_SEQS]
    is_gc = (seqs == 1) | (seqs == 2)
    gc = is_gc.sum(axis=1)
    print(f"Total: {len(seqs)} sequences. GC: min={gc.min()}, max={gc.max()}, mean={gc.mean():.2f}, sd={gc.std():.2f}")
    # Histogram
    counts = np.bincount(gc, minlength=130)
    print("GC histogram:")
    for g in GC_VALS:
        print(f"  {g}: {counts[g]}")
    rng.shuffle(seqs)
    with open(out_path, "w") as f:
        for row in seqs:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
