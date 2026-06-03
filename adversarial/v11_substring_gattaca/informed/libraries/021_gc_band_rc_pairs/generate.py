"""
E21: gc_band_rc_pairs

25k random IID with GC ∈ [90, 110], plus their 25k reverse complements.
RC preserves GC (C↔G, A↔T, so GC count unchanged). Tests if RC pairing
adds value when layered on the GC band optimum.
"""
import numpy as np
import os

SEED = 0
N_HALF = 25_000
SEQ_LEN = 200
GC_LOW, GC_HIGH = 90, 110

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])
# RC: A↔T (0↔3), C↔G (1↔2)
RC_MAP = np.array([3, 2, 1, 0], dtype=np.int8)

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    BATCH = 100_000
    selected = []
    total_drawn = 0
    while sum(len(s) for s in selected) < N_HALF:
        batch = rng.integers(0, 4, size=(BATCH, SEQ_LEN), dtype=np.int8)
        total_drawn += BATCH
        is_gc = (batch == 1) | (batch == 2)
        gc = is_gc.sum(axis=1)
        mask = (gc >= GC_LOW) & (gc <= GC_HIGH)
        selected.append(batch[mask])
    fwd = np.vstack(selected)[:N_HALF]
    # compute RC: reverse + complement
    rc = RC_MAP[fwd[:, ::-1]]
    all_seqs = np.vstack([fwd, rc])
    print(f"Forward: {fwd.shape}, RC: {rc.shape}, total: {all_seqs.shape}")
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
