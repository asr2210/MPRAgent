"""
E19: gc_and_perbase_band

IID random rejection to:
- per-seq GC ∈ [90, 110]  (E15's optimum)
- each base count (A, T, C, G) ∈ [40, 60]

Tests if eval has per-base filtering beyond GC.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
GC_LOW, GC_HIGH = 90, 110
BASE_LOW, BASE_HIGH = 40, 60

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    BATCH = 200_000
    selected = []
    total_drawn = 0
    n_after_gc = 0
    while sum(len(s) for s in selected) < N_SEQS:
        batch = rng.integers(0, 4, size=(BATCH, SEQ_LEN), dtype=np.int8)
        total_drawn += BATCH
        # per-position is 0=A, 1=C, 2=G, 3=T
        counts = np.zeros((BATCH, 4), dtype=np.int32)
        for b in range(4):
            counts[:, b] = (batch == b).sum(axis=1)
        gc = counts[:, 1] + counts[:, 2]
        gc_mask = (gc >= GC_LOW) & (gc <= GC_HIGH)
        base_mask = np.all((counts >= BASE_LOW) & (counts <= BASE_HIGH), axis=1)
        n_after_gc += int(gc_mask.sum())
        combined = gc_mask & base_mask
        selected.append(batch[combined])
        kept = sum(len(s) for s in selected)
        print(f"  drew {total_drawn}, kept {kept}, gc_pass={n_after_gc}, combined_pass_rate={kept/total_drawn:.3f}")
    seqs = np.vstack(selected)[:N_SEQS]
    # verify
    counts = np.zeros((len(seqs), 4), dtype=np.int32)
    for b in range(4):
        counts[:, b] = (seqs == b).sum(axis=1)
    print(f"A range: {counts[:,0].min()}-{counts[:,0].max()}")
    print(f"C range: {counts[:,1].min()}-{counts[:,1].max()}")
    print(f"G range: {counts[:,2].min()}-{counts[:,2].max()}")
    print(f"T range: {counts[:,3].min()}-{counts[:,3].max()}")
    gc = counts[:, 1] + counts[:, 2]
    print(f"GC range: {gc.min()}-{gc.max()}, mean={gc.mean():.2f}, sd={gc.std():.2f}")
    rng.shuffle(seqs)
    with open(out_path, "w") as f:
        for row in seqs:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
