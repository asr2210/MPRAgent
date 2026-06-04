"""
E20: gc_gaussian_weighted

IID random with importance acceptance: accept each sequence with
probability exp(-(GC-100)^2 / (2*sigma^2)), sigma=5. Soft taper instead
of hard cutoff. Tests if the eval's GC distribution is a smooth Gaussian
vs E15's hard-restricted uniform.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
SIGMA = 5.0

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
        gc = is_gc.sum(axis=1).astype(np.float64)
        prob = np.exp(-(gc - 100.0) ** 2 / (2 * SIGMA * SIGMA))
        u = rng.random(BATCH)
        mask = u < prob
        selected.append(batch[mask])
        kept = sum(len(s) for s in selected)
        print(f"  drew {total_drawn}, kept {kept}, accept_rate={kept/total_drawn:.3f}")
    seqs = np.vstack(selected)[:N_SEQS]
    is_gc = (seqs == 1) | (seqs == 2)
    gc = is_gc.sum(axis=1)
    print(f"GC range: {gc.min()}-{gc.max()}, mean={gc.mean():.2f}, sd={gc.std():.2f}")
    # histogram
    hist, edges = np.histogram(gc, bins=np.arange(70, 131))
    print("GC histogram (peaks):")
    for i in range(len(hist)):
        if hist[i] > 0:
            print(f"  {edges[i]}: {hist[i]}")
    rng.shuffle(seqs)
    with open(out_path, "w") as f:
        for row in seqs:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
