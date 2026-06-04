"""
E25: gc_band_no_homopolymer

E15's GC band [90, 110] + reject sequences with any 10-bp homopolymer
run (10 consecutive identical bases). Tests if rare local-bias outliers
are noise.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
GC_LOW, GC_HIGH = 90, 110
HOMO_LEN = 10

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])

def has_homopolymer(batch, k):
    """Detect if any sequence in batch has run of k+ identical bases."""
    # diff: 0 if adjacent equal, else nonzero
    # consecutive equals: run of (k-1) zeros in diff means run of k identical bases
    n, l = batch.shape
    has_run = np.zeros(n, dtype=bool)
    # cumulative count of consecutive equals starting from left
    # simpler: convolve a "1 if equal to previous" indicator over window of (k-1)
    is_same = (batch[:, 1:] == batch[:, :-1]).astype(np.int8)
    # rolling sum of window k-1
    window = k - 1
    if window <= 0:
        return np.ones(n, dtype=bool)
    csum = np.cumsum(is_same, axis=1)
    # rolling sum of length `window` over is_same
    # sum from i to i+window-1 = csum[i+window-1] - csum[i-1] (with csum[-1] = 0)
    # If any rolling sum == window, that's a run of k+ identical
    if is_same.shape[1] < window:
        return has_run
    rolling = csum[:, window-1:] - np.concatenate([np.zeros((n, 1), dtype=np.int64), csum[:, :-window]], axis=1)
    has_run = (rolling >= window).any(axis=1)
    return has_run

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    BATCH = 100_000
    selected = []
    total_drawn = 0
    n_after_gc = 0
    n_after_homo = 0
    while sum(len(s) for s in selected) < N_SEQS:
        batch = rng.integers(0, 4, size=(BATCH, SEQ_LEN), dtype=np.int8)
        total_drawn += BATCH
        is_gc = (batch == 1) | (batch == 2)
        gc = is_gc.sum(axis=1)
        gc_mask = (gc >= GC_LOW) & (gc <= GC_HIGH)
        n_after_gc += int(gc_mask.sum())
        # filter by homopolymer
        kept_post_gc = batch[gc_mask]
        homo = has_homopolymer(kept_post_gc, HOMO_LEN)
        no_homo = ~homo
        n_after_homo += int(no_homo.sum())
        selected.append(kept_post_gc[no_homo])
        print(f"  drew {total_drawn}, gc_pass={n_after_gc}, homo_pass={n_after_homo}, kept={sum(len(s) for s in selected)}")
    seqs = np.vstack(selected)[:N_SEQS]
    is_gc = (seqs == 1) | (seqs == 2)
    gc = is_gc.sum(axis=1)
    print(f"GC: min={gc.min()}, max={gc.max()}, mean={gc.mean():.2f}, sd={gc.std():.2f}")
    # Verify no homopolymer
    homo = has_homopolymer(seqs, HOMO_LEN)
    print(f"Has homopolymer ≥{HOMO_LEN}: {homo.sum()} (should be 0)")
    rng.shuffle(seqs)
    with open(out_path, "w") as f:
        for row in seqs:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
