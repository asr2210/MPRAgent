"""Experiment 005 — random hg38 genomic windows at GC ∈ [0.45, 0.55].

Sample 50k 200bp windows uniformly at random from autosomes (chr1-22, X, Y),
then filter for ACGT-only and GC content in [0.45, 0.55]. Mostly intergenic
and intronic — biologically "unselected" but composition-matched.

Hypothesis: tests whether the curation status of natural genomic sequence
matters at this scale. cCRE GC-matched scored 0.392; gc_50 random scored
0.397. Where does random natural genome land?
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)
from utils.seqlib import extract, write_sequences, SEQLEN, fasta

OUT = os.path.join(THIS, "sequences_0.txt")
SEED = 1
N = 50000
GC_MIN = 0.45
GC_MAX = 0.55


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    fa = fasta()
    # Use primary chromosomes only
    main_chroms = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]
    chroms = [c for c in main_chroms if c in fa]
    sizes = np.array([len(fa[c]) for c in chroms], dtype=np.float64)
    probs = sizes / sizes.sum()
    print(f"Using {len(chroms)} chromosomes, total {int(sizes.sum()):,} bp")

    seqs = []
    seen = set()
    tries = 0
    rejected_gc = 0
    rejected_n = 0
    batch = 5000
    while len(seqs) < N:
        ci = rng.choice(len(chroms), size=batch, p=probs)
        starts = rng.integers(0, sizes[ci].astype(np.int64) - SEQLEN, dtype=np.int64)
        for k in range(batch):
            tries += 1
            chrom = chroms[ci[k]]
            center = int(starts[k]) + SEQLEN // 2
            seq = extract(chrom, center)
            if seq is None:
                rejected_n += 1
                continue
            if seq in seen:
                continue
            g = sum(1 for c in seq if c in "GC") / SEQLEN
            if g < GC_MIN or g > GC_MAX:
                rejected_gc += 1
                continue
            seen.add(seq)
            seqs.append(seq)
            if len(seqs) >= N:
                break

    print(f"Tried {tries:,} extractions: N/edge rejects {rejected_n:,}, GC rejects {rejected_gc:,}, collected {len(seqs):,}")
    rng.shuffle(seqs)
    write_sequences(seqs, OUT)
    gcs = [sum(1 for c in s if c in "GC")/SEQLEN for s in seqs]
    print(f"Final GC: mean={np.mean(gcs):.3f}, std={np.std(gcs):.3f}")
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
