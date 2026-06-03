"""
003_encode_ccre
===============
50k 200bp windows centered on midpoints of randomly-sampled ENCODE v4
cCREs (~2.35M total). Different annotation than DHS — typed by enhancer
(dELS/pELS), promoter (PLS), CTCF (CA-CTCF), etc.

Hypothesis: if the eval data was built around ENCODE cCRE-like elements
(as Agarwal et al. lentiMPRA uses), this should beat DHS-random. If it
also scores ~0.04, the issue isn't about which annotation we use; it's
that any standard genomic library can't hit the eval's hidden structure.
"""
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
CCRE_PATH = DATA_DIR / "ccre_v4.bed"
HG38_PATH = DATA_DIR / "hg38.fa"

CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def load_ccres():
    recs = []
    with open(CCRE_PATH) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            if chrom not in CANONICAL:
                continue
            start, end = int(parts[1]), int(parts[2])
            mid = (start + end) // 2
            recs.append((chrom, mid))
    return recs


def main():
    rng = np.random.default_rng(SEED)
    print("Loading cCREs...", file=sys.stderr)
    recs = load_ccres()
    print(f"  {len(recs)} canonical cCREs", file=sys.stderr)
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    out_path = Path(__file__).parent / "sequences_0.txt"
    seqs = []
    perm = rng.permutation(len(recs))
    pi = 0
    while len(seqs) < N_SEQS and pi < len(perm):
        chrom, mid = recs[perm[pi]]
        pi += 1
        L = len(fa[chrom])
        start = max(0, mid - HALF)
        end = start + SEQ_LEN
        if end > L:
            end = L
            start = end - SEQ_LEN
        if start < 0:
            continue
        s = str(fa[chrom][start:end])
        if len(s) != SEQ_LEN or set(s) - set("ACGT"):
            continue
        seqs.append(s)
        if len(seqs) % 10_000 == 0:
            print(f"  {len(seqs)}/{N_SEQS}", file=sys.stderr)

    assert len(seqs) == N_SEQS
    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {len(seqs)} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
