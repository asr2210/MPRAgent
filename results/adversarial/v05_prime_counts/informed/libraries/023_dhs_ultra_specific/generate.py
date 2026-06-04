"""
023_dhs_ultra_specific
======================
50k DHS summit-centered windows weighted by 1/numsamples (vs Exp 004's
1/sqrt(numsamples)). MORE aggressive specificity bias — sites present
in only 1 sample weighted 100x more than sites in 100 samples.

Hypothesis: if eval_01 rewards tissue-specific elements, stricter
specificity weighting should beat dhs_specific (0.0487).
"""
import gzip
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
DHS_PATH = DATA_DIR / "dhs_index.txt.gz"
HG38_PATH = DATA_DIR / "hg38.fa"
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def load_dhs():
    chroms, summits, nsamples = [], [], []
    with gzip.open(DHS_PATH, "rt") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] not in CANONICAL:
                continue
            chroms.append(parts[0])
            summits.append(int(parts[6]))
            nsamples.append(int(parts[5]))
    return chroms, np.array(summits), np.array(nsamples, dtype=np.int32)


def main():
    rng = np.random.default_rng(SEED)
    chroms, summits, nsamples = load_dhs()
    w = 1.0 / nsamples.astype(np.float64)
    w /= w.sum()
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    pool = rng.choice(len(chroms), size=int(N_SEQS * 1.3), replace=False, p=w)
    seqs = []
    pi = 0
    while len(seqs) < N_SEQS and pi < len(pool):
        i = pool[pi]
        pi += 1
        chrom = chroms[i]
        summit = int(summits[i])
        L = len(fa[chrom])
        start = max(0, summit - HALF)
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
    assert len(seqs) == N_SEQS, len(seqs)
    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
