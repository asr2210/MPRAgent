"""
009_dhs_specific_rc
===================
Same as 004_dhs_specific but every sequence is REVERSE-COMPLEMENTED.
Tests whether the prepare.py pipeline is strand-symmetric. If yes, score
should match 004 (~0.049). If sensitive to orientation, score will differ.

Important because Agarwal et al. (2025) explicitly showed promoters have
strand-orientation bias in lentiMPRA — the eval may also reflect that.
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
COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


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
    w = 1.0 / np.sqrt(nsamples.astype(np.float64))
    w /= w.sum()
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    out_path = Path(__file__).parent / "sequences_0.txt"
    seqs = []
    pool = rng.choice(len(chroms), size=int(N_SEQS * 1.3), replace=False, p=w)
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
        seqs.append(rc(s))  # <-- reverse complement!
    assert len(seqs) == N_SEQS
    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS} RC'd sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
