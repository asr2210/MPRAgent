"""
010_dhs_high_signal
===================
50k 200bp DHS summit-centered windows from the TOP 10% of DHS sites by
mean_signal (strongest accessibility = most active regulatory elements).
~360k candidates.

Hypothesis: highly-active regulatory elements provide stronger / cleaner
activity signal in the MPRA, helping the model learn better. If the
ceiling at ~0.05 is due to weak signal in random DHS, this should break it.
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
    chroms, summits, signals = [], [], []
    with gzip.open(DHS_PATH, "rt") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] not in CANONICAL:
                continue
            chroms.append(parts[0])
            summits.append(int(parts[6]))
            signals.append(float(parts[4]))
    return chroms, np.array(summits), np.array(signals)


def main():
    rng = np.random.default_rng(SEED)
    chroms, summits, signals = load_dhs()
    print(f"Loaded {len(chroms)} DHS sites", file=sys.stderr)
    # Top 10% by mean_signal
    threshold = np.quantile(signals, 0.90)
    keep_idx = np.where(signals >= threshold)[0]
    print(f"Top 10% threshold={threshold:.3f}, n={len(keep_idx)}", file=sys.stderr)
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    out_path = Path(__file__).parent / "sequences_0.txt"
    seqs = []
    perm = rng.permutation(keep_idx)
    pi = 0
    while len(seqs) < N_SEQS and pi < len(perm):
        i = perm[pi]
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
    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
