"""
004_dhs_specific
================
50k 200bp DHS summit-centered sequences, weighted to favor cell-type-specific
elements (low numsamples). Sampling weight = 1/numsamples (cell-type-specific
elements get oversampled). Approximates the topic-weighted sampling that the
instructions.md baseline `dhs_topic` claims to use.

Hypothesis: cell-type-specific elements may carry stronger / more distinct
motifs than tissue-invariant ones, providing more learnable signal even if
eval doesn't directly reward DHS. If eval_01 score is similar to dhs_random
(0.041), specificity doesn't matter here.
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
            chrom = parts[0]
            if chrom not in CANONICAL:
                continue
            chroms.append(chrom)
            summits.append(int(parts[6]))
            nsamples.append(int(parts[5]))
    return chroms, np.array(summits), np.array(nsamples, dtype=np.int32)


def main():
    rng = np.random.default_rng(SEED)
    print("Loading DHS...", file=sys.stderr)
    chroms, summits, nsamples = load_dhs()
    print(f"  {len(chroms)} records", file=sys.stderr)
    # Weight inversely with number of biosamples (favor cell-type-specific).
    w = 1.0 / np.sqrt(nsamples.astype(np.float64))
    w /= w.sum()
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    out_path = Path(__file__).parent / "sequences_0.txt"
    seqs = []
    # Oversample to handle Ns.
    pool = rng.choice(len(chroms), size=int(N_SEQS * 1.3), replace=False, p=w)
    pi = 0
    while len(seqs) < N_SEQS and pi < len(pool):
        i = pool[pi]
        pi += 1
        chrom, summit = chroms[i], int(summits[i])
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
        if len(seqs) % 10_000 == 0:
            print(f"  {len(seqs)}/{N_SEQS}", file=sys.stderr)

    assert len(seqs) == N_SEQS
    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {len(seqs)}", file=sys.stderr)


if __name__ == "__main__":
    main()
