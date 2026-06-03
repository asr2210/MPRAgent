"""
002_dhs_random_summit
=====================
50,000 200bp sequences extracted from hg38, each centered on a randomly-
sampled DHS summit (Meuleman 2020, 3.59M index elements).

Hypothesis (post Exp 001 surprise): the cliff between synthetic (~0.04 in
this env) and DHS-derived libraries should be huge. Real genomic context is
necessary. This experiment establishes the dhs_random baseline (instructions.md
table claims dhs_random eval_01 = 0.7089).

Output: sequences_0.txt (50000 lines, each 200 chars from {A,C,G,T}).
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

CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def load_dhs_records():
    """Returns list of (chrom, summit) for canonical-chrom DHS sites."""
    recs = []
    with gzip.open(DHS_PATH, "rt") as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            if chrom not in CANONICAL_CHROMS:
                continue
            summit = int(parts[6])
            recs.append((chrom, summit))
    return recs


def main():
    rng = np.random.default_rng(SEED)
    print("Loading DHS index...", file=sys.stderr)
    recs = load_dhs_records()
    print(f"  {len(recs)} canonical-chrom DHS records", file=sys.stderr)

    print("Opening hg38...", file=sys.stderr)
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    out_path = Path(__file__).parent / "sequences_0.txt"
    seqs = []
    # Oversample to account for N-containing windows we'll discard.
    target = N_SEQS
    pool_idx = rng.permutation(len(recs))
    pi = 0
    while len(seqs) < target and pi < len(pool_idx):
        i = pool_idx[pi]
        pi += 1
        chrom, summit = recs[i]
        L = len(fa[chrom])
        start = max(0, summit - HALF)
        end = start + SEQ_LEN
        if end > L:
            end = L
            start = end - SEQ_LEN
        if start < 0:
            continue
        s = str(fa[chrom][start:end])
        if len(s) != SEQ_LEN:
            continue
        if set(s) - set("ACGT"):
            continue
        seqs.append(s)
        if len(seqs) % 10_000 == 0:
            print(f"  {len(seqs)}/{target}", file=sys.stderr)

    assert len(seqs) == N_SEQS, f"only collected {len(seqs)} valid windows"

    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {len(seqs)} sequences to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
