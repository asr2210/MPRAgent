"""
024_tss_proximal
================
50k 200bp windows centered on GENCODE basic transcript TSS positions.
Promoter regions are motif-dense and known drivers of MPRA activity.
Filter to protein-coding genes only for highest signal.
"""
import gzip
import re
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
TSS_PATH = DATA_DIR / "gencode_basic_tss.gtf.gz"
HG38_PATH = DATA_DIR / "hg38.fa"
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def load_tss():
    recs = []
    seen = set()
    pat = re.compile(r'transcript_type "([^"]+)"')
    with gzip.open(TSS_PATH, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            chrom = parts[0]
            if chrom not in CANONICAL:
                continue
            attr = parts[8]
            m = pat.search(attr)
            if not m or m.group(1) != "protein_coding":
                continue
            strand = parts[6]
            tss = int(parts[3]) - 1 if strand == "+" else int(parts[4]) - 1
            key = (chrom, tss, strand)
            if key in seen:
                continue
            seen.add(key)
            recs.append((chrom, tss))
    return recs


def main():
    rng = np.random.default_rng(SEED)
    print("Loading TSS...", file=sys.stderr)
    recs = load_tss()
    print(f"  {len(recs)} unique protein-coding TSS", file=sys.stderr)
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    perm = rng.permutation(len(recs))
    seqs = []
    pi = 0
    while len(seqs) < N_SEQS and pi < len(perm):
        chrom, tss = recs[perm[pi]]
        pi += 1
        L = len(fa[chrom])
        start = max(0, tss - HALF)
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

    # If we ran out, recycle with replacement
    if len(seqs) < N_SEQS:
        print(f"  only got {len(seqs)} unique; filling by random pick", file=sys.stderr)
        while len(seqs) < N_SEQS:
            chrom, tss = recs[rng.integers(0, len(recs))]
            L = len(fa[chrom])
            start = max(0, tss - HALF)
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
    assert len(seqs) == N_SEQS
    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
