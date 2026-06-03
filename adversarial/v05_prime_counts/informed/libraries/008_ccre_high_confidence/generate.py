"""
008_ccre_high_confidence
========================
50k 200bp windows from high-confidence regulatory cCRE classes only:
  - PLS (promoter-like, 47k)
  - pELS (proximal enhancer-like, 249k)
  - CA-TF (chromatin accessible + TF binding, 26k)
  - CA-H3K4me3 (chromatin accessible + active promoter mark, 79k)
  - TF (TF binding, 105k)
~506k regulatory elements (vs 2.35M total cCREs). Excludes the more
ambiguous distal-only (dELS) and CTCF-only elements.

Hypothesis: if eval rewards stronger/cleaner regulatory signal, restricting
to high-confidence regulatory cCREs should beat full cCRE (0.047). If eval
is insensitive to regulatory cleanliness, this won't help.
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
HIGH_CONF_CLASSES = {"PLS", "pELS", "CA-TF", "CA-H3K4me3", "TF"}


def main():
    rng = np.random.default_rng(SEED)
    print("Loading high-confidence cCREs...", file=sys.stderr)
    recs = []
    with open(CCRE_PATH) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom, cls = parts[0], parts[5]
            if chrom not in CANONICAL or cls not in HIGH_CONF_CLASSES:
                continue
            mid = (int(parts[1]) + int(parts[2])) // 2
            recs.append((chrom, mid))
    print(f"  {len(recs)} high-confidence regulatory cCREs", file=sys.stderr)
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
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
