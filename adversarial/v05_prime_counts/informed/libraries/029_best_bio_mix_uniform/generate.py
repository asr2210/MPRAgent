"""
029_best_bio_mix_uniform
========================
33/67 mix where bio half is 50/50 mix of two best bio sources:
  - 8333 dhs_specific (1/sqrt weighting)
  - 8333 cCRE_high_confidence (PLS+pELS+CA-TF+CA-H3K4me3+TF)
  - 33333 pure uniform random
Total 50k.

Tests if combining two top bio sources in the 33% bio slot beats
single-source DHS at 33% (best mean_r so far: 0.0488 with Exp 019).
"""
import gzip
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
N_DHS = 8_333
N_CCRE = 8_334
N_UNIFORM = 50_000 - N_DHS - N_CCRE
N_SEQS = 50_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
DHS_PATH = DATA_DIR / "dhs_index.txt.gz"
CCRE_PATH = DATA_DIR / "ccre_v4.bed"
HG38_PATH = DATA_DIR / "hg38.fa"
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
HIGH_CONF = {"PLS", "pELS", "CA-TF", "CA-H3K4me3", "TF"}


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


def load_ccre_high():
    recs = []
    with open(CCRE_PATH) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom, cls = parts[0], parts[5]
            if chrom not in CANONICAL or cls not in HIGH_CONF:
                continue
            mid = (int(parts[1]) + int(parts[2])) // 2
            recs.append((chrom, mid))
    return recs


def extract(fa, chrom, center):
    L = len(fa[chrom])
    start = max(0, center - HALF)
    end = start + SEQ_LEN
    if end > L:
        end = L
        start = end - SEQ_LEN
    if start < 0:
        return None
    s = str(fa[chrom][start:end])
    if len(s) != SEQ_LEN or set(s) - set("ACGT"):
        return None
    return s


def main():
    rng = np.random.default_rng(SEED)
    chroms, summits, nsamples = load_dhs()
    w = 1.0 / np.sqrt(nsamples.astype(np.float64))
    w /= w.sum()
    ccre = load_ccre_high()
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    dhs_seqs = []
    pool = rng.choice(len(chroms), size=int(N_DHS * 1.5), replace=False, p=w)
    for i in pool:
        if len(dhs_seqs) >= N_DHS:
            break
        s = extract(fa, chroms[i], int(summits[i]))
        if s:
            dhs_seqs.append(s)
    assert len(dhs_seqs) == N_DHS

    ccre_seqs = []
    perm = rng.permutation(len(ccre))
    for i in perm:
        if len(ccre_seqs) >= N_CCRE:
            break
        chrom, mid = ccre[i]
        s = extract(fa, chrom, mid)
        if s:
            ccre_seqs.append(s)
    assert len(ccre_seqs) == N_CCRE

    bases = np.array(list("ACGT"))
    idx = rng.integers(0, 4, size=(N_UNIFORM, SEQ_LEN))
    uniform_seqs = ["".join(row) for row in bases[idx]]

    all_seqs = dhs_seqs + ccre_seqs + uniform_seqs
    rng.shuffle(all_seqs)
    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in all_seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
