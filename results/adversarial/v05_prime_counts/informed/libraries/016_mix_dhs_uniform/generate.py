"""
016_mix_dhs_uniform
===================
25k dhs_specific + 25k pure uniform random. Combines:
  - dhs_specific (best on eval_01 0.049, best mean_r 0.045)
  - pure_uniform (best on eval_08 0.124, tied mean_r 0.045)

Hypothesis: 50/50 mix may break the mean_r ceiling because the bio half
satisfies eval_01-07/09-14 and the uniform half drives eval_08 sky-high.
"""
import gzip
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
HALF_N = 25_000
N_SEQS = HALF_N * 2
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
    w = 1.0 / np.sqrt(nsamples.astype(np.float64))
    w /= w.sum()
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    # 25k DHS specific
    dhs_seqs = []
    pool = rng.choice(len(chroms), size=int(HALF_N * 1.5), replace=False, p=w)
    for i in pool:
        if len(dhs_seqs) >= HALF_N:
            break
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
        dhs_seqs.append(s)
    assert len(dhs_seqs) == HALF_N, len(dhs_seqs)

    # 25k uniform random
    bases = np.array(list("ACGT"))
    idx = rng.integers(0, 4, size=(HALF_N, SEQ_LEN))
    uniform_seqs = ["".join(row) for row in bases[idx]]

    all_seqs = dhs_seqs + uniform_seqs
    rng.shuffle(all_seqs)
    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in all_seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
