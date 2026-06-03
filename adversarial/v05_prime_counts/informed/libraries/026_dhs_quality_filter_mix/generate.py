"""
026_dhs_quality_filter_mix
==========================
33/67 mix where bio half is dhs_specific FILTERED to numsamples<=10 AND
mean_signal in top 50%. Combines specificity AND signal-strength criteria.

Tests if "best regulatory elements" (specific AND strong) beat plain
dhs_specific in the bio anchor.
"""
import gzip
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
N_BIO = 16_667
N_UNIFORM = 50_000 - N_BIO
N_SEQS = 50_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
DHS_PATH = DATA_DIR / "dhs_index.txt.gz"
HG38_PATH = DATA_DIR / "hg38.fa"
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def load_dhs():
    chroms, summits, signals, nsamples = [], [], [], []
    with gzip.open(DHS_PATH, "rt") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] not in CANONICAL:
                continue
            chroms.append(parts[0])
            summits.append(int(parts[6]))
            signals.append(float(parts[4]))
            nsamples.append(int(parts[5]))
    return (chroms, np.array(summits), np.array(signals),
            np.array(nsamples, dtype=np.int32))


def main():
    rng = np.random.default_rng(SEED)
    chroms, summits, signals, nsamples = load_dhs()
    sig_med = np.median(signals)
    keep = (nsamples <= 10) & (signals >= sig_med)
    print(f"Filter: nsamp<=10 & sig>={sig_med:.2f} -> {keep.sum()} sites",
          file=sys.stderr)
    idx = np.where(keep)[0]
    sub_chroms = [chroms[i] for i in idx]
    sub_summits = summits[idx]
    sub_n = nsamples[idx]
    w = 1.0 / np.sqrt(sub_n.astype(np.float64))
    w /= w.sum()
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    bio_seqs = []
    pool = rng.choice(len(idx), size=int(N_BIO * 1.4), replace=False, p=w)
    for i in pool:
        if len(bio_seqs) >= N_BIO:
            break
        chrom = sub_chroms[i]
        summit = int(sub_summits[i])
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
        bio_seqs.append(s)
    assert len(bio_seqs) == N_BIO, len(bio_seqs)

    bases = np.array(list("ACGT"))
    idx_u = rng.integers(0, 4, size=(N_UNIFORM, SEQ_LEN))
    uniform_seqs = ["".join(row) for row in bases[idx_u]]

    all_seqs = bio_seqs + uniform_seqs
    rng.shuffle(all_seqs)
    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in all_seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
