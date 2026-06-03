"""
027_dhs_uniform_gc_matched
==========================
Same as Exp 019 (33dhs/67uniform) but uniform half uses GC content
matched to human regulatory DNA (~46% GC, vs 50% uniform).

Tests: does eval_08 reward "random sequence" specifically (GC-neutral)
or is it just "non-genomic"? If uniform_GC45 still scores high eval_08,
composition is fine. If it drops, eval_08 wants 50/50 base distribution.
"""
import gzip
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
N_DHS = 16_667
N_UNIFORM = 50_000 - N_DHS
N_SEQS = 50_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
GC_FRAC = 0.46
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

    dhs_seqs = []
    pool = rng.choice(len(chroms), size=int(N_DHS * 1.5), replace=False, p=w)
    for i in pool:
        if len(dhs_seqs) >= N_DHS:
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
    assert len(dhs_seqs) == N_DHS

    # GC-matched uniform: p(G)=p(C)=GC/2, p(A)=p(T)=(1-GC)/2
    p = np.array([(1 - GC_FRAC) / 2, GC_FRAC / 2, GC_FRAC / 2, (1 - GC_FRAC) / 2])
    bases = np.array(list("ACGT"))
    idx = rng.choice(4, size=(N_UNIFORM, SEQ_LEN), p=p)
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
