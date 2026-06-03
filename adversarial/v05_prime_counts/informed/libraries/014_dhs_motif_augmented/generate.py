"""
014_dhs_motif_augmented
=======================
Take dhs_specific sequences (best library so far, 0.049) and PLANT
strong JASPAR motifs at random positions in each. 4 motifs per sequence.

Hypothesis: motif content explained ~0.012 of dhs_specific score
(Exp 013). Adding MORE motifs on top of an already-motif-rich
DHS background should push eval_01 above 0.05 for the first time.
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
N_MOTIFS_PER_SEQ = 4
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
DHS_PATH = DATA_DIR / "dhs_index.txt.gz"
HG38_PATH = DATA_DIR / "hg38.fa"
JASPAR_PATH = DATA_DIR / "jaspar2024_core.txt"
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def parse_jaspar(path):
    pfms = []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith(">"):
            rows = []
            for j in range(1, 5):
                row = lines[i + j].strip()
                start = row.index("[") + 1
                end = row.index("]")
                rows.append([float(x) for x in row[start:end].split()])
            mat = np.array(rows, dtype=np.float64)
            col_sums = mat.sum(axis=0, keepdims=True)
            col_sums[col_sums == 0] = 1.0
            mat = (mat + 0.01) / (col_sums + 0.04)
            if 6 <= mat.shape[1] <= 16:
                pfms.append(mat)
            i += 5
        else:
            i += 1
    return pfms


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
    bases = np.array(["A", "C", "G", "T"])
    pfms = parse_jaspar(JASPAR_PATH)
    print(f"  {len(pfms)} PFMs", file=sys.stderr)
    chroms, summits, nsamples = load_dhs()
    w = 1.0 / np.sqrt(nsamples.astype(np.float64))
    w /= w.sum()
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    pool = rng.choice(len(chroms), size=int(N_SEQS * 1.3), replace=False, p=w)
    seqs = []
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
        # plant N_MOTIFS_PER_SEQ motifs
        seq_arr = list(s)
        occ = []
        chosen = rng.choice(len(pfms), size=N_MOTIFS_PER_SEQ, replace=True)
        for k in chosen:
            pfm = pfms[k]
            ml = pfm.shape[1]
            inst = "".join(rng.choice(bases, p=pfm[:, c]) for c in range(ml))
            for _try in range(40):
                st = int(rng.integers(0, SEQ_LEN - ml + 1))
                en = st + ml
                if all(en <= a or st >= b for a, b in occ):
                    break
            else:
                st = int(rng.integers(0, SEQ_LEN - ml + 1))
                en = st + ml
            occ.append((st, en))
            for off, ch in enumerate(inst):
                seq_arr[st + off] = ch
        seqs.append("".join(seq_arr))
        if len(seqs) % 10_000 == 0:
            print(f"  {len(seqs)}/{N_SEQS}", file=sys.stderr)
    assert len(seqs) == N_SEQS
    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
