#!/usr/bin/env python3
"""
007_dhs_bimodal_signal — 25k DHS with highest mean_signal + 25k with lowest.
Maximum biological-activity contrast in one library. Tests whether v14's
model responds to activity-level differences at all.
"""
import gzip
from pathlib import Path
import random

from pyfaidx import Fasta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "data"
DHS = DATA / "DHS_Index_hg38.txt.gz"
GENOME = DATA / "hg38.fa"
OUT = HERE / "sequences_0.txt"

L = 200
N_HALF = 25_000
SEED = 0
ACGT = set("ACGT")


def main():
    rows = []
    with gzip.open(DHS, "rt") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            if "_" in chrom or chrom == "chrM":
                continue
            try:
                summit = int(parts[6])
                signal = float(parts[4])
            except ValueError:
                continue
            rows.append((signal, chrom, summit))
    print(f"Loaded {len(rows)} DHS elements")

    # Sort by signal
    rows.sort(key=lambda r: r[0])
    # Bottom and top pools — over-grab in case of N rejections
    bot_pool = rows[: N_HALF * 2]
    top_pool = rows[-N_HALF * 2:]

    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)

    def extract(pool, n_target, label):
        out = []
        for sig, chrom, summit in pool:
            if len(out) >= n_target:
                break
            s = summit - L // 2
            e = s + L
            if s < 0:
                continue
            try:
                clen = len(fa[chrom])
            except KeyError:
                continue
            if e > clen:
                continue
            seq = str(fa[chrom][s:e]).upper()
            if len(seq) != L or not set(seq).issubset(ACGT):
                continue
            out.append(seq)
        print(f"  {label}: {len(out)}/{n_target}")
        return out

    high = extract(top_pool, N_HALF, "high")
    low = extract(bot_pool, N_HALF, "low")
    assert len(high) == N_HALF and len(low) == N_HALF
    # Interleave
    out = []
    for h, l in zip(high, low):
        out.append(h)
        out.append(l)
    rng = random.Random(SEED)
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} sequences -> {OUT}")


if __name__ == "__main__":
    main()
