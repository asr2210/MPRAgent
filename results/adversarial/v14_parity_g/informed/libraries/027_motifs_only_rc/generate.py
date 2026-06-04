#!/usr/bin/env python3
"""
027_ccre_only_rc — Ablation: only cCRE (5 classes, 5k each) + RC.
25k unique cCRE + 25k RC = 50k. Pairs with 019 (synthetic-only+RC) to
isolate the contribution of real-genomic vs synthetic in flagship v3.
"""
import random
from collections import defaultdict
from pathlib import Path

from pyfaidx import Fasta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CCRE = ROOT / "data" / "GRCh38-cCREs.bed"
GENOME = ROOT / "data" / "hg38.fa"
OUT = HERE / "sequences_0.txt"

L = 200
PER_CLS = 5_000
SEED = 0
ACGT_SET = set("ACGT")
RC = str.maketrans("ACGT", "TGCA")

CLASSES = ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3")


def base_class(label):
    for p in label.split(","):
        if p in CLASSES:
            return p
    return None


def main():
    rng = random.Random(SEED)
    by_class = defaultdict(list)
    with open(CCRE) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            if "_" in chrom or chrom == "chrM":
                continue
            cls = base_class(parts[-1])
            if cls is None:
                continue
            try:
                start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue
            by_class[cls].append((chrom, (start + end) // 2))
    for cls in by_class:
        rng.shuffle(by_class[cls])

    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)
    parents = []
    for cls in CLASSES:
        taken = 0
        for chrom, mid in by_class[cls]:
            if taken >= PER_CLS:
                break
            s = mid - L // 2
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
            if len(seq) != L or not set(seq).issubset(ACGT_SET):
                continue
            parents.append(seq)
            taken += 1
        print(f"  {cls}: {taken}")

    out = []
    for s in parents:
        out.append(s)
        out.append(s.translate(RC)[::-1])
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} cCRE+RC -> {OUT}")


if __name__ == "__main__":
    main()
