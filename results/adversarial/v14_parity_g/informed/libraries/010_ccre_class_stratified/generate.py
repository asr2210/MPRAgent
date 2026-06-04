#!/usr/bin/env python3
"""
010_ccre_class_stratified — Equal-proportion sample across the major
ENCODE cCRE V3 functional classes. Balances regulatory element types so
the model sees comparable counts of each.

Classes (canonical, ignoring CTCF-bound suffix):
  PLS, pELS, dELS, CTCF-only, DNase-H3K4me3

5 classes × 10,000 sequences = 50,000.

Rationale (principled-design phase): a real-world MPRA model that generalizes
across regulatory contexts needs balanced exposure to promoters, proximal
enhancers, distal enhancers, CTCF-bound insulators, and broad-H3K4me3
elements. Natural cCRE frequencies are heavily skewed toward dELS — stratifying
prevents that bias from dominating training.
"""
import random
from collections import defaultdict
from pathlib import Path

from pyfaidx import Fasta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "data"
CCRE = DATA / "GRCh38-cCREs.bed"
GENOME = DATA / "hg38.fa"
OUT = HERE / "sequences_0.txt"

L = 200
PER_CLASS = 10_000
SEED = 0
ACGT = set("ACGT")


def base_class(label):
    parts = label.split(",")
    for p in parts:
        if p in ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"):
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
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            cls = base_class(parts[-1])
            if cls is None:
                continue
            mid = (start + end) // 2
            by_class[cls].append((chrom, mid))

    for cls in by_class:
        rng.shuffle(by_class[cls])
        print(f"  {cls}: {len(by_class[cls])} candidates")

    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)
    out = []
    for cls in ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"):
        taken = 0
        for chrom, mid in by_class[cls]:
            if taken >= PER_CLASS:
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
            if len(seq) != L or not set(seq).issubset(ACGT):
                continue
            out.append(seq)
            taken += 1
        print(f"  {cls}: took {taken}")

    assert len(out) == PER_CLASS * 5, f"got {len(out)}"
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} -> {OUT}")


if __name__ == "__main__":
    main()
