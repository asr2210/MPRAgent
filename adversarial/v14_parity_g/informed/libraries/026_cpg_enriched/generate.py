#!/usr/bin/env python3
"""
026_cpg_enriched — cCRE-PLS (promoter-like) sequences that fall inside or
adjacent to CpG islands, identified by elevated observed/expected CpG ratio
in the 200bp window. 50k high-CpG cCRE sequences.

Rationale: CpG islands mark active promoters; their characteristic sequence
composition (high GC, high CpG dinucleotide frequency) is a strong feature
that no other library tested has focused on.
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
N_SEQ = 50_000
SEED = 0
ACGT_SET = set("ACGT")


def cpg_oe(seq):
    # Observed/expected CpG: (CpG count * L) / (C count * G count)
    c = seq.count("C")
    g = seq.count("G")
    if c == 0 or g == 0:
        return 0.0
    cpg = sum(1 for i in range(len(seq) - 1) if seq[i] == "C" and seq[i+1] == "G")
    return cpg * len(seq) / (c * g)


def gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)


def base_class(label):
    for p in label.split(","):
        if p in ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"):
            return p
    return None


def main():
    rng = random.Random(SEED)
    # Take ALL cCREs (any class) — CpG-island filter will pick promoter-like
    coords = []
    with open(CCRE) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            if "_" in chrom or chrom == "chrM":
                continue
            try:
                start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue
            coords.append((chrom, (start + end) // 2))
    rng.shuffle(coords)
    print(f"  candidate cCRE pool: {len(coords)}")

    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)
    out = []
    for chrom, mid in coords:
        if len(out) >= N_SEQ:
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
        # CpG-island-like criterion: GC > 0.55 and CpG O/E > 0.65
        if gc(seq) > 0.55 and cpg_oe(seq) > 0.65:
            out.append(seq)
    print(f"  passed CpG-island filter: {len(out)}")
    if len(out) < N_SEQ:
        # If we ran out of candidates, fill with relaxed filter
        for chrom, mid in coords:
            if len(out) >= N_SEQ:
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
            if seq in out:
                continue
            out.append(seq)

    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out[:N_SEQ]:
            f.write(s + "\n")
    print(f"Wrote {N_SEQ} CpG-enriched cCRE sequences -> {OUT}")


if __name__ == "__main__":
    main()
