#!/usr/bin/env python3
"""
025_flagship_v3 — Motifs embedded INSIDE real cCRE backbones, plus
Sharpr poles, plus RC. Hybridizes 020 (motif-in-cCRE) + 019 (RC) +
011 (Sharpr inclusion).

Composition (50k total):
  15k cCRE-backbone + 4-8 motif inserts (50-vocab)
  15k RC of those motif-in-cCRE
   5k Sharpr top
   5k Sharpr bottom
   5k cCRE unmodified (PLS+dELS mix)
   5k RC of cCRE unmodified
"""
import random
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
from pyfaidx import Fasta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CCRE = ROOT / "data" / "GRCh38-cCREs.bed"
SHARPR = ROOT / "data" / "sharpr_train.hdf5"
GENOME = ROOT / "data" / "hg38.fa"
OUT = HERE / "sequences_0.txt"

L = 200
SEED = 0
ACGT_SET = set("ACGT")
ACGT = "ACGT"
SHARPR_L = 145
PAD_LEFT = 27
PAD_RIGHT = 28

N_CCRE_BB = 15_000
N_CCRE_UNMOD = 5_000
N_SHARPR_TOP = 5_000
N_SHARPR_BOT = 5_000

MOTIFS_RAW = [
    "TGACTCA", "TGAGTCA", "TGACGTCA", "GGGACTTTCC", "GGGAATTCCC",
    "GGGGCGGGG", "GGGGGCGGGG", "CACGTG", "CAGGTG", "CAGCTG",
    "CACGCG", "CCAATCA", "ATTGGC", "GATAAG", "AGATAA",
    "GGAAGT", "AGGAAG", "CGGAAG", "TTGCGCAAT", "ATTGCGCAA",
    "GTTAATNATTAAC", "TTAATGA", "TGTGGTTT", "AACAAAG", "AGGTCA",
    "TGACCT", "AGAACA", "TGTTCT", "GTAAACA", "TGTTTAC",
    "AACAATG", "CATTGTT", "TAATTA", "TAATGG", "CCATTA",
    "GGGGAGGGG", "CCCCGCCC", "ATTTGCAT", "ATGCAAAT", "AANAGTGT",
    "ACACTTNNT", "TGASTCAGCA", "TTCNNNGAA", "TTCCGGGAA", "CAGGAAG",
    "ACCGGAAG", "TGTGGAAA", "TTTCCACA", "CATATG", "GCCNNNGGC",
]
IUPAC = {"N": "ACGT", "S": "GC", "R": "AG", "Y": "CT", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}
RC = str.maketrans("ACGT", "TGCA")

CCRE_CLASSES = ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3")


def resolve_iupac(s, rng):
    return "".join(rng.choice(IUPAC[c]) if c in IUPAC else c for c in s)


def base_class(label):
    for p in label.split(","):
        if p in CCRE_CLASSES:
            return p
    return None


def insert_motifs(seq_str, rng, motifs_clean):
    seq = list(seq_str)
    n_motifs = rng.randint(4, 8)
    occupied = []
    for _ in range(n_motifs):
        m = rng.choice(motifs_clean)
        for _ in range(20):
            start = rng.randint(0, L - len(m))
            end = start + len(m)
            if all(end <= s or start >= e for s, e in occupied):
                for i, ch in enumerate(m):
                    seq[start + i] = ch
                occupied.append((start, end))
                break
    return "".join(seq)


def main():
    rng = random.Random(SEED)
    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)
    motifs_clean = [resolve_iupac(m, rng) for m in MOTIFS_RAW]

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

    # Build a pool of cCRE genomic sequences (mix of all classes, balanced)
    needed = N_CCRE_BB + N_CCRE_UNMOD
    per_cls = needed // len(CCRE_CLASSES) + 100  # slack
    pool = []
    for cls in CCRE_CLASSES:
        taken = 0
        for chrom, mid in by_class[cls]:
            if taken >= per_cls:
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
            pool.append(seq)
            taken += 1
    rng.shuffle(pool)
    print(f"  cCRE pool: {len(pool)}")

    out = []
    # 15k cCRE-backbone + motifs
    bb = pool[:N_CCRE_BB]
    bb_motif = [insert_motifs(s, rng, motifs_clean) for s in bb]
    out += bb_motif
    out += [s.translate(RC)[::-1] for s in bb_motif]
    print(f"  motif-in-cCRE: {len(bb_motif)} + RC")

    # 5k cCRE unmodified
    unmod = pool[N_CCRE_BB:N_CCRE_BB + N_CCRE_UNMOD]
    out += unmod
    out += [s.translate(RC)[::-1] for s in unmod]
    print(f"  cCRE unmod: {len(unmod)} + RC")

    # 5k+5k Sharpr poles
    with h5py.File(SHARPR, "r") as f:
        Y = f["Y/output"][:]
        mean_act = Y.mean(axis=1)
        order = np.argsort(mean_act)
        sel = np.concatenate([order[:N_SHARPR_BOT], order[-N_SHARPR_TOP:]])
        sel_sorted = np.sort(sel)
        X_sel = f["X/sequence"][sel_sorted]
    idx = X_sel.argmax(axis=-1)
    has_base = X_sel.sum(axis=-1) > 0
    sharpr_seqs = []
    for i in range(idx.shape[0]):
        chars = [ACGT[idx[i, j]] if has_base[i, j] else rng.choice(ACGT)
                 for j in range(SHARPR_L)]
        left = "".join(rng.choices(ACGT, k=PAD_LEFT))
        right = "".join(rng.choices(ACGT, k=PAD_RIGHT))
        sharpr_seqs.append(left + "".join(chars) + right)
    out += sharpr_seqs
    print(f"  Sharpr poles: {len(sharpr_seqs)}")

    assert len(out) == 50_000, f"got {len(out)}"
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} flagship-v3 -> {OUT}")


if __name__ == "__main__":
    main()
