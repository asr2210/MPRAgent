#!/usr/bin/env python3
"""
023_flagship_v1 — Principled flagship library v1. Built from defensible
ML/biology priors, regardless of (uninformative) v14 metric.

Composition (50k total):
  - 10k cCRE PLS (promoters)
  - 10k cCRE dELS (distal enhancers)
  - 5k cCRE pELS (proximal enhancers)
  - 2.5k cCRE CTCF-only (insulators)
  - 2.5k cCRE DNase-H3K4me3 (broad chromatin)
  - 5k Sharpr-MPRA top by activity (real strong enhancers)
  - 5k Sharpr-MPRA bottom by activity (real silencers / weak elements)
  - 5k synthetic motif-packed (50-TF vocab, 4-8 motifs)
  - 5k reverse-complements of above synthetic for strand invariance

Rationale:
  * Major real-element class coverage (cCRE breadth) — generalization
    across regulatory contexts.
  * Real measured activity poles (Sharpr) — anchors for activity prediction.
  * Synthetic motif-packed sequences — controlled exposure to TF grammar
    above the natural occurrence rate.
  * Strand augmentation — strand-invariant motif representations.
  * No replication, no extreme distributional bias.
"""
import gzip
import random
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
from pyfaidx import Fasta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "data"
CCRE = DATA / "GRCh38-cCREs.bed"
SHARPR = DATA / "sharpr_train.hdf5"
GENOME = DATA / "hg38.fa"
OUT = HERE / "sequences_0.txt"

L = 200
SEED = 0
ACGT_SET = set("ACGT")
ACGT = "ACGT"
SHARPR_L = 145
PAD_LEFT = 27
PAD_RIGHT = 28

# Composition
N_PLS = 10_000
N_DELS = 10_000
N_PELS = 5_000
N_CTCF = 2_500
N_K4 = 2_500
N_SHARPR_TOP = 5_000
N_SHARPR_BOT = 5_000
N_MOTIF = 5_000
N_MOTIF_RC = 5_000

MOTIFS = [
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


def resolve_iupac(s, rng):
    return "".join(rng.choice(IUPAC[c]) if c in IUPAC else c for c in s)


def bg(rng, n):
    return "".join(rng.choices(ACGT, k=n))


def base_class(label):
    for p in label.split(","):
        if p in ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"):
            return p
    return None


def load_ccre_pool(rng):
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
    return by_class


def extract(coords, fa, n, label):
    out = []
    for chrom, mid in coords:
        if len(out) >= n:
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
        out.append(seq)
    print(f"  {label}: {len(out)}/{n}")
    return out


def load_sharpr_pairs(rng, n_top, n_bot):
    with h5py.File(SHARPR, "r") as f:
        Y = f["Y/output"][:]
        mean_act = Y.mean(axis=1)
        order = np.argsort(mean_act)
        bot_idx = order[:n_bot]
        top_idx = order[-n_top:]
        all_idx = np.concatenate([bot_idx, top_idx])
        sorted_idx = np.sort(all_idx)
        pos = {orig: i for i, orig in enumerate(sorted_idx)}
        X_sel = f["X/sequence"][sorted_idx]
    top_in_sel = [pos[i] for i in top_idx]
    bot_in_sel = [pos[i] for i in bot_idx]
    idx = X_sel.argmax(axis=-1)
    has_base = X_sel.sum(axis=-1) > 0

    def decode(positions):
        out = []
        for p in positions:
            chars = [ACGT[idx[p, j]] if has_base[p, j] else rng.choice(ACGT)
                     for j in range(SHARPR_L)]
            left = "".join(rng.choices(ACGT, k=PAD_LEFT))
            right = "".join(rng.choices(ACGT, k=PAD_RIGHT))
            out.append(left + "".join(chars) + right)
        return out

    return decode(top_in_sel), decode(bot_in_sel)


def motif_packed(rng, motifs_clean):
    n_motifs = rng.randint(4, 8)
    seq = list(bg(rng, L))
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
    motifs_clean = [resolve_iupac(m, rng) for m in MOTIFS]

    print("Loading cCRE classes...")
    by_class = load_ccre_pool(rng)
    out = []
    out += extract(by_class["PLS"], fa, N_PLS, "cCRE-PLS")
    out += extract(by_class["dELS"], fa, N_DELS, "cCRE-dELS")
    out += extract(by_class["pELS"], fa, N_PELS, "cCRE-pELS")
    out += extract(by_class["CTCF-only"], fa, N_CTCF, "cCRE-CTCF-only")
    out += extract(by_class["DNase-H3K4me3"], fa, N_K4, "cCRE-DNase-H3K4me3")

    print("Loading Sharpr poles...")
    s_top, s_bot = load_sharpr_pairs(rng, N_SHARPR_TOP, N_SHARPR_BOT)
    out += s_top + s_bot
    print(f"  Sharpr top: {len(s_top)}, bot: {len(s_bot)}")

    print("Generating synthetic motif-packed...")
    motifs_seqs = [motif_packed(rng, motifs_clean) for _ in range(N_MOTIF)]
    out += motifs_seqs
    out += [s.translate(RC)[::-1] for s in motifs_seqs[:N_MOTIF_RC]]
    print(f"  motif+RC: {N_MOTIF + N_MOTIF_RC}")

    assert len(out) == 50_000, f"got {len(out)}"
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} flagship-v1 sequences -> {OUT}")


if __name__ == "__main__":
    main()
