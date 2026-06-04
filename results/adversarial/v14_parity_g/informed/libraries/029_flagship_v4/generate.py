#!/usr/bin/env python3
"""
029_flagship_v4 — Final flagship principled library.

Composition (50k total, every component RC-augmented):

  10k cCRE multi-class (PLS 3k, dELS 3k, pELS 1.5k, CTCF-only 1.25k, K4 1.25k)
  10k RC of above

   5k Sharpr top-activity (padded 145→200 with random flanks)
   5k RC of above
   5k Sharpr bottom-activity
   5k RC of above

   2.5k synthetic motif-packed (50-TF vocab, 4-8 motifs)
   2.5k RC of above

   2.5k motif-packed inside real cCRE backbones (mix of all 5 classes)
   2.5k RC of above

Design rationale (post-30-experiments findings):
  - v14 metric is uninformative; no library scored outside ±0.005 (noise).
  - Same recipe gives different scores across seeds: variance ~0.001 between
    019 (+0.0025) and 022 (+0.0013). Score is pure sample noise.
  - Therefore optimize for biological coverage that *would* train a good
    MPRA model in any genuine setting:
      * regulatory element breadth (5 cCRE classes covers promoters,
        proximal+distal enhancers, insulators, broad-H3K4me3 marks)
      * measured activity anchors (Sharpr poles)
      * combinatorial TF grammar (50-TF synthetic library)
      * authentic local context for motif inserts (motif-in-cCRE hybrid)
      * strand augmentation (RC) for invariance
  - No replication, no extreme distributional bias.
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

# Per-class cCRE budgets (parents; doubled by RC)
CCRE_BUDGET = {"PLS": 3_000, "dELS": 3_000, "pELS": 1_500,
               "CTCF-only": 1_250, "DNase-H3K4me3": 1_250}
N_SHARPR_TOP = 5_000
N_SHARPR_BOT = 5_000
N_SYNTH_MOTIF = 2_500
N_MOTIF_IN_CCRE = 2_500

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
CLASSES = ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3")


def resolve_iupac(s, rng):
    return "".join(rng.choice(IUPAC[c]) if c in IUPAC else c for c in s)


def base_class(label):
    for p in label.split(","):
        if p in CLASSES:
            return p
    return None


def bg(rng, n):
    return "".join(rng.choices(ACGT, k=n))


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


def synth_motif(rng, motifs_clean):
    return insert_motifs(bg(rng, L), rng, motifs_clean)


def main():
    rng = random.Random(SEED)
    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)
    motifs_clean = [resolve_iupac(m, rng) for m in MOTIFS_RAW]

    # cCRE pool by class
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

    def extract(coords, n):
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
        return out

    ccre_parents = []
    for cls, n in CCRE_BUDGET.items():
        got = extract(by_class[cls], n)
        ccre_parents += got
        print(f"  cCRE-{cls}: {len(got)}")

    # Backbones for motif-in-cCRE (extra pool, balanced across classes)
    per_cls_bb = N_MOTIF_IN_CCRE // len(CLASSES) + 50
    bb_pool = []
    for cls in CLASSES:
        # skip parents already taken
        skip = CCRE_BUDGET[cls]
        bb_pool += extract(by_class[cls][skip:], per_cls_bb)
    rng.shuffle(bb_pool)
    bb_pool = bb_pool[:N_MOTIF_IN_CCRE]
    print(f"  motif-in-cCRE backbones: {len(bb_pool)}")

    # Sharpr top/bot
    with h5py.File(SHARPR, "r") as f:
        Y = f["Y/output"][:]
        mean_act = Y.mean(axis=1)
        order = np.argsort(mean_act)
        sel = np.concatenate([order[:N_SHARPR_BOT], order[-N_SHARPR_TOP:]])
        sel_sorted = np.sort(sel)
        pos = {orig: i for i, orig in enumerate(sel_sorted)}
        X_sel = f["X/sequence"][sel_sorted]
    top_positions = [pos[i] for i in order[-N_SHARPR_TOP:]]
    bot_positions = [pos[i] for i in order[:N_SHARPR_BOT]]
    idx = X_sel.argmax(axis=-1)
    has_base = X_sel.sum(axis=-1) > 0

    def decode(positions):
        seqs = []
        for p in positions:
            chars = [ACGT[idx[p, j]] if has_base[p, j] else rng.choice(ACGT)
                     for j in range(SHARPR_L)]
            left = "".join(rng.choices(ACGT, k=PAD_LEFT))
            right = "".join(rng.choices(ACGT, k=PAD_RIGHT))
            seqs.append(left + "".join(chars) + right)
        return seqs

    sharpr_top = decode(top_positions)
    sharpr_bot = decode(bot_positions)
    print(f"  Sharpr top/bot: {len(sharpr_top)}/{len(sharpr_bot)}")

    # Synth motif library
    synth_seqs = [synth_motif(rng, motifs_clean) for _ in range(N_SYNTH_MOTIF)]
    motif_in_bb = [insert_motifs(s, rng, motifs_clean) for s in bb_pool]

    def with_rc(seqs):
        return seqs + [s.translate(RC)[::-1] for s in seqs]

    out = []
    out += with_rc(ccre_parents)
    out += with_rc(sharpr_top)
    out += with_rc(sharpr_bot)
    out += with_rc(synth_seqs)
    out += with_rc(motif_in_bb)

    print(f"  total: {len(out)}")
    assert len(out) == 50_000, f"got {len(out)}"
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} flagship-v4 -> {OUT}")


if __name__ == "__main__":
    main()
