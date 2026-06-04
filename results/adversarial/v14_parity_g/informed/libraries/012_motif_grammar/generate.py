#!/usr/bin/env python3
"""
012_motif_grammar — Synthetic sequences with controlled (motif_A, motif_B,
spacing, orientation) combinatorial coverage. Every ordered TF-motif pair
is represented ~255 times across diverse spacings and flanks.

Rationale: an MPRA model that learns TF cooperativity needs to see many
combinations of (which two motifs, how far apart, in what orientation).
Random co-occurrence undersamples rare pairs; explicit combinatorial
construction gives uniform coverage.
"""
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "sequences_0.txt"

L = 200
N_SEQ = 50_000
SEED = 0

MOTIFS = [
    "TGACTCA",      # AP-1
    "TGACGTCA",     # CREB/ATF
    "GGGACTTTCC",   # NF-kB
    "GGGGCGGGG",    # SP1
    "CACGTG",       # E-box
    "CCAATCA",      # NF-Y
    "GATAAG",       # GATA
    "AGGAAG",       # ETS
    "TTGCGCAAT",    # C/EBP
    "CAGCTG",       # E-box-like
    "GGGGAGGGG",    # KLF
    "AGGTCA",       # nuclear receptor
    "TGTGGTTT",     # HNF1
    "AACAAAG",      # HNF4
]
ACGT = "ACGT"
RC = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(RC)[::-1]


def bg(rng, n):
    return "".join(rng.choices(ACGT, k=n))


def main():
    rng = random.Random(SEED)
    n_motifs = len(MOTIFS)
    pairs = [(a, b) for a in range(n_motifs) for b in range(n_motifs)]
    reps_per_pair = N_SEQ // len(pairs)  # 50000 / 196 = 255
    leftover = N_SEQ - reps_per_pair * len(pairs)

    with open(OUT, "w") as f:
        out_count = 0
        for ai, bi in pairs:
            for _ in range(reps_per_pair):
                ma = MOTIFS[ai]
                mb = MOTIFS[bi]
                # random orientations
                if rng.random() < 0.5:
                    ma = rc(ma)
                if rng.random() < 0.5:
                    mb = rc(mb)
                # random spacing: 5-50bp between motifs
                gap = rng.randint(5, 50)
                inner = ma + bg(rng, gap) + mb
                flank_total = L - len(inner)
                if flank_total < 0:
                    # truncate gap if too long
                    gap = max(5, L - len(ma) - len(mb) - 10)
                    inner = ma + bg(rng, gap) + mb
                    flank_total = L - len(inner)
                left = rng.randint(0, flank_total)
                right = flank_total - left
                seq = bg(rng, left) + inner + bg(rng, right)
                assert len(seq) == L
                f.write(seq + "\n")
                out_count += 1
        # fill leftover with random pair
        for _ in range(leftover):
            ai, bi = rng.choice(pairs)
            ma = MOTIFS[ai]
            mb = MOTIFS[bi]
            gap = rng.randint(5, 50)
            inner = ma + bg(rng, gap) + mb
            flank_total = L - len(inner)
            left = rng.randint(0, flank_total)
            right = flank_total - left
            seq = bg(rng, left) + inner + bg(rng, right)
            f.write(seq + "\n")
            out_count += 1

    print(f"Wrote {out_count} motif-grammar sequences -> {OUT}")


if __name__ == "__main__":
    main()
