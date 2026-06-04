#!/usr/bin/env python3
"""
003_motif_packed — 50,000 synthetic 200bp sequences each containing several
strong TF binding sites embedded in random GC-50 background.

Test: does packing each sequence with strong known motifs let the model learn
a motif→activity mapping under v14's apparently brief training? If MPRA
activity is driven by motif content (the standard hypothesis), this library
should give the model very clean signal — each motif's coefficient can be
learned from many copies across the 50k sequences.

We use a curated set of strong consensus motifs for TFs known to be active in
K562, HepG2, and SK-N-SH plus tissue-agnostic factors.
"""
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "sequences_0.txt"

L = 200
N_SEQ = 50_000
SEED = 0

# Strong consensus motifs (single best instance). Mix of ubiquitous and
# cell-type-enriched factors so model gets diverse signal.
MOTIFS = [
    "TGACTCA",      # AP-1 (FOS/JUN) — ubiquitous
    "TGACGTCA",     # CREB/ATF
    "GGGACTTTCC",   # NF-kB
    "GGGGCGGGG",    # SP1
    "TATAAA",       # TATA box
    "CACGTG",       # E-box (MYC/MAX)
    "CCAATCA",      # NF-Y / CCAAT
    "GATAAG",       # GATA1 (K562 erythroid)
    "AGGAAG",       # ETS family
    "TTGCGCAAT",    # C/EBP (HepG2 liver)
    "TGAACTTTGAACT",# HNF4 (HepG2)
    "GTTAATCATTAAC",# HNF1 (HepG2)
    "CAGCTG",       # E-box variant / NEUROD1 (SK-N-SH)
    "GGGGAGGGG",    # KLF family
    "AACAAAG",      # OCT/POU
    "AGGTCA",       # nuclear receptor half-site
    "TGTGGTTT",     # FOXA
]

ACGT = "ACGT"


def random_bg(rng, n):
    return "".join(rng.choices(ACGT, k=n))


def make_sequence(rng):
    # Choose 4-8 motif instances, place them at non-overlapping random
    # positions in a random background.
    n_motifs = rng.randint(4, 8)
    chosen = [rng.choice(MOTIFS) for _ in range(n_motifs)]
    # Sort by length to pack greedily
    rng.shuffle(chosen)
    # Find placements
    placements = []
    # Try many times to fit
    occupied = []
    for m in chosen:
        for _ in range(20):
            start = rng.randint(0, L - len(m))
            end = start + len(m)
            ok = all(end <= s or start >= e for s, e in occupied)
            if ok:
                placements.append((start, m))
                occupied.append((start, end))
                break
    # Fill background
    seq = list(random_bg(rng, L))
    for start, m in placements:
        for i, ch in enumerate(m):
            seq[start + i] = ch
    return "".join(seq)


def main():
    rng = random.Random(SEED)
    with open(OUT, "w") as f:
        for _ in range(N_SEQ):
            f.write(make_sequence(rng) + "\n")
    print(f"Wrote {N_SEQ} motif-packed sequences -> {OUT}")


if __name__ == "__main__":
    main()
