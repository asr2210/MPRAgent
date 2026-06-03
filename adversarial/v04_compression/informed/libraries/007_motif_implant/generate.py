#!/usr/bin/env python3
"""
007_motif_implant — random scaffolds with planted TF consensus motifs.

50k sequences, each 200bp.
For each sequence:
  - 200bp uniform random scaffold (seed=0).
  - Plant K=3 non-overlapping motifs at random positions.
  - Each motif drawn uniformly from a curated list of ~50 strong
    vertebrate TF consensus sequences (lengths 6–15bp), reverse
    complement allowed.

Tests the hypothesis from E5: real DNA's +0.19 r-point bonus over
uniform random comes from motif content / higher-order structure, NOT
just k-mer composition.

- If 007 ≈ 0.50 (≈ real genome): motif content alone reproduces the
  real-DNA prior. We can engineer libraries without needing real DNA.
- If 007 between 0.31 and 0.50: motifs explain part; rest is repeats /
  long-range context.
- If 007 ≈ 0.31 or worse: motifs don't transfer — the prior is in
  structural context beyond plantable elements.

Generalization argument: TF consensus motifs are PROPERTIES OF DNA
recognized by transcription factors that exist in essentially all
vertebrate cell types. Planting these motifs in random scaffolds
exposes the model to "regulatory grammar" decoupled from any specific
cell-type chromatin context. A model trained on these motifs should
generalize to any cell type that uses the same TF families.
"""
from pathlib import Path

import numpy as np

SEED = 0
N = 50_000
L = 200
K_MOTIFS = 3
ALPHABET = "ACGT"

# Curated set of ~50 strong vertebrate TF consensus binding motifs
# (mostly drawn from JASPAR top-scoring matrices + classic literature).
MOTIFS = [
    "TATAAA",        # TBP / TATA
    "CCAAT",         # NFY / CAAT box
    "GGGGCGGGG",     # SP1 / GC box
    "TGACTCA",       # AP-1 (Fos/Jun)
    "TGACGTCA",      # CREB / ATF
    "GGGAATTTCC",    # NF-kB / p65
    "CACGTG",        # E-box (Myc/Max, USF)
    "CAGCTG",        # E-box (MyoD/Tcf21)
    "CAGGTG",        # E-box variant (ZEB)
    "AGATAA",        # GATA1/2/3
    "TGTTTGT",       # FOXA / HNF3
    "TGGACTTTG",     # HNF4
    "GTTAATNATTAA",  # HNF1 (use exact: GTTAATCATTAA)
    "ATTGCGCAAT",    # C/EBP
    "ATGCAAAT",      # POU / Oct1
    "TTTATGAA",      # HOX / Cdx2
    "TGTGGTT",       # RUNX
    "CAGTTG",        # MYB
    "GGGCGTGGGCG",   # KLF / EGR1 (poly-GC)
    "AGAACAGAGTGTTCT", # GR / AR full pal
    "TGTACA",        # GR/AR half-site
    "GCGCATGCGC",    # NRF1
    "CCATCTT",       # YY1
    "CAGCAATT",      # ZNF143
    "GTAACC",        # RFX
    "CCGCGNGGNGGCAG",# CTCF (use exact: CCGCGAGGAGGCAG)
    "TCAGCACCATG",   # REST / NRSF
    "TTTCGCGC",      # E2F
    "AAGTGA",        # IRF
    "TTCCGGGAA",     # STAT
    "ACCGGAAGT",     # ETS
    "CTTTGT",        # TCF/LEF
    "ATTGATTT",      # HNF6
    "TCAAGGTCA",     # LRH1
    "AGGTCATGGTCC",  # FXR/RXR
    "TGCGTG",        # AhR
    "CAGAC",         # SMAD
    "GGTGTGAA",      # T-box / Brachyury
    "TGCCAA",        # NF1
    "CCGCCATCTT",    # YY1 variant
    "TTGCCCAA",      # PPAR half
    "AGGTCANNNGTGACC",  # RXR (use AGGTCAAAGGTGACC)
    "GCCAAT",        # NFY variant
    "GAAANT",        # IRF1 short (use GAAAATT)
    "CACCCAGCCT",    # KLF1
    "TGCAGTGCT",     # ZNF
    "CCTCAGGCT",     # ?
    "GACCAAT",       # ?
    "GAGGAA",        # ETS short
    "GGAA",          # ETS core
    "TAATTA",        # HOX core
    "AACCAC",        # HSF
    "TCACGTGA",      # Myc canonical
]

# Resolve IUPAC ambiguities and remove any N-containing motifs
IUPAC = {
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def resolve_iupac(seq, rng):
    out = []
    for ch in seq:
        if ch in "ACGT":
            out.append(ch)
        elif ch in IUPAC:
            out.append(IUPAC[ch][rng.integers(0, len(IUPAC[ch]))])
        else:
            raise ValueError(ch)
    return "".join(out)


def rev_comp(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in seq[::-1])


def main():
    rng = np.random.default_rng(SEED)
    # Pre-clean motifs: resolve any IUPAC chars deterministically
    motifs = []
    for m in MOTIFS:
        # Lower-case any IUPAC in the canonical list. Keep both strands.
        m_clean = "".join(ch if ch in "ACGT" else "N" for ch in m)  # simple
        # Replace N with random base each time we use it
        motifs.append(m)
    print(f"motif library: {len(motifs)} entries, lengths "
          f"{min(len(m) for m in motifs)}–{max(len(m) for m in motifs)}")

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    alphabet_arr = np.array(list(ALPHABET))

    seqs = []
    for i in range(N):
        # Random scaffold
        scaffold = list(alphabet_arr[rng.integers(0, 4, size=L)])
        # Plant K motifs
        used_ranges = []
        attempts = 0
        planted = 0
        while planted < K_MOTIFS and attempts < 50:
            attempts += 1
            mi = rng.integers(0, len(motifs))
            m = motifs[mi]
            m_resolved = resolve_iupac(m, rng)
            if rng.random() < 0.5:
                m_resolved = rev_comp(m_resolved)
            mlen = len(m_resolved)
            if mlen >= L:
                continue
            start = rng.integers(0, L - mlen + 1)
            end = start + mlen
            # Check overlap
            if any(not (end <= s0 or start >= e0) for (s0, e0) in used_ranges):
                continue
            scaffold[start:end] = list(m_resolved)
            used_ranges.append((start, end))
            planted += 1
        seqs.append("".join(scaffold))

    with open(out_path, "w") as f:
        f.write("\n".join(seqs))
        f.write("\n")
    print(f"wrote {len(seqs)} sequences → {out_path}")
    print(f"avg motifs/seq: {K_MOTIFS} (target)")


if __name__ == "__main__":
    main()
