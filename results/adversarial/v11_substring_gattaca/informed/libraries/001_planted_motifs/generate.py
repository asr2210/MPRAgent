"""
E1: planted_motifs

Random 50% GC backbone with planted canonical TF motifs from diverse TF
families. Tests whether explicit motif content improves on random_uniform.

Library design:
- 50,000 sequences, 200 bp
- 80% sequences receive 1-4 planted motifs at random positions / orientations
- Motif set covers diverse TF families (chromatin, lineage, basal, immediate-early, etc.)
"""

import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200

rng = np.random.default_rng(SEED)

# Curated canonical TF binding motifs (consensus / IUPAC simplified to DNA).
# Drawn from JASPAR-style consensus sequences widely cited in regulatory genomics.
# Each entry is the forward-strand consensus; reverse complement applied randomly.
MOTIFS = [
    # Insulator / architectural
    "CCGCGNGGNGGCAG",       # CTCF (M00139-like core)
    "GCCACCTGCTG",          # E-box variant (USF1/MYC)
    # Hematopoietic / lineage
    "AGATAAG",              # GATA1/2/3
    "TGASTCA",              # AP-1 (TPA response, FOS/JUN family)
    "GGGRNNYYCC",           # NFKB
    "TTGTTT",                # FOXA / FOXO
    "TGACGTCA",             # CREB / ATF
    # Pluripotency / development
    "ATGCAAAT",             # OCT4 (POU)
    "CATTGT",               # SOX (HMG)
    "CAGGTG",               # E-box (MYC, ASCL)
    # Basal / housekeeping
    "GGGCGG",               # SP1
    "GGGGGCG",              # KLF
    "GCCACGTGAC",           # NRF1
    "CCATATAAG",            # YY1 (DPE-like)
    "TATAAAA",              # TBP / TATA box
    "TTCCGGGAA",            # ETS / ELK1
    # Hepatic / lineage-specific
    "TGTTTGY",              # HNF1 (Y=C/T)
    "RTAAACA",              # FOXA1/HNF3
    # Neural / TEAD / signaling
    "GGAATG",               # TEAD
    "TTTCNNTTTC",           # STAT
    "GAAANNGAAA",           # IRF
    "RRRCWWGYYY",           # p53 (half-site)
    "TTTSSCGC",             # E2F
    "GCGGGGGCG",            # EGR1 / SP-related
    "CACCT",                # ZEB / SNAI
    "TGACGTCA",             # CRE (CREB) again
    "CCAAT",                # NF-Y / CCAAT box
    "TGTGGT",               # RUNX
    "TGACAG",               # MEIS / PBX
    "CTTTGT",               # LEF/TCF (Wnt)
    "GTTGCCATGGCAAC",       # RFX (X-box)
    "AGGTCANNNNAGGTCA",     # nuclear receptor DR4
    "AGAACANNNTGTTCT",      # GR/AR half-sites
    "CAGCTG",               # bHLH MyoD
]

# IUPAC random expansion table for ambiguous letters
IUPAC = {
    "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"], "N": ["A", "C", "G", "T"],
    "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"],
    "K": ["G", "T"], "M": ["A", "C"],
    "B": ["C", "G", "T"], "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"],
}

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def realize_motif(motif: str, rng) -> str:
    """Convert IUPAC consensus to concrete DNA by sampling ambiguous positions."""
    return "".join(rng.choice(IUPAC[b]) for b in motif)


def reverse_complement(seq: str) -> str:
    return "".join(COMP[b] for b in reversed(seq))


def make_backbone(rng) -> list:
    """Random 50% GC backbone — uniform over A/C/G/T."""
    return list(rng.choice(["A", "C", "G", "T"], size=SEQ_LEN))


def plant_motifs_in_seq(rng) -> str:
    seq = make_backbone(rng)
    n_motifs = int(rng.integers(1, 5))  # 1..4 motifs
    for _ in range(n_motifs):
        motif_template = MOTIFS[rng.integers(0, len(MOTIFS))]
        realized = realize_motif(motif_template, rng)
        if rng.random() < 0.5:
            realized = reverse_complement(realized)
        if len(realized) > SEQ_LEN:
            realized = realized[:SEQ_LEN]
        max_pos = SEQ_LEN - len(realized)
        pos = int(rng.integers(0, max_pos + 1))
        for i, b in enumerate(realized):
            seq[pos + i] = b
    return "".join(seq)


def make_random_seq(rng) -> str:
    return "".join(rng.choice(["A", "C", "G", "T"], size=SEQ_LEN))


def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    sequences = []
    # 80% with motifs, 20% pure random backbone
    n_motif = int(0.8 * N_SEQS)
    n_random = N_SEQS - n_motif
    for _ in range(n_motif):
        sequences.append(plant_motifs_in_seq(rng))
    for _ in range(n_random):
        sequences.append(make_random_seq(rng))
    # Shuffle so motif/random are interleaved
    idx = rng.permutation(N_SEQS)
    sequences = [sequences[i] for i in idx]

    # Validate
    assert len(sequences) == N_SEQS
    for s in sequences:
        assert len(s) == SEQ_LEN
        assert set(s).issubset({"A", "C", "G", "T"})

    with open(out_path, "w") as f:
        for s in sequences:
            f.write(s + "\n")
    print(f"Wrote {len(sequences)} sequences to {out_path}")


if __name__ == "__main__":
    main()
