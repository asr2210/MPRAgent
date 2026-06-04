"""
Experiment 026: Bio-realistic dinucleotide backbone + 8 motifs overlap.

Generate backbone using mammalian dinucleotide frequencies (approximate
hg38 genome-wide): low CpG (~0.01), elevated TT/AA, balanced GC ~41%.
Then insert 8 motifs (overlap) from 289-PFM cell-type pool.

Tests whether the surrogate's training data implies a preferred backbone
composition — if so, bio-realistic backbone might trigger different
features and break plateau.
"""

import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
JASPAR = ROOT / "data" / "JASPAR2024_vertebrates.txt"
OUT = Path(__file__).parent / "sequences_0.txt"

N_TOTAL = 50_000
LEN = 200
N_MOTIFS = 8
SEED = 0
BASES = np.array(list("ACGT"))

# Approximate hg38 dinucleotide frequencies (normalized rows = transition probs
# from first base to second). Source: typical mammalian genome stats.
# Indexing: rows = first base (A,C,G,T), cols = second base
DINUC_PROBS = np.array([
    # A     C     G     T
    [0.327, 0.176, 0.240, 0.257],  # from A
    [0.354, 0.252, 0.048, 0.346],  # from C (low CG)
    [0.282, 0.193, 0.247, 0.278],  # from G
    [0.213, 0.211, 0.252, 0.324],  # from T
])
BASE_FREQS = np.array([0.295, 0.205, 0.205, 0.295])  # genome avg, ~41% GC

TARGET_TOKENS = [
    "GATA1", "GATA2", "GATA3", "KLF1", "KLF4", "KLF15", "NFE2", "MAF::NFE2",
    "MYB", "TAL1", "RUNX1", "SPI1",
    "HNF4A", "HNF4G", "HNF1A", "HNF1B", "FOXA1", "FOXA2", "FOXA3",
    "CEBPA", "CEBPB", "CEBPD", "CEBPG", "ONECUT1", "ONECUT2", "ONECUT3",
    "RXRA", "NR1H4", "PPARA",
    "NEUROD1", "NEUROG1", "NEUROG2", "ASCL1", "ASCL2", "OLIG1", "OLIG2", "OLIG3",
    "SOX2", "SOX10", "SOX21", "POU3F1", "POU3F2", "POU3F3", "POU3F4",
    "ISL2", "MEF2C", "MYCN", "PHOX2A", "PHOX2B", "REST", "RFX1", "RFX3", "RFX5",
    "PAX3", "PAX6",
    "SP1", "SP2", "SP3", "CTCF", "TBP", "NFYA", "NFYB", "NFYC",
    "FOS", "JUN", "FOSL", "JUND", "BACH",
    "ELK1", "ELK4", "ETV", "GABPA",
    "E2F1", "E2F4", "E2F6",
    "MYC", "MAX", "MAZ",
    "ATF", "CREB", "USF1", "USF2",
    "YY1", "EGR1", "NRF1",
]


def parse_jaspar_filtered(path: Path) -> list[tuple[str, np.ndarray]]:
    text = path.read_text()
    pfms = []
    blocks = re.split(r"\n>", "\n" + text.lstrip())
    targets_upper = [t.upper() for t in TARGET_TOKENS]
    for blk in blocks:
        blk = blk.strip()
        if not blk:
            continue
        lines = blk.splitlines()
        m = re.match(r"(\S+)\s+(\S+)", lines[0].strip())
        if not m:
            continue
        name = m.group(2).upper()
        if not any(t in name for t in targets_upper):
            continue
        rows = {}
        for ln in lines[1:5]:
            mb = re.match(r"\s*([ACGT])\s*\[\s*(.*?)\s*\]\s*$", ln)
            if not mb:
                break
            rows[mb.group(1)] = [float(x) for x in mb.group(2).split()]
        if set(rows.keys()) != set("ACGT"):
            continue
        L = len(rows["A"])
        if L < 4 or L > 30:
            continue
        mat = np.stack([rows["A"], rows["C"], rows["G"], rows["T"]], axis=0) + 0.1
        mat = mat / mat.sum(axis=0, keepdims=True)
        pfms.append((name, mat))
    return pfms


def sample_from_pfm(pfm: np.ndarray, rng: np.random.Generator) -> str:
    cols = [rng.choice(4, p=pfm[:, j]) for j in range(pfm.shape[1])]
    return "".join(BASES[c] for c in cols)


def generate_dinuc_backbone(rng: np.random.Generator, length: int) -> list:
    """Markov-chain backbone with mammalian dinucleotide transitions."""
    out = [int(rng.choice(4, p=BASE_FREQS))]
    for _ in range(length - 1):
        nxt = int(rng.choice(4, p=DINUC_PROBS[out[-1]]))
        out.append(nxt)
    return [BASES[c] for c in out]


def main() -> None:
    rng = np.random.default_rng(SEED)
    pfms = parse_jaspar_filtered(JASPAR)
    print(f"Loaded {len(pfms)} targeted PFMs")

    seqs = []
    for i in range(N_TOTAL):
        seq = generate_dinuc_backbone(rng, LEN)
        for _ in range(N_MOTIFS):
            _, pfm = pfms[int(rng.integers(0, len(pfms)))]
            mlen = pfm.shape[1]
            if mlen > LEN:
                continue
            pos = int(rng.integers(0, LEN - mlen + 1))
            motif = sample_from_pfm(pfm, rng)
            for j, ch in enumerate(motif):
                seq[pos + j] = ch

        seq_str = "".join(seq)
        assert len(seq_str) == LEN
        assert set(seq_str).issubset(set("ACGT"))
        seqs.append(seq_str)

    gc = np.array([(s.count("G") + s.count("C")) / LEN for s in seqs[:5000]])
    print(f"Wrote {len(seqs)} to {OUT}. GC: mean={gc.mean():.3f}, std={gc.std():.3f}")
    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")


if __name__ == "__main__":
    main()
