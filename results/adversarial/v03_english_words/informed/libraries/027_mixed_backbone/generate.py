"""
Experiment 027: MIXED BACKBONE library — hedge across backbone composition.

25k seqs with random-uniform backbone + 8 motifs overlap (favors K562/HepG2)
25k seqs with mammalian dinuc backbone + 8 motifs overlap (favors SK-N-SH)
All motifs drawn from 289-PFM cell-type pool.

Tests T18: SK-N-SH oracle responds to backbone composition statistics.
If the lever-per-cell-type trade-off averages favorably at library level,
mixed backbone gives net mean_r improvement.

Compare:
- exp 020 (all random-uniform): mean_r 0.4284, SK-N-SH 0.057
- exp 026 (all mammalian dinuc): mean_r 0.4153, SK-N-SH 0.062
- exp 027 target: mean_r > 0.428 with SK-N-SH > 0.058
"""

import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
JASPAR = ROOT / "data" / "JASPAR2024_vertebrates.txt"
OUT = Path(__file__).parent / "sequences_0.txt"

N_TOTAL = 50_000
N_DINUC = 25_000
N_RANDOM = N_TOTAL - N_DINUC
LEN = 200
N_MOTIFS = 8
SEED = 0
BASES = np.array(list("ACGT"))

DINUC_PROBS = np.array([
    [0.327, 0.176, 0.240, 0.257],
    [0.354, 0.252, 0.048, 0.346],
    [0.282, 0.193, 0.247, 0.278],
    [0.213, 0.211, 0.252, 0.324],
])
BASE_FREQS = np.array([0.295, 0.205, 0.205, 0.295])

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


def dinuc_backbone(rng: np.random.Generator, length: int) -> list:
    out = [int(rng.choice(4, p=BASE_FREQS))]
    for _ in range(length - 1):
        nxt = int(rng.choice(4, p=DINUC_PROBS[out[-1]]))
        out.append(nxt)
    return [BASES[c] for c in out]


def random_backbone(rng: np.random.Generator, length: int) -> list:
    return list(BASES[rng.integers(0, 4, size=length)])


def insert_motifs(seq: list, pfms: list, rng: np.random.Generator, n: int) -> None:
    for _ in range(n):
        _, pfm = pfms[int(rng.integers(0, len(pfms)))]
        mlen = pfm.shape[1]
        if mlen > LEN:
            continue
        pos = int(rng.integers(0, LEN - mlen + 1))
        motif = sample_from_pfm(pfm, rng)
        for j, ch in enumerate(motif):
            seq[pos + j] = ch


def main() -> None:
    rng = np.random.default_rng(SEED)
    pfms = parse_jaspar_filtered(JASPAR)
    print(f"Loaded {len(pfms)} targeted PFMs")

    seqs = []
    for _ in range(N_RANDOM):
        seq = random_backbone(rng, LEN)
        insert_motifs(seq, pfms, rng, N_MOTIFS)
        seqs.append("".join(seq))
    for _ in range(N_DINUC):
        seq = dinuc_backbone(rng, LEN)
        insert_motifs(seq, pfms, rng, N_MOTIFS)
        seqs.append("".join(seq))

    rng.shuffle(seqs)
    assert len(seqs) == N_TOTAL
    for s in seqs[:100]:
        assert len(s) == LEN
        assert set(s).issubset(set("ACGT"))

    gc = np.array([(s.count("G") + s.count("C")) / LEN for s in seqs[:5000]])
    print(f"Wrote {len(seqs)} to {OUT}. GC: mean={gc.mean():.3f}, std={gc.std():.3f}")
    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")


if __name__ == "__main__":
    main()
