"""
Experiment 023: Combined design — 10 motifs/seq with overlap + 60/40 pool weighting.

Synthesis of:
- T16 (exp 022): density 10 with overlap peaks SK-N-SH at 0.061
- T13 (exp 017): broader pool (60/40 cell-type/full JASPAR) helps SK-N-SH
- T11/T15: PFM stochasticity, motif presence as the feature

Each insertion: 60% draw from 289-PFM cell-type pool, 40% from full 2344 JASPAR.
10 motifs/seq, overlap allowed.

Hypothesis: stack the two SK-N-SH-friendly levers without hurting
K562/HepG2 too much. If mean_r > 0.43 → combined design wins.
"""

import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
JASPAR = ROOT / "data" / "JASPAR2024_vertebrates.txt"
OUT = Path(__file__).parent / "sequences_0.txt"

N_TOTAL = 50_000
LEN = 200
N_MOTIFS_PER_SEQ = 10
SEED = 0
BASES = np.array(list("ACGT"))
P_CELLTYPE_DRAW = 0.6

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


def parse_jaspar_all(path: Path) -> list[tuple[str, np.ndarray]]:
    text = path.read_text()
    pfms = []
    blocks = re.split(r"\n>", "\n" + text.lstrip())
    for blk in blocks:
        blk = blk.strip()
        if not blk:
            continue
        lines = blk.splitlines()
        m = re.match(r"(\S+)\s+(\S+)", lines[0].strip())
        if not m:
            continue
        name = m.group(2).upper()
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


def main() -> None:
    rng = np.random.default_rng(SEED)
    all_pfms = parse_jaspar_all(JASPAR)
    targets_u = [t.upper() for t in TARGET_TOKENS]
    ct_indices = [i for i, (n, _) in enumerate(all_pfms)
                  if any(t in n for t in targets_u)]
    full_indices = list(range(len(all_pfms)))
    print(f"Loaded {len(all_pfms)} total PFMs, {len(ct_indices)} cell-type-targeted")

    seqs = []
    for i in range(N_TOTAL):
        backbone = rng.integers(0, 4, size=LEN)
        seq = list(BASES[backbone])
        for _ in range(N_MOTIFS_PER_SEQ):
            if rng.random() < P_CELLTYPE_DRAW:
                idx = ct_indices[int(rng.integers(0, len(ct_indices)))]
            else:
                idx = full_indices[int(rng.integers(0, len(full_indices)))]
            _, pfm = all_pfms[idx]
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

    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {len(seqs)} to {OUT}")


if __name__ == "__main__":
    main()
