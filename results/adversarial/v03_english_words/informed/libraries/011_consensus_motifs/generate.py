"""
Experiment 011: Cell-type-targeted CONSENSUS motifs (no PFM stochasticity).

Same design as exp 008 (best so far at 0.4283) but replaces stochastic PFM
sampling with the consensus (argmax) string at every position. Tests whether
PFM noise was diluting the surrogate's signal — if pure consensus inserts
yield cleaner per-TF features, mean_r should rise. If PFM noise was actually
helping (more diverse motif instances → better generalization), this should
underperform 008.
"""

import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
JASPAR = ROOT / "data" / "JASPAR2024_vertebrates.txt"
OUT = Path(__file__).parent / "sequences_0.txt"

N_TOTAL = 50_000
LEN = 200
N_MOTIFS_PER_SEQ = 3
SEED = 0
BASES = np.array(list("ACGT"))

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


def consensus_from_pfm(pfm: np.ndarray) -> str:
    cols = [int(np.argmax(pfm[:, j])) for j in range(pfm.shape[1])]
    return "".join(BASES[c] for c in cols)


def main() -> None:
    rng = np.random.default_rng(SEED)
    pfms = parse_jaspar_filtered(JASPAR)
    print(f"Loaded {len(pfms)} targeted PFMs")
    # Pre-compute consensus strings (deterministic, same each call)
    consensus = [(name, consensus_from_pfm(pfm)) for name, pfm in pfms]

    seqs = []
    for i in range(N_TOTAL):
        backbone = rng.integers(0, 4, size=LEN)
        seq = list(BASES[backbone])

        positions_taken: list[tuple[int, int]] = []
        inserted = 0
        attempts = 0
        while inserted < N_MOTIFS_PER_SEQ and attempts < 100:
            attempts += 1
            idx = int(rng.integers(0, len(consensus)))
            _, motif = consensus[idx]
            mlen = len(motif)
            if mlen > LEN:
                continue
            pos = int(rng.integers(0, LEN - mlen + 1))
            overlap = any(not (pos + mlen <= a or pos >= b) for a, b in positions_taken)
            if overlap and attempts < 30:
                continue
            for j, ch in enumerate(motif):
                seq[pos + j] = ch
            positions_taken.append((pos, pos + mlen))
            inserted += 1

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
