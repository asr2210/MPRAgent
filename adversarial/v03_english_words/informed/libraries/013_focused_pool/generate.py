"""
Experiment 013: Focused 50-motif pool — only canonical cell-type TFs.

Same design as exp 008 (3 motifs/seq, random uniform backbone) but with a
NARROWER motif pool. Drop ambiguous/rare/universal TFs from the 289-motif
set; keep only canonical drivers for K562, HepG2, SK-N-SH.

Hypothesis: 289 PFMs included TFs that don't materially affect the surrogate
prediction (universals like SP1, TBP, NFY, etc., plus rare variants).
Removing them concentrates training signal on the TFs that matter.
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

FOCUSED_TOKENS = [
    "GATA1", "GATA2", "KLF1", "NFE2", "TAL1", "RUNX1", "SPI1", "MYB",
    "HNF4A", "HNF1A", "FOXA1", "FOXA2", "CEBPA", "CEBPB", "ONECUT1", "RXRA",
    "NR1H4", "PPARA",
    "NEUROD1", "NEUROG2", "ASCL1", "OLIG2", "SOX2", "SOX10", "POU3F2",
    "ISL2", "MEF2C", "PHOX2B", "REST", "RFX3", "PAX6",
]


def parse_jaspar_filtered(path: Path, tokens: list[str]) -> list[tuple[str, np.ndarray]]:
    text = path.read_text()
    pfms = []
    blocks = re.split(r"\n>", "\n" + text.lstrip())
    tokens_u = [t.upper() for t in tokens]
    for blk in blocks:
        blk = blk.strip()
        if not blk:
            continue
        lines = blk.splitlines()
        m = re.match(r"(\S+)\s+(\S+)", lines[0].strip())
        if not m:
            continue
        name = m.group(2).upper()
        if not any(t in name for t in tokens_u):
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


def main() -> None:
    rng = np.random.default_rng(SEED)
    pfms = parse_jaspar_filtered(JASPAR, FOCUSED_TOKENS)
    print(f"Loaded {len(pfms)} focused PFMs")
    print(f"Names: {sorted({n for n, _ in pfms})}")

    seqs = []
    for i in range(N_TOTAL):
        backbone = rng.integers(0, 4, size=LEN)
        seq = list(BASES[backbone])
        positions_taken: list[tuple[int, int]] = []
        inserted = 0
        attempts = 0
        while inserted < N_MOTIFS_PER_SEQ and attempts < 100:
            attempts += 1
            _, pfm = pfms[int(rng.integers(0, len(pfms)))]
            mlen = pfm.shape[1]
            if mlen > LEN:
                continue
            pos = int(rng.integers(0, LEN - mlen + 1))
            overlap = any(not (pos + mlen <= a or pos >= b) for a, b in positions_taken)
            if overlap and attempts < 40:
                continue
            motif = sample_from_pfm(pfm, rng)
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
