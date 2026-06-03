"""
Experiment 024: HYBRID library — 25k seqs of exp 020 design + 25k of exp 022.

Test library-level multiplexing:
- 25k sequences: 8 motifs/seq overlap (exp 020, K562/HepG2 best)
- 25k sequences: 10 motifs/seq overlap (exp 022, SK-N-SH best)

Same 289-PFM pool throughout. Tests whether the surrogate can learn
from BOTH sub-libraries and pick up the cell-type-specific advantages
of each.

If exp 024 > 0.4284 (020's score): library multiplexing works.
If 024 ≈ 0.4283: no benefit, surrogate averages over the library.
"""

import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
JASPAR = ROOT / "data" / "JASPAR2024_vertebrates.txt"
OUT = Path(__file__).parent / "sequences_0.txt"

N_TOTAL = 50_000
LEN = 200
N_PER_DESIGN = 25_000
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


def sample_from_pfm(pfm: np.ndarray, rng: np.random.Generator) -> str:
    cols = [rng.choice(4, p=pfm[:, j]) for j in range(pfm.shape[1])]
    return "".join(BASES[c] for c in cols)


def build_seq(rng: np.random.Generator, pfms: list, n_motifs: int) -> str:
    backbone = rng.integers(0, 4, size=LEN)
    seq = list(BASES[backbone])
    for _ in range(n_motifs):
        _, pfm = pfms[int(rng.integers(0, len(pfms)))]
        mlen = pfm.shape[1]
        if mlen > LEN:
            continue
        pos = int(rng.integers(0, LEN - mlen + 1))
        motif = sample_from_pfm(pfm, rng)
        for j, ch in enumerate(motif):
            seq[pos + j] = ch
    return "".join(seq)


def main() -> None:
    rng = np.random.default_rng(SEED)
    pfms = parse_jaspar_filtered(JASPAR)
    print(f"Loaded {len(pfms)} targeted PFMs")

    seqs = []
    # 25k from exp 020 design: 8 motifs/seq overlap
    for _ in range(N_PER_DESIGN):
        seqs.append(build_seq(rng, pfms, 8))
    # 25k from exp 022 design: 10 motifs/seq overlap
    for _ in range(N_PER_DESIGN):
        seqs.append(build_seq(rng, pfms, 10))

    # Shuffle so sub-libraries are intermixed
    rng.shuffle(seqs)

    assert len(seqs) == N_TOTAL
    for s in seqs:
        assert len(s) == LEN
        assert set(s).issubset(set("ACGT"))

    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {len(seqs)} to {OUT}")


if __name__ == "__main__":
    main()
