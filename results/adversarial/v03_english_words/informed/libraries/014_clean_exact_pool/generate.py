"""
Experiment 014: CLEAN exact-match pool — only the intended human TFs.

Exp 013 was polluted by substring match (MYB → 99 plant MYBs). This rebuilds
with strict TF name matching to test whether smaller-but-clean beats
larger-but-mixed (exp 008, 289 PFMs).

Specifically: human cell-type TFs + a curated set of universals that exp 008
benefited from.
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

# Exact TF names (or recognized JASPAR variants). Match if name == token or
# name starts with f"{token}(" or name contains f"::{token}" or f"{token}::"
EXACT_NAMES = {
    # K562 / erythroid
    "GATA1", "GATA2", "GATA3", "KLF1", "KLF4", "NFE2", "MYB", "TAL1",
    "RUNX1", "SPI1", "LMO2", "FLI1", "MAF::NFE2", "GATA1::TAL1",
    # HepG2 / hepatocyte
    "HNF4A", "HNF4G", "HNF1A", "HNF1B", "FOXA1", "FOXA2", "FOXA3",
    "CEBPA", "CEBPB", "CEBPD", "CEBPG", "ONECUT1", "ONECUT2", "ONECUT3",
    "RXRA", "NR1H4", "PPARA", "NR1H2::RXRA", "PPARG::RXRA",
    # SK-N-SH / neural
    "NEUROD1", "NEUROD2", "NEUROG1", "NEUROG2", "ASCL1", "ASCL2",
    "OLIG1", "OLIG2", "OLIG3", "SOX2", "SOX10", "SOX21",
    "POU3F1", "POU3F2", "POU3F3", "POU3F4", "ISL1", "ISL2",
    "MEF2C", "MEF2A", "MYCN", "PHOX2A", "PHOX2B", "REST", "RFX1", "RFX3",
    "PAX3", "PAX6", "POU2F1::SOX2", "POU5F1::SOX2",
    # Universal regulators that exp 008 included
    "SP1", "SP2", "SP3", "CTCF", "TBP", "NFYA", "NFYB", "NFYC",
    "FOS", "JUN", "FOSL1", "FOSL2", "JUND", "JUNB", "BACH1", "BACH2",
    "ELK1", "ELK4", "ETV1", "ETV2", "ETV4", "ETV5", "ETV6", "GABPA",
    "E2F1", "E2F4", "E2F6", "MYC", "MAX", "MAZ",
    "ATF1", "ATF2", "ATF3", "ATF4", "CREB1", "USF1", "USF2",
    "YY1", "EGR1", "NRF1",
}


def name_matches(name: str) -> bool:
    n = name.upper()
    if n in EXACT_NAMES:
        return True
    # Match dimer variants where one half is in our list
    if "::" in n:
        parts = n.split("::")
        return any(p in EXACT_NAMES for p in parts)
    # Match variant suffixes like GATA1(var.2) — but the JASPAR format uses
    # MA0035.4 GATA1 — name field is just GATA1 so the simple set check above
    # covers it. No extra logic needed.
    return False


def parse_jaspar_filtered(path: Path) -> list[tuple[str, np.ndarray]]:
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
        name = m.group(2)
        if not name_matches(name):
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
        pfms.append((name.upper(), mat))
    return pfms


def sample_from_pfm(pfm: np.ndarray, rng: np.random.Generator) -> str:
    cols = [rng.choice(4, p=pfm[:, j]) for j in range(pfm.shape[1])]
    return "".join(BASES[c] for c in cols)


def main() -> None:
    rng = np.random.default_rng(SEED)
    pfms = parse_jaspar_filtered(JASPAR)
    print(f"Loaded {len(pfms)} CLEAN PFMs")
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
