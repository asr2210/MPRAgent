"""
Experiment 010: Cell-type-CLUSTERED motif libraries.

Instead of mixing K562/HepG2/SK-N-SH/universal motifs in each sequence (exp 008),
each sequence carries motifs of a SINGLE cell type. This produces stronger
per-cell-type signatures that the surrogate may exploit.

Library composition (50000 total):
- 14000 sequences: only K562 motifs (4 per seq) — erythroid TFs
- 14000 sequences: only HepG2 motifs (4 per seq) — hepatocyte TFs
- 14000 sequences: only SK-N-SH motifs (4 per seq) — neural TFs
- 8000 sequences: random uniform (no motifs) — baseline/diversity floor

Hypothesis: per-cell-type clustering gives the surrogate cleaner training
signal for each cell type. Should boost individual cell-type correlations and
mean_r. Particularly hopeful for SK-N-SH (currently 0.06 floor) if neural-only
clusters are sufficient signal to break that floor.
"""

import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
JASPAR = ROOT / "data" / "JASPAR2024_vertebrates.txt"
OUT = Path(__file__).parent / "sequences_0.txt"

N_TOTAL = 50_000
LEN = 200
SEED = 0
BASES = np.array(list("ACGT"))

K562_TOKENS = ["GATA1", "GATA2", "GATA3", "KLF1", "KLF4", "NFE2", "MAF::NFE2",
               "MYB", "TAL1", "RUNX1", "SPI1", "LMO2"]
HEPG2_TOKENS = ["HNF4A", "HNF4G", "HNF1A", "HNF1B", "FOXA1", "FOXA2", "FOXA3",
                "CEBPA", "CEBPB", "CEBPD", "CEBPG", "ONECUT1", "ONECUT2", "ONECUT3",
                "RXRA", "NR1H4", "PPARA"]
SKNSH_TOKENS = ["NEUROD1", "NEUROG1", "NEUROG2", "ASCL1", "ASCL2", "OLIG1", "OLIG2",
                "SOX2", "SOX10", "POU3F1", "POU3F2", "POU3F3", "POU3F4",
                "ISL2", "MEF2C", "MYCN", "PHOX2A", "PHOX2B", "REST", "RFX1", "RFX3",
                "PAX3", "PAX6"]


def parse_jaspar_for_tokens(path: Path, tokens: list[str]) -> list[np.ndarray]:
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
        pfms.append(mat)
    return pfms


def sample_from_pfm(pfm: np.ndarray, rng: np.random.Generator) -> str:
    cols = [rng.choice(4, p=pfm[:, j]) for j in range(pfm.shape[1])]
    return "".join(BASES[c] for c in cols)


def build_with_motifs(rng: np.random.Generator, pfms: list[np.ndarray], n_motifs: int) -> str:
    backbone = rng.integers(0, 4, size=LEN)
    seq = list(BASES[backbone])
    positions_taken: list[tuple[int, int]] = []
    inserted = 0
    attempts = 0
    while inserted < n_motifs and attempts < 100:
        attempts += 1
        pfm = pfms[int(rng.integers(0, len(pfms)))]
        mlen = pfm.shape[1]
        if mlen > LEN:
            continue
        pos = int(rng.integers(0, LEN - mlen + 1))
        overlap = any(not (pos + mlen <= a or pos >= b) for a, b in positions_taken)
        if overlap and attempts < 30:
            continue
        motif = sample_from_pfm(pfm, rng)
        for j, ch in enumerate(motif):
            seq[pos + j] = ch
        positions_taken.append((pos, pos + mlen))
        inserted += 1
    return "".join(seq)


def main() -> None:
    rng = np.random.default_rng(SEED)
    k562_pfms = parse_jaspar_for_tokens(JASPAR, K562_TOKENS)
    hepg2_pfms = parse_jaspar_for_tokens(JASPAR, HEPG2_TOKENS)
    sknsh_pfms = parse_jaspar_for_tokens(JASPAR, SKNSH_TOKENS)
    print(f"PFMs: K562={len(k562_pfms)}, HepG2={len(hepg2_pfms)}, SK-N-SH={len(sknsh_pfms)}")

    parts = []
    n_each = 14000
    for _ in range(n_each):
        parts.append(build_with_motifs(rng, k562_pfms, 4))
    for _ in range(n_each):
        parts.append(build_with_motifs(rng, hepg2_pfms, 4))
    for _ in range(n_each):
        parts.append(build_with_motifs(rng, sknsh_pfms, 4))
    n_random = N_TOTAL - 3 * n_each
    for _ in range(n_random):
        b = rng.integers(0, 4, size=LEN)
        parts.append("".join(BASES[b]))

    # Shuffle the library so cell-type blocks aren't adjacent
    rng.shuffle(parts)

    assert len(parts) == N_TOTAL
    for s in parts:
        assert len(s) == LEN
        assert set(s).issubset(set("ACGT"))

    with open(OUT, "w") as f:
        for s in parts:
            f.write(s + "\n")
    print(f"Wrote {len(parts)} to {OUT}")


if __name__ == "__main__":
    main()
