"""
Experiment 028: Strong-PFMs-only library (high info content filter).

Filter the 289 cell-type-targeted PFMs to those with high mean
information content (lowest per-column entropy). Keeps the most
"crisp" / unambiguous motifs — the surrogate may train better on
sharply-defined motif patterns vs noisy ones.

Design: same as exp 020 (8 motifs/seq, overlap, uniform random
backbone), but draws from the filtered pool.

Tests T19 indirectly: if the plateau is metric-imposed (anti-correlated
cell-type optima), filtering pool quality shouldn't change much. If
mean_r jumps, pool quality is the real lever and prior experiments were
diluted by weak PFMs.
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
TOP_FRAC = 0.3  # keep top 30% by information content
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


def mean_info_content(pfm: np.ndarray) -> float:
    """Per-column mean info content (bits). 2.0 = perfect, 0 = uniform."""
    eps = 1e-12
    entropy = -np.sum(pfm * np.log2(pfm + eps), axis=0)  # per column
    ic = 2.0 - entropy
    return float(ic.mean())


def sample_from_pfm(pfm: np.ndarray, rng: np.random.Generator) -> str:
    cols = [rng.choice(4, p=pfm[:, j]) for j in range(pfm.shape[1])]
    return "".join(BASES[c] for c in cols)


def main() -> None:
    rng = np.random.default_rng(SEED)
    all_pfms = parse_jaspar_filtered(JASPAR)
    print(f"Loaded {len(all_pfms)} targeted PFMs")

    scored = [(name, pfm, mean_info_content(pfm)) for name, pfm in all_pfms]
    scored.sort(key=lambda t: -t[2])  # descending IC
    k = max(10, int(len(scored) * TOP_FRAC))
    strong = scored[:k]
    print(f"Kept top {k} by IC. IC range: {strong[0][2]:.3f} (best) -> "
          f"{strong[-1][2]:.3f} (cutoff). Dropped IC<={scored[k][2]:.3f}")
    print(f"Top 5: {[(n, round(ic, 2)) for n, _, ic in strong[:5]]}")

    pfms = [(n, p) for n, p, _ in strong]

    seqs = []
    for _ in range(N_TOTAL):
        backbone = rng.integers(0, 4, size=LEN)
        seq = list(BASES[backbone])
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
