"""
Experiment 018: Per-sequence cell-type-balanced motifs.

Each sequence contains EXACTLY:
  - 1 K562-pool motif (erythroid TFs)
  - 1 HepG2-pool motif (hepatocyte TFs)
  - 1 SK-N-SH-pool motif (neural TFs)

Total: 3 motifs/seq, but with guaranteed coverage of all 3 cell types per
training sample. Random uniform backbone.

Tests T13: forcing per-sequence cell-type balance gives each oracle a
clear signal in every training example. Compare to 008 where motifs were
uniformly sampled from a merged pool (so a given sequence might have
mostly K562-relevant motifs, mostly HepG2-relevant, etc.).
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

K562_TOKENS = ["GATA1", "GATA2", "GATA3", "KLF1", "KLF4", "KLF15", "NFE2",
               "MAF::NFE2", "MYB", "TAL1", "RUNX1", "SPI1", "LMO2", "FLI1"]
HEPG2_TOKENS = ["HNF4A", "HNF4G", "HNF1A", "HNF1B", "FOXA1", "FOXA2", "FOXA3",
                "CEBPA", "CEBPB", "CEBPD", "CEBPG", "ONECUT1", "ONECUT2",
                "ONECUT3", "RXRA", "NR1H4", "PPARA"]
SKNSH_TOKENS = ["NEUROD1", "NEUROG1", "NEUROG2", "ASCL1", "ASCL2", "OLIG1",
                "OLIG2", "OLIG3", "SOX2", "SOX10", "SOX21", "POU3F1", "POU3F2",
                "POU3F3", "POU3F4", "ISL2", "MEF2C", "MYCN", "PHOX2A",
                "PHOX2B", "REST", "RFX1", "RFX3", "RFX5", "PAX3", "PAX6"]


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


def insert_motif(seq: list, pfm: np.ndarray, rng: np.random.Generator,
                 positions_taken: list[tuple[int, int]]) -> bool:
    mlen = pfm.shape[1]
    if mlen > LEN:
        return False
    for _ in range(50):
        pos = int(rng.integers(0, LEN - mlen + 1))
        overlap = any(not (pos + mlen <= a or pos >= b) for a, b in positions_taken)
        if overlap:
            continue
        motif = sample_from_pfm(pfm, rng)
        for j, ch in enumerate(motif):
            seq[pos + j] = ch
        positions_taken.append((pos, pos + mlen))
        return True
    # fallback: allow overlap
    pos = int(rng.integers(0, LEN - mlen + 1))
    motif = sample_from_pfm(pfm, rng)
    for j, ch in enumerate(motif):
        seq[pos + j] = ch
    positions_taken.append((pos, pos + mlen))
    return True


def main() -> None:
    rng = np.random.default_rng(SEED)
    k562_pfms = parse_jaspar_for_tokens(JASPAR, K562_TOKENS)
    hepg2_pfms = parse_jaspar_for_tokens(JASPAR, HEPG2_TOKENS)
    sknsh_pfms = parse_jaspar_for_tokens(JASPAR, SKNSH_TOKENS)
    print(f"PFMs: K562={len(k562_pfms)}, HepG2={len(hepg2_pfms)}, SK-N-SH={len(sknsh_pfms)}")

    seqs = []
    for i in range(N_TOTAL):
        backbone = rng.integers(0, 4, size=LEN)
        seq = list(BASES[backbone])
        positions_taken: list[tuple[int, int]] = []
        # Insert one from each cell-type pool
        k_pfm = k562_pfms[int(rng.integers(0, len(k562_pfms)))]
        h_pfm = hepg2_pfms[int(rng.integers(0, len(hepg2_pfms)))]
        s_pfm = sknsh_pfms[int(rng.integers(0, len(sknsh_pfms)))]
        insert_motif(seq, k_pfm, rng, positions_taken)
        insert_motif(seq, h_pfm, rng, positions_taken)
        insert_motif(seq, s_pfm, rng, positions_taken)

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
