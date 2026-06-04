"""
Experiment 025: Per-cell-type DENSE clustered library.

15k seqs with 8 K562-pool motifs (overlap allowed)
15k seqs with 8 HepG2-pool motifs (overlap allowed)
15k seqs with 8 SK-N-SH-pool motifs (overlap allowed)
5k random uniform (baseline diversity)

Tests T16 (per-cell-type density preference) with full cell-type
dedication per sequence at sweet-spot density. Compare to:
- Exp 010 (cell-type clustered at density 3): mean_r 0.4217
- Exp 020 (mixed at density 8): mean_r 0.4284

If exp 025 > 0.4284: per-cell-type dedicated dense beats mixed dense.
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
N_PER_CT = 15_000
N_RANDOM = N_TOTAL - 3 * N_PER_CT
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


def build_dense_seq(rng: np.random.Generator, pfms: list, n_motifs: int) -> str:
    backbone = rng.integers(0, 4, size=LEN)
    seq = list(BASES[backbone])
    for _ in range(n_motifs):
        pfm = pfms[int(rng.integers(0, len(pfms)))]
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
    k562 = parse_jaspar_for_tokens(JASPAR, K562_TOKENS)
    hepg2 = parse_jaspar_for_tokens(JASPAR, HEPG2_TOKENS)
    sknsh = parse_jaspar_for_tokens(JASPAR, SKNSH_TOKENS)
    print(f"PFMs: K562={len(k562)}, HepG2={len(hepg2)}, SK-N-SH={len(sknsh)}")

    seqs = []
    for _ in range(N_PER_CT):
        seqs.append(build_dense_seq(rng, k562, N_MOTIFS))
    for _ in range(N_PER_CT):
        seqs.append(build_dense_seq(rng, hepg2, N_MOTIFS))
    for _ in range(N_PER_CT):
        seqs.append(build_dense_seq(rng, sknsh, N_MOTIFS))
    for _ in range(N_RANDOM):
        b = rng.integers(0, 4, size=LEN)
        seqs.append("".join(BASES[b]))

    rng.shuffle(seqs)
    assert len(seqs) == N_TOTAL

    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {len(seqs)} to {OUT}")


if __name__ == "__main__":
    main()
