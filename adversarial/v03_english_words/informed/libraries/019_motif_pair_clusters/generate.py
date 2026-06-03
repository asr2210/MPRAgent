"""
Experiment 019: Motif pair clusters at biological spacing.

Real enhancers cluster TF motifs within ~50bp windows. Test whether the
surrogate uses motif CO-OCCURRENCE / proximity as a feature.

Design: each sequence contains 2 PAIRS of motifs. Within a pair, motifs
are spaced 5-20bp apart (center-to-center). Between pairs, distance is
unconstrained (just non-overlapping). 4 motifs/seq total, all from the
289 cell-type pool.

If exp 019 > 0.43: spatial syntax matters, plateau broken.
If exp 019 ≈ 0.42: surrogate doesn't use spacing.
"""

import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
JASPAR = ROOT / "data" / "JASPAR2024_vertebrates.txt"
OUT = Path(__file__).parent / "sequences_0.txt"

N_TOTAL = 50_000
LEN = 200
N_PAIRS = 2
PAIR_GAP_MIN = 5
PAIR_GAP_MAX = 20
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


def try_insert_pair(seq: list, pfms: list, rng: np.random.Generator,
                    positions_taken: list[tuple[int, int]]) -> bool:
    """Insert two motifs spaced PAIR_GAP_MIN..PAIR_GAP_MAX apart (gap between them)."""
    for _ in range(50):
        i1 = int(rng.integers(0, len(pfms)))
        i2 = int(rng.integers(0, len(pfms)))
        m1 = pfms[i1][1]
        m2 = pfms[i2][1]
        l1, l2 = m1.shape[1], m2.shape[1]
        gap = int(rng.integers(PAIR_GAP_MIN, PAIR_GAP_MAX + 1))
        pair_len = l1 + gap + l2
        if pair_len > LEN:
            continue
        pos = int(rng.integers(0, LEN - pair_len + 1))
        a1, b1 = pos, pos + l1
        a2, b2 = pos + l1 + gap, pos + l1 + gap + l2
        overlap = any(not (b1 <= a or a1 >= b) or not (b2 <= a or a2 >= b)
                      for a, b in positions_taken)
        if overlap:
            continue
        s1 = sample_from_pfm(m1, rng)
        s2 = sample_from_pfm(m2, rng)
        for j, ch in enumerate(s1):
            seq[a1 + j] = ch
        for j, ch in enumerate(s2):
            seq[a2 + j] = ch
        positions_taken.append((a1, b1))
        positions_taken.append((a2, b2))
        return True
    return False


def main() -> None:
    rng = np.random.default_rng(SEED)
    pfms = parse_jaspar_filtered(JASPAR)
    print(f"Loaded {len(pfms)} targeted PFMs")

    seqs = []
    pairs_placed_hist = []
    for i in range(N_TOTAL):
        backbone = rng.integers(0, 4, size=LEN)
        seq = list(BASES[backbone])
        positions_taken: list[tuple[int, int]] = []
        placed = 0
        for _ in range(N_PAIRS):
            if try_insert_pair(seq, pfms, rng, positions_taken):
                placed += 1
        pairs_placed_hist.append(placed)

        seq_str = "".join(seq)
        assert len(seq_str) == LEN
        assert set(seq_str).issubset(set("ACGT"))
        seqs.append(seq_str)

    pairs_placed_hist = np.array(pairs_placed_hist)
    print(f"Pairs placed: mean={pairs_placed_hist.mean():.2f}, "
          f"counts={[int((pairs_placed_hist == k).sum()) for k in range(N_PAIRS + 1)]}")

    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {len(seqs)} to {OUT}")


if __name__ == "__main__":
    main()
