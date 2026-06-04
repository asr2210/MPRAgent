"""
Experiment 007: gc_50 backbone + 6 JASPAR motifs per sequence (high density).

Exp 006 showed 3-motif insertion adds ~1% over random uniform. Test if doubling
motif density (6 motifs/seq) gives proportionally more, or if there's a soft
ceiling on motif coverage. Also test a slightly improved backbone: gc_50
(exactly 50% GC by base-shuffling) instead of plain random uniform. gc_50
baseline (0.4243) was slightly above random_uniform (0.4228).

Hypothesis:
- If exp 007 > exp 006 by >1%: density matters, more motifs = better
- If exp 007 ≈ exp 006: saturation at 3 motifs
- If exp 007 < exp 006: too much motif content disrupts the surrogate's
  ability to track them individually
"""

import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
JASPAR = ROOT / "data" / "JASPAR2024_vertebrates.txt"
OUT = Path(__file__).parent / "sequences_0.txt"

N_TOTAL = 50_000
LEN = 200
N_MOTIFS_PER_SEQ = 6
SEED = 0
BASES = np.array(list("ACGT"))


def parse_jaspar(path: Path) -> list[np.ndarray]:
    text = path.read_text()
    pfms = []
    blocks = re.split(r"\n>", "\n" + text.lstrip())
    for blk in blocks:
        blk = blk.strip()
        if not blk:
            continue
        lines = blk.splitlines()
        rows = {}
        for ln in lines[1:5]:
            m = re.match(r"\s*([ACGT])\s*\[\s*(.*?)\s*\]\s*$", ln)
            if not m:
                break
            base, vals = m.group(1), m.group(2)
            rows[base] = [float(x) for x in vals.split()]
        if set(rows.keys()) != set("ACGT"):
            continue
        L = len(rows["A"])
        if L < 4 or L > 30:
            continue
        mat = np.stack([rows["A"], rows["C"], rows["G"], rows["T"]], axis=0)
        mat = mat + 0.1
        mat = mat / mat.sum(axis=0, keepdims=True)
        pfms.append(mat)
    return pfms


def sample_from_pfm(pfm: np.ndarray, rng: np.random.Generator) -> str:
    L = pfm.shape[1]
    cols = [rng.choice(4, p=pfm[:, j]) for j in range(L)]
    return "".join(BASES[c] for c in cols)


def gc50_backbone(rng: np.random.Generator) -> list[str]:
    # 100 G/C + 100 A/T, then shuffle. Half within G/C pick C or G; half within A/T pick A or T.
    half = LEN // 2
    gc_bases = rng.choice(["G", "C"], size=half)
    at_bases = rng.choice(["A", "T"], size=half)
    seq = list(gc_bases) + list(at_bases)
    rng.shuffle(seq)
    return seq


def main() -> None:
    rng = np.random.default_rng(SEED)
    pfms = parse_jaspar(JASPAR)
    print(f"Parsed {len(pfms)} PFMs")

    seqs = []
    for i in range(N_TOTAL):
        seq = gc50_backbone(rng)

        positions_taken: list[tuple[int, int]] = []
        inserted = 0
        attempts = 0
        while inserted < N_MOTIFS_PER_SEQ and attempts < 100:
            attempts += 1
            pfm_idx = rng.integers(0, len(pfms))
            pfm = pfms[pfm_idx]
            mlen = pfm.shape[1]
            if mlen > LEN:
                continue
            pos = rng.integers(0, LEN - mlen + 1)
            overlap = any(not (pos + mlen <= a or pos >= b) for a, b in positions_taken)
            if overlap and attempts < 30:
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

    gcs = np.array([(s.count("G") + s.count("C")) / LEN for s in seqs[:5000]])
    avg_motif_cov = np.mean([sum(b - a for a, b in [(0, 0)]) for _ in range(1)])
    print(f"Wrote {len(seqs)} seqs. GC: mean={gcs.mean():.3f} std={gcs.std():.3f}")


if __name__ == "__main__":
    main()
