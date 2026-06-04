"""
Experiment 006: Random uniform backbone + JASPAR motif insertions.

Tests whether regulatory grammar (real TF binding sites) added to random
sequences can push performance above the random_uniform ceiling (0.42) in
this pipeline.

Library construction:
- Each 200bp sequence: random uniform ACGT backbone
- Inject 3 motifs per sequence, at randomly chosen positions, sampled
  from JASPAR2024 vertebrate PFMs (2346 motifs).
- Each motif instance is drawn stochastically from its PFM (position
  probabilities), not the consensus.

Hypothesis: if the pipeline rewards informative regulatory content, this
will beat 0.42. If pipeline is at hard ceiling, will match or fall below.

Prediction: ~0.42 (no improvement). I expect the surrogate model to be too
shallow to exploit discrete motif insertions in random backgrounds. Useful
NEGATIVE result either way: confirms whether motifs help here.
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


def parse_jaspar(path: Path) -> list[np.ndarray]:
    """Return list of (4, L) probability matrices, one per motif."""
    text = path.read_text()
    pfms = []
    blocks = re.split(r"\n>", "\n" + text.lstrip())
    for blk in blocks:
        blk = blk.strip()
        if not blk:
            continue
        lines = blk.splitlines()
        # First line: ">ID Name" or "ID Name"
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
        # Normalize column-wise to probabilities (with pseudocount)
        mat = mat + 0.1
        mat = mat / mat.sum(axis=0, keepdims=True)
        pfms.append(mat)
    return pfms


def sample_from_pfm(pfm: np.ndarray, rng: np.random.Generator) -> str:
    """Sample one motif instance from the PFM."""
    L = pfm.shape[1]
    cols = []
    for j in range(L):
        cols.append(rng.choice(4, p=pfm[:, j]))
    return "".join(BASES[c] for c in cols)


def main() -> None:
    rng = np.random.default_rng(SEED)
    pfms = parse_jaspar(JASPAR)
    print(f"Parsed {len(pfms)} PFMs from JASPAR")

    seqs = []
    for i in range(N_TOTAL):
        # backbone
        backbone = rng.integers(0, 4, size=LEN)
        seq = list(BASES[backbone])

        # insert N_MOTIFS_PER_SEQ motifs at random positions (non-overlapping)
        positions_taken: list[tuple[int, int]] = []
        attempts = 0
        inserted = 0
        while inserted < N_MOTIFS_PER_SEQ and attempts < 50:
            attempts += 1
            pfm_idx = rng.integers(0, len(pfms))
            pfm = pfms[pfm_idx]
            mlen = pfm.shape[1]
            if mlen > LEN:
                continue
            pos = rng.integers(0, LEN - mlen + 1)
            # Allow overlap if too many attempts; soft non-overlap
            if any(not (pos + mlen <= a or pos >= b) for a, b in positions_taken):
                if attempts < 20:
                    continue  # retry
                # else allow overlap
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

    # Diagnostics
    gcs = np.array([(s.count("G") + s.count("C")) / LEN for s in seqs[:5000]])
    print(f"GC: mean={gcs.mean():.3f} std={gcs.std():.3f}")


if __name__ == "__main__":
    main()
