"""
Experiment 006: Dirichlet(0.5) backbones with JASPAR motif instances inserted.

Tests whether explicit TF motif content adds learning signal that the
model can use, beyond what diverse composition provides.

For each sequence:
- Generate Dirichlet(0.5) iid 200bp backbone (current best at 0.1395)
- Sample 2 motifs uniformly from JASPAR CORE 2024 vertebrates (879 motifs)
- For each motif, sample a binding instance from its PFM
- Insert at random non-overlapping positions

If exp 006 > 0.1395: model uses motif/k-mer features → push motif-rich
If exp 006 ≈ 0.1395: noise; motifs don't help → stick with pure dirichlet
If exp 006 < 0.1395: motif insertion overrides composition diversity → bad
"""
import re
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[2] / "data"
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
N_MOTIFS_PER_SEQ = 2
BASES = np.array(list("ACGT"))


def parse_meme(path):
    motifs = []
    with open(path) as f:
        text = f.read()
    blocks = text.split("MOTIF ")[1:]
    for blk in blocks:
        lines = blk.splitlines()
        header = lines[0].strip()
        name = header.split()[0]
        # find "letter-probability matrix" line
        for i, ln in enumerate(lines):
            if "letter-probability matrix" in ln:
                m = re.search(r"w=\s*(\d+)", ln)
                w = int(m.group(1))
                mat = []
                for j in range(w):
                    row = list(map(float, lines[i + 1 + j].split()))
                    mat.append(row)
                mat = np.asarray(mat)
                mat = mat / mat.sum(axis=1, keepdims=True)  # normalize
                motifs.append((name, mat))  # (w, 4)
                break
    return motifs


def main():
    rng = np.random.default_rng(SEED)
    motifs = parse_meme(DATA / "jaspar2024_vertebrates_meme.txt")
    print(f"loaded {len(motifs)} motifs; widths: "
          f"min={min(m[1].shape[0] for m in motifs)} "
          f"max={max(m[1].shape[0] for m in motifs)}")

    backbone_probs = rng.dirichlet((0.5, 0.5, 0.5, 0.5), size=N)

    seqs = []
    for i in range(N):
        # backbone
        seq = list(BASES[rng.choice(4, size=LEN, p=backbone_probs[i])])

        # insert motifs
        used = []
        for _ in range(N_MOTIFS_PER_SEQ):
            name, pwm = motifs[rng.integers(len(motifs))]
            w = pwm.shape[0]
            if w > LEN:
                continue
            # sample motif instance from PWM
            inst = [BASES[rng.choice(4, p=pwm[k])] for k in range(w)]
            # find non-overlapping position
            for _try in range(20):
                pos = rng.integers(0, LEN - w + 1)
                if not any(p <= pos < p + ww or pos <= p < pos + w
                           for p, ww in used):
                    used.append((pos, w))
                    seq[pos:pos + w] = inst
                    break
        seqs.append("".join(seq))

    assert all(len(s) == 200 for s in seqs)
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} seqs to {OUT}")


if __name__ == "__main__":
    main()
