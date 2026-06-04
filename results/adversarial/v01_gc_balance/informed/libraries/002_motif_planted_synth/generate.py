"""Exp 002 — Motif-planted synthetic library.

Theory test: H0 says regulatory grammar density drives library utility for
generalisation. If pure synthetic with JASPAR motifs planted at random
positions substantially outperforms pure random (synth_oracle baseline = 0.684),
then motif identity is what the model is learning. If it barely improves,
context/arrangement matters more.

Design:
  - 50,000 sequences of 200bp
  - Each sequence: uniform-random backbone, then plant 3-8 motifs at random
    non-overlapping positions, motifs drawn uniformly from JASPAR vert 2024.
  - Motif instances sampled stochastically from each PWM (not consensus).
  - Random strand orientation per motif (revcomp 50% of the time) so we
    cover both half-sites.
"""
import re
import numpy as np
from pathlib import Path

ROOT = Path(__file__).parent
DATA = ROOT.parent.parent / "data"
OUT = ROOT / "sequences_0.txt"
JASPAR = DATA / "JASPAR_vert.txt"

N, L = 50_000, 200
SEED = 2
MIN_MOTIFS, MAX_MOTIFS = 3, 8  # inclusive

rng = np.random.default_rng(SEED)


def parse_jaspar(path):
    motifs = []
    rows, header = {}, None
    def flush():
        if header and len(rows) == 4:
            counts = np.stack([rows[b] for b in 'ACGT']).astype(float) + 0.01
            pwm = counts / counts.sum(axis=0, keepdims=True)
            motifs.append((header[0], header[1], pwm))
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                flush()
                parts = line[1:].split('\t')
                header = (parts[0], parts[1] if len(parts) > 1 else parts[0])
                rows = {}
            else:
                m = re.match(r'([ACGT])\s*\[(.*)\]', line)
                if m:
                    rows[m.group(1)] = np.array([float(x) for x in m.group(2).split()])
        flush()
    return motifs


def sample_from_pwm(pwm, rng):
    cdf = pwm.cumsum(axis=0)
    u = rng.random(pwm.shape[1])
    idx = (u[None, :] >= cdf).sum(axis=0)
    return "".join("ACGT"[i] for i in idx)


COMP = str.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s.translate(COMP)[::-1]


def main():
    motifs = parse_jaspar(JASPAR)
    print(f"loaded {len(motifs)} JASPAR motifs")
    lens = np.array([m[2].shape[1] for m in motifs])
    print(f"motif length: min={lens.min()}, median={int(np.median(lens))}, max={lens.max()}")

    bases = np.array(list("ACGT"))
    backbone = bases[rng.integers(0, 4, size=(N, L))]  # (N, L) chars
    seqs = ["".join(row) for row in backbone]

    n_motifs_per_seq = rng.integers(MIN_MOTIFS, MAX_MOTIFS + 1, size=N)
    motif_choices = rng.integers(0, len(motifs), size=int(n_motifs_per_seq.sum()))
    strands = rng.random(int(n_motifs_per_seq.sum())) < 0.5

    cursor = 0
    out_seqs = []
    for i in range(N):
        seq = list(seqs[i])
        k = n_motifs_per_seq[i]
        placed_intervals = []  # list of (start, end)
        attempts = 0
        for j in range(k):
            mi = motif_choices[cursor]
            strand = strands[cursor]
            cursor += 1
            _, _, pwm = motifs[mi]
            ml = pwm.shape[1]
            if ml >= L:
                continue
            # try up to a few times to find a non-overlapping position
            for _try in range(6):
                pos = rng.integers(0, L - ml + 1)
                end = pos + ml
                if any(not (end <= s or pos >= e) for s, e in placed_intervals):
                    continue
                inst = sample_from_pwm(pwm, rng)
                if strand:
                    inst = revcomp(inst)
                for offset, ch in enumerate(inst):
                    seq[pos + offset] = ch
                placed_intervals.append((pos, end))
                break
        out_seqs.append("".join(seq))

    with OUT.open('w') as f:
        for s in out_seqs:
            assert len(s) == L
            f.write(s + "\n")
    print(f"wrote {N} sequences x {L}bp to {OUT}")


if __name__ == "__main__":
    main()
