#!/usr/bin/env python3
"""
Experiment 003 — JASPAR motif-embedded synthetic.

50,000 sequences = 200bp random ACGT backgrounds with 4-8 JASPAR 2024 CORE
vertebrate non-redundant TF motifs embedded at random positions, sampled
uniformly from the 879 motifs, random strand.

For each motif instance we sample bases from the PWM (column-wise probability),
so each instance is a plausible match to the PWM rather than the consensus.
This gives the model a wide range of motif strengths/variants.
"""
import os
import re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
JASPAR = os.path.abspath(os.path.join(HERE, "..", "..", "data", "JASPAR2024_vertebrate.txt"))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
SEQ_LEN = 200
N_MOTIFS_PER_SEQ_LO = 4
N_MOTIFS_PER_SEQ_HI = 8  # inclusive
ALPHABET = np.array(list("ACGT"))
BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def parse_jaspar(path):
    """Parse a JASPAR PFM file, return list of (name, pwm 4×L) tuples."""
    pwms = []
    with open(path) as fh:
        cur_name = None
        cur_rows = {}
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_name is not None and len(cur_rows) == 4:
                    pwms.append((cur_name, _rows_to_pwm(cur_rows)))
                cur_name = line[1:].split()[0]
                cur_rows = {}
            else:
                m = re.match(r"^([ACGT])\s*\[\s*(.*?)\s*\]\s*$", line)
                if not m:
                    continue
                base = m.group(1)
                counts = [int(x) for x in m.group(2).split()]
                cur_rows[base] = counts
        if cur_name is not None and len(cur_rows) == 4:
            pwms.append((cur_name, _rows_to_pwm(cur_rows)))
    return pwms


def _rows_to_pwm(rows):
    L = len(rows["A"])
    mat = np.zeros((4, L), dtype=np.float64)
    for i, b in enumerate("ACGT"):
        mat[i] = rows[b]
    # Per-column normalize with pseudocount.
    mat = mat + 0.5
    mat = mat / mat.sum(axis=0, keepdims=True)
    return mat


def sample_motif_seq(pwm, rng, rc=False):
    """Sample a sequence from a PWM (per-column draw)."""
    L = pwm.shape[1]
    out = np.empty(L, dtype="<U1")
    for c in range(L):
        out[c] = ALPHABET[rng.choice(4, p=pwm[:, c])]
    s = "".join(out)
    if rc:
        s = "".join(COMPLEMENT[b] for b in s[::-1])
    return s


def main():
    rng = np.random.default_rng(SEED)
    pwms = parse_jaspar(JASPAR)
    print(f"loaded {len(pwms)} JASPAR PWMs")

    # Filter to motifs that can fit (we cap at 30bp for safety; longer motifs are rare in CORE).
    pwms = [(n, p) for n, p in pwms if 4 <= p.shape[1] <= 30]
    print(f"after length filter: {len(pwms)} PWMs")

    seqs = []
    n_motif_idx = np.arange(len(pwms))
    for i in range(N_TOTAL):
        # Random background.
        bg = ALPHABET[rng.choice(4, size=SEQ_LEN)]
        # Number of motifs for this sequence.
        n_emb = rng.integers(N_MOTIFS_PER_SEQ_LO, N_MOTIFS_PER_SEQ_HI + 1)
        # Pick motifs (without replacement to keep variety high within sequence).
        motif_ids = rng.choice(len(pwms), size=n_emb, replace=False)
        # Pick non-overlapping positions greedily.
        placed = []  # (start, end)
        for mid in motif_ids:
            pwm = pwms[mid][1]
            L = pwm.shape[1]
            # Try a few positions to find a non-overlapping one.
            for _try in range(20):
                start = int(rng.integers(0, SEQ_LEN - L + 1))
                end = start + L
                clash = any(not (end <= s or start >= e) for s, e in placed)
                if not clash:
                    rc = bool(rng.integers(0, 2))
                    motif_seq = sample_motif_seq(pwm, rng, rc=rc)
                    bg[start:end] = list(motif_seq)
                    placed.append((start, end))
                    break
        seqs.append("".join(bg))

    assert len(seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"

    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
