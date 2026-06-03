"""
Experiment 005: random backbone with planted JASPAR TF motifs.

Hypothesis: HepG2 r is stuck at ~0.56 for random sequences because random
backbones rarely contain strong TF binding motifs. Planting 5-10 known
cell-line-relevant motifs per sequence should give the model rich, learnable
HepG2 (and possibly SK-N-SH) signal.

If HepG2 jumps above 0.6 → motifs are the lever.
If unchanged → oracle response is not motif-driven; need another approach.
"""
import os
import re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
JASPAR = os.path.join(ROOT, "data", "JASPAR2024_CORE_vertebrates.meme")
OUT = os.path.join(HERE, "sequences_0.txt")

N_SEQ = 50_000
WINDOW = 200
SEED = 1
N_MOTIFS_PER_SEQ = 8  # plant 8 motifs per 200bp sequence


def load_jaspar(path):
    motifs = {}
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("MOTIF "):
            parts = line.split()
            mid = parts[1]
            name = parts[2] if len(parts) > 2 else mid
            while i < len(lines) and not lines[i].startswith("letter-probability"):
                i += 1
            hdr = lines[i].split()
            w_idx = hdr.index("w=")
            w = int(hdr[w_idx + 1])
            rows = []
            i += 1
            for _ in range(w):
                rows.append([float(x) for x in lines[i].split()])
                i += 1
            motifs[mid] = (name, np.array(rows))
        else:
            i += 1
    return motifs


def sample_pwm(pwm, rng):
    bases = np.array(list("ACGT"))
    # Normalize each row to sum to 1 (some PFMs have tiny rounding error)
    pwm_norm = pwm / pwm.sum(axis=1, keepdims=True)
    return "".join(bases[rng.choice(4, p=p)] for p in pwm_norm)


def revcomp(s):
    comp = str.maketrans("ACGT", "TGCA")
    return s.translate(comp)[::-1]


def main():
    motifs = load_jaspar(JASPAR)
    motif_keys = list(motifs.keys())
    print(f"Loaded {len(motifs)} JASPAR motifs")

    # Use a balanced selection across all TFs (no specific cell-line bias).
    # This tests "does ANY motif diversity help" — not "do specific TFs help".
    # Filter to motifs of width 5..20 (typical TFBS lengths)
    usable = [k for k, (_, p) in motifs.items() if 5 <= p.shape[0] <= 20]
    print(f"Usable motifs (width 5-20): {len(usable)}")

    rng = np.random.default_rng(SEED)
    bases = np.array(list("ACGT"))

    with open(OUT, "w") as f:
        for _ in range(N_SEQ):
            # Random backbone
            seq = list(bases[rng.integers(0, 4, WINDOW)])

            # Sample N_MOTIFS_PER_SEQ motifs and plant them at random positions
            chosen = rng.choice(len(usable), size=N_MOTIFS_PER_SEQ, replace=False)
            for mi in chosen:
                _, pwm = motifs[usable[mi]]
                w = pwm.shape[0]
                instance = sample_pwm(pwm, rng)
                if rng.random() < 0.5:
                    instance = revcomp(instance)
                # random non-overlapping placement attempt
                pos = rng.integers(0, WINDOW - w + 1)
                for j, b in enumerate(instance):
                    seq[pos + j] = b
            f.write("".join(seq) + "\n")
    print(f"Wrote {N_SEQ} sequences")


if __name__ == "__main__":
    main()
