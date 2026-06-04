"""
012_dhs_all_components
======================
50k 200bp DHS summit-centered windows, evenly distributed across ALL 16
NMF tissue components (3,125 per component). Maximum tissue-type breadth.

The opposite of 011 (which focused on only the 3 cell types we have measured
and HURT). If broad regulatory coverage helps cross-eval generalization,
this should beat dhs_specific (0.049).
"""
import gzip
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
SEQ_LEN = 200
HALF = SEQ_LEN // 2
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
DHS_PATH = DATA_DIR / "dhs_index.txt.gz"
HG38_PATH = DATA_DIR / "hg38.fa"
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
COMPONENTS = [
    "Primitive / embryonic", "Neural", "Stromal B", "Lymphoid",
    "Placental / trophoblast", "Musculoskeletal", "Cancer / epithelial",
    "Myeloid / erythroid", "Organ devel. / renal", "Tissue invariant",
    "Digestive", "Renal / cancer", "Cardiac", "Pulmonary devel.",
    "Vascular / endothelial", "Stromal A",
]
PER_COMP = 50_000 // len(COMPONENTS)  # 3125
N_SEQS = PER_COMP * len(COMPONENTS)    # 50000


def main():
    rng = np.random.default_rng(SEED)
    by_comp = {c: [] for c in COMPONENTS}
    with gzip.open(DHS_PATH, "rt") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] not in CANONICAL:
                continue
            comp = parts[9]
            if comp in by_comp:
                by_comp[comp].append((parts[0], int(parts[6])))
    for c, lst in by_comp.items():
        print(f"  {c}: {len(lst)} elements", file=sys.stderr)
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    all_seqs = []
    for comp, lst in by_comp.items():
        perm = rng.permutation(len(lst))
        seqs = []
        pi = 0
        while len(seqs) < PER_COMP and pi < len(perm):
            chrom, summit = lst[perm[pi]]
            pi += 1
            L = len(fa[chrom])
            start = max(0, summit - HALF)
            end = start + SEQ_LEN
            if end > L:
                end = L
                start = end - SEQ_LEN
            if start < 0:
                continue
            s = str(fa[chrom][start:end])
            if len(s) != SEQ_LEN or set(s) - set("ACGT"):
                continue
            seqs.append(s)
        print(f"  {comp}: kept {len(seqs)}", file=sys.stderr)
        all_seqs.extend(seqs)

    assert len(all_seqs) == N_SEQS, len(all_seqs)
    rng.shuffle(all_seqs)
    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in all_seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
