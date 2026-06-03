"""
Experiment 004: DHS Neural component only.

Diagnostic: does massively over-representing neural sequences move SK-N-SH r
off ~0? Meuleman 2020's "Neural" component has 461k DHS elements that load
strongly on the neural NMF axis.

If SK-N-SH r becomes nonzero → cell-type-targeted regulatory sequences are
the lever. If still ~0 → the pipeline cannot learn SK-N-SH regardless.
"""
import gzip
import os
import numpy as np
import twobitreader

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DHS = os.path.join(ROOT, "data", "DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz")
TWOBIT = os.path.join(ROOT, "data", "hg38.2bit")
OUT = os.path.join(HERE, "sequences_0.txt")

WINDOW = 200
N_SEQ = 50_000
SEED = 1
CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

rows = []
with gzip.open(DHS, "rt") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[0] in CHROMS and p[9] == "Neural":
            try:
                rows.append((p[0], int(p[6])))
            except ValueError:
                pass
print(f"Neural DHS rows: {len(rows):,}")

rng = np.random.default_rng(SEED)
idxs = rng.choice(len(rows), size=int(N_SEQ * 1.3), replace=False)
genome = twobitreader.TwoBitFile(TWOBIT)
out = []
half = WINDOW // 2
for i in idxs:
    chrom, summit = rows[i]
    start = max(0, summit - half)
    seq = genome[chrom][start:start + WINDOW].upper()
    if len(seq) == WINDOW and all(c in "ACGT" for c in seq):
        out.append(seq)
out = out[:N_SEQ]
assert len(out) == N_SEQ
with open(OUT, "w") as f:
    for s in out:
        f.write(s + "\n")
print(f"Wrote {N_SEQ}")
