"""
Experiment 002 — dhs_random reproduction.

Sample 50,000 sequences uniformly at random from the Meuleman 2020 DHS Index
(~3.6M elements). Extract 200bp around each summit.

Rationale: dhs_random in instructions.md gets eval_01 = 0.7089. If v08 reproduces
this number, then v08 behaves like the instructions.md baselines and my cCRE
methodology was the issue. If v08 gives ~0 here, then v08 is genuinely different
and the instructions.md baselines do not transfer.

This is a calibration probe.
"""
import os, gzip, random
import pyfaidx

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
SEQ_LEN = 200
TARGET_N = 50_000

# Read DHS index. Columns: seqname start end identifier mean_signal numsamples summit
#                           core_start core_end component
dhs = []
with gzip.open(os.path.join(DATA, "DHS_Index_hg38.txt.gz"), "rt") as f:
    header = f.readline()
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom = parts[0]
        summit = int(parts[6])
        if "_" in chrom or chrom in ("chrM", "chrEBV"):
            continue
        dhs.append((chrom, summit))
print(f"loaded {len(dhs)} DHS elements")

random.shuffle(dhs)
fa = pyfaidx.Fasta(os.path.join(DATA, "hg38.fa"))

def get_seq(chrom, summit):
    s = summit - SEQ_LEN // 2
    e = s + SEQ_LEN
    if s < 0 or e > len(fa[chrom]):
        return None
    seq = str(fa[chrom][s:e]).upper()
    if len(seq) != SEQ_LEN or "N" in seq or not set(seq).issubset(set("ACGT")):
        return None
    return seq

seqs = []
i = 0
while len(seqs) < TARGET_N and i < len(dhs):
    chrom, summit = dhs[i]
    i += 1
    seq = get_seq(chrom, summit)
    if seq is not None:
        seqs.append(seq)
print(f"extracted {len(seqs)} sequences from {i} attempts")
assert len(seqs) == TARGET_N

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
