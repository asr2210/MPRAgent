"""
E018: DHS restricted to chr 7, 9, 13, 21, X (parallel to E008 for DHS).

E008 (Gosai test chroms) gave +0.013 over Gosai random. Does the same
chrom restriction help on DHS, or is it Gosai-specific?

- If E018 ≈ DHS random (0.318): test-chrom boost is Gosai-specific quirk
  (maybe correlated with motif content of Gosai's chrom distribution).
- If E018 boosts to ≈ 0.33: chrom restriction is source-agnostic;
  reinforces eval-leak hypothesis.
"""
import os, gzip, numpy as np
from pyfaidx import Fasta

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
DHS = f"{ROOT}/data/DHS_Index_hg38.txt.gz"
FA = f"{ROOT}/data/hg38.fa"
SEED = 1
N_SEQS = 50_000
LEN = 200
HALF = LEN // 2
TEST_CHROMS = {"chr7", "chr9", "chr13", "chr21", "chrX"}

rng = np.random.default_rng(SEED)

print("Loading DHS test chroms...", flush=True)
recs = []
with gzip.open(DHS, "rt") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[0] not in TEST_CHROMS: continue
        recs.append((p[0], int(p[6])))
print(f"DHS sites in test chroms: {len(recs)}", flush=True)

fa = Fasta(FA, sequence_always_upper=True, as_raw=True)
chroms = set(fa.keys())
BASES = "ACGT"
order = rng.permutation(len(recs))
final = []
for ix in order:
    chrom, summit = recs[ix]
    if chrom not in chroms: continue
    s = summit - HALF; e = summit + HALF
    if s < 0 or e > len(fa[chrom]): continue
    seq = fa[chrom][s:e].upper()
    if len(seq) != LEN: continue
    nc = sum(1 for b in seq if b not in BASES)
    if nc > LEN * 0.1: continue
    if nc > 0:
        seq = "".join(b if b in BASES else BASES[rng.integers(0, 4)] for b in seq)
    final.append(seq)
    if len(final) >= N_SEQS: break

assert len(final) == N_SEQS, f"Got {len(final)}"

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
