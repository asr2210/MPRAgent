"""
E011: Compound boosters — Gosai test chroms (7,9,13,21,X) AND loud (|fc|>=1)
AND clean (lfcSE < median in ≥2 cells).

Tests whether the small boosters from E008 (chrom filter +0.013) and
E005 (quality filter ~0) compound. If we get ~0.34+, possibly the
chrom filter does carry real signal. If we get ~0.32, the quality
filter actively neutralizes some of the chrom benefit.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000
TEST_CHROMS = {"7", "9", "13", "21", "X"}

rng = np.random.default_rng(SEED)

rows = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if c[1] not in TEST_CHROMS: continue
        if any(b not in "ACGT" for b in c[11]): continue
        try:
            rows.append((c[11], float(c[5]), float(c[6]), float(c[7]),
                         float(c[8]), float(c[9]), float(c[10])))
        except: continue
print(f"Test-chr pool: {len(rows)}", flush=True)

arr = np.array([r[1:] for r in rows])
seqs = [r[0] for r in rows]

med_kse, med_hse, med_sse = np.median(arr[:,3]), np.median(arr[:,4]), np.median(arr[:,5])

max_abs = np.maximum.reduce([np.abs(arr[:,0]), np.abs(arr[:,1]), np.abs(arr[:,2])])
se_ok = ((arr[:,3] < med_kse).astype(int) +
         (arr[:,4] < med_hse).astype(int) +
         (arr[:,5] < med_sse).astype(int)) >= 2

mask = (max_abs >= 1.0) & se_ok
n_pass = mask.sum()
print(f"Pass loud+clean: {n_pass}", flush=True)

valid = np.where(mask)[0]
if len(valid) < N_SEQS:
    print(f"Relaxing |fc| threshold")
    for thresh in [0.8, 0.6, 0.4, 0.0]:
        mask = (max_abs >= thresh) & se_ok
        valid = np.where(mask)[0]
        print(f"  thresh={thresh}: {len(valid)}")
        if len(valid) >= N_SEQS: break

assert len(valid) >= N_SEQS

sel = rng.choice(valid, N_SEQS, replace=False)
final = [seqs[i] for i in sel]
out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
