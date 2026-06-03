"""
E025: Top 25K SNR test-chr Gosai + their 25K reverse complements = 50K.

Stacks on E022 (0.341) by adding RC augmentation for strand-invariance.
E010 showed RC aug alone is neutral, but combined with the best library
it might help (more independent label noise on same coordinates).

If E025 > 0.341: RC aug + top-SNR-test-chr genuinely compounds.
If ≈ 0.341: RC adds nothing.
If < 0.341: replication (E007 style) hurts.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_UNIQUE = 25_000
TEST_CHROMS = {"7", "9", "13", "21", "X"}
COMP = str.maketrans("ACGT", "TGCA")

def rc(s): return s.translate(COMP)[::-1]

rng = np.random.default_rng(SEED)

print("Loading...", flush=True)
rows = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if c[1] not in TEST_CHROMS: continue
        if any(b not in "ACGT" for b in c[11]): continue
        try:
            kfc, hfc, sfc = float(c[5]), float(c[6]), float(c[7])
            kse, hse, sse = float(c[8]), float(c[9]), float(c[10])
        except: continue
        snr = max(abs(kfc)/max(kse, 0.05),
                  abs(hfc)/max(hse, 0.05),
                  abs(sfc)/max(sse, 0.05))
        rows.append((c[11], snr))
print(f"Pool: {len(rows)}", flush=True)

rows.sort(key=lambda r: -r[1])
top = rows[:N_UNIQUE]
unique = [r[0] for r in top]
rc_seqs = [rc(s) for s in unique]
final = unique + rc_seqs
rng.shuffle(final)
assert len(final) == 50_000
for s in final[:3]:
    assert len(s) == 200 and all(b in "ACGT" for b in s)

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
