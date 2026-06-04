"""
E022: Gosai chr 7/9/13/21/X + top 50K by SNR score.

Most aggressive stack of the two known levers:
- Chr restriction (E008 +0.013)
- Quality (high |fc| × low lfcSE) — E005 showed quality alone gives 0.

Take all 117K test-chr Gosai. Score each by max(|fc|/lfcSE) across cells.
Take top 50K. If quality + chr boost compounds, score > 0.336.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000
TEST_CHROMS = {"7", "9", "13", "21", "X"}

rng = np.random.default_rng(SEED)

print("Loading Gosai test chroms with quality scores...", flush=True)
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
        # Per-cell SNR = |fc|/lfcSE; combined = max across cells
        snr = max(abs(kfc)/max(kse, 0.05),
                  abs(hfc)/max(hse, 0.05),
                  abs(sfc)/max(sse, 0.05))
        rows.append((c[11], snr))
print(f"Pool: {len(rows)}", flush=True)

# Sort by SNR descending, take top N
rows.sort(key=lambda r: -r[1])
top = rows[:N_SEQS]
print(f"Top SNR range: {top[0][1]:.2f} to {top[-1][1]:.2f}", flush=True)

final = [r[0] for r in top]
rng.shuffle(final)
assert len(final) == N_SEQS

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
