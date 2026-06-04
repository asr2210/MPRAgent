"""
E026: Cell-balanced top-SNR within Gosai test chroms.

E022's top-SNR (max across cells) picks sequences where ANY cell is loud.
This may oversample SKNSH-dominant. E026 picks top 16.7K per cell type
by per-cell SNR (|fc|/lfcSE), ensuring balanced K562/HepG2/SKNSH
representation.

If E026 > 0.341: cell-balance + chr + SNR stacks for additional lift.
If ≈ 0.341: per-cell SNR is similar to max-SNR; balance is irrelevant.
If < 0.341: forcing per-cell balance loses some top-quality sequences.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_TOTAL = 50_000
N_PER_CELL = N_TOTAL // 3  # 16666 + remainder
TEST_CHROMS = {"7", "9", "13", "21", "X"}

rng = np.random.default_rng(SEED)

print("Loading Gosai test chroms...", flush=True)
rows = []  # (seq, kfc/kse, hfc/hse, sfc/sse)
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
        rows.append((c[11],
                     abs(kfc)/max(kse, 0.05),
                     abs(hfc)/max(hse, 0.05),
                     abs(sfc)/max(sse, 0.05)))
print(f"Pool: {len(rows)}", flush=True)

# Sort by each cell separately, take top N per cell, dedup
def topN(sortkey, n):
    sorted_rows = sorted(rows, key=lambda r: -sortkey(r))
    return [r[0] for r in sorted_rows[:n]]

k_top = topN(lambda r: r[1], N_PER_CELL + 5000)
h_top = topN(lambda r: r[2], N_PER_CELL + 5000)
s_top = topN(lambda r: r[3], N_PER_CELL + 5000)

# Dedup while interleaving
seen = set()
final = []
i = 0
while len(final) < N_TOTAL:
    for src in (k_top, h_top, s_top):
        if i < len(src) and src[i] not in seen:
            seen.add(src[i])
            final.append(src[i])
            if len(final) >= N_TOTAL: break
    i += 1
    if i > len(k_top) and i > len(h_top) and i > len(s_top):
        break

print(f"Got {len(final)} unique (after dedup interleaving)", flush=True)
# Pad with any remaining test-chr sequences if short
if len(final) < N_TOTAL:
    remaining = [r[0] for r in rows if r[0] not in seen]
    rng.shuffle(remaining)
    for s in remaining:
        final.append(s)
        if len(final) >= N_TOTAL: break

rng.shuffle(final)
assert len(final) == N_TOTAL

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
