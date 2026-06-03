"""
E013: Motif-engineered — Gosai backbones with planted TF motifs.

Take 50K random Gosai backbones, inject 2 strong cell-type-specific
consensus motifs into each at random positions (random orientation).
Motifs cover K562, HepG2, SKNSH primary TFs.

Hypothesis: if pipeline is capacity-limited on K562/HepG2 because their
motif signal is too sparse/subtle in Gosai, AMPLIFYING motif content
should boost per-cell ceilings. If not, the bottleneck is downstream
(model architecture, label noise on those cells).

Motifs (consensus, from literature):
- K562 erythroid: GATA1=GATAAG, KLF1=CCACGCCC
- HepG2 liver: HNF4A=AGGTCAAAGGTCA, HNF1A=GTTAATNATTAAC, FOXA1=TGTTTAC
- SKNSH neural: ASCL1=CAGCTG, NEUROD1=NCAGCTG, REST=TTCAGCACC
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000
LEN = 200

rng = np.random.default_rng(SEED)

# Cell-type-specific motifs (consensus, with N expanded to random)
MOTIFS = {
    "GATA1":   "GATAAG",
    "KLF1":    "CCACGCCC",
    "HNF4A":   "AGGTCAAAGGTCA",
    "HNF1A":   "GTTAATCATTAAC",  # N→C
    "FOXA1":   "TGTTTAC",
    "ASCL1":   "CAGCTG",
    "NEUROD1": "GCAGCTG",
    "REST":    "TTCAGCACC",
}
NAMES = list(MOTIFS.keys())
COMP = str.maketrans("ACGT", "TGCA")

def rc(s): return s.translate(COMP)[::-1]

print("Loading Gosai backbones...", flush=True)
seqs = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if any(b not in "ACGT" for b in c[11]): continue
        seqs.append(c[11])
print(f"Pool: {len(seqs)}", flush=True)

idx = rng.permutation(len(seqs))[:N_SEQS]
backbones = [seqs[i] for i in idx]

final = []
for bb in backbones:
    s = list(bb)
    # Pick 2 random motifs (with replacement OK across categories)
    picks = rng.choice(NAMES, size=2, replace=False)
    used_ranges = []
    for name in picks:
        m = MOTIFS[name]
        if rng.random() < 0.5:
            m = rc(m)
        # Pick a position that doesn't overlap previous insertions
        for _ in range(20):
            pos = int(rng.integers(0, LEN - len(m) + 1))
            if not any(pos < r2 and pos + len(m) > r1 for r1, r2 in used_ranges):
                break
        for i, b in enumerate(m):
            s[pos + i] = b
        used_ranges.append((pos, pos + len(m)))
    final.append("".join(s))

rng.shuffle(final)
assert len(final) == N_SEQS
for s in final[:5]:
    assert len(s) == LEN and all(b in "ACGT" for b in s)

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
