"""Exp 011 — motif-content-aware DHS selection.

Hypothesis (H4 prediction): DHSs containing many DISTINCT TF motifs ("regulatory
grammar dense") are richer training material per-sequence than motif-sparse DHSs.

Method:
1. Subsample 200k DHSs from the non-label-aligned pool (12 NMF topics, same as
   exp 005).
2. Take a stride sample of ~80 JASPAR motifs (one every ~11) for fast scan.
3. Convert each sequence to one-hot, scan each motif at all positions with
   log-odds PWM, record per-motif MAX hit score.
4. Score each DHS = count of motifs whose best hit exceeds 80% of motif's
   max possible score (= "distinct strong motifs present").
5. Sample top 50k DHSs by this diversity score.

If eval_01 lifts above 0.6752 (exp 005) → diversity of regulatory programs WITHIN
DHS is a meaningful selection axis. If flat or worse → the model already extracts
enough grammar from natural DHS variance.
"""
import re
import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta

ROOT = Path(__file__).parent
DATA = ROOT.parent.parent / "data"
OUT = ROOT / "sequences_0.txt"

N, L = 50_000, 200
HALF = L // 2
SEED = 11
PRESCAN = 200_000  # subsample before scoring
MOTIF_STRIDE = 11  # ~80 motifs from 879
THRESH_FRAC = 0.80  # call a hit when score >= 80% of motif's max

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

# ------------------ load DHS pool ------------------
print("loading DHS index...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component"])
df = df[~df["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)
print(f"  pool size after edge filter: {len(df):,}")

idx = rng.choice(len(df), size=PRESCAN, replace=False)
sub = df.iloc[idx].reset_index(drop=True)
print(f"  prescan subsample: {len(sub):,}")

# ------------------ load sequences ------------------
print("extracting sequences...")
seqs = []
keep_idx = []
for i, r in enumerate(sub.itertuples(index=False)):
    seq = str(fa[r.seqname][r.win_start:r.win_end]).upper()
    if "N" in seq or len(seq) != L:
        continue
    seqs.append(seq)
    keep_idx.append(i)
print(f"  N-free sequences: {len(seqs):,}")

# ------------------ one-hot encode ------------------
def onehot(seqs):
    """list[str] -> (N, 4, L) float32"""
    base_to_i = np.array([4]*256, dtype=np.int8)
    base_to_i[ord('A')] = 0
    base_to_i[ord('C')] = 1
    base_to_i[ord('G')] = 2
    base_to_i[ord('T')] = 3
    arr = np.frombuffer(("".join(seqs)).encode(), dtype=np.uint8)
    idx = base_to_i[arr].reshape(len(seqs), L)  # (N, L)
    oh = np.zeros((len(seqs), 4, L), dtype=np.float32)
    rows = np.arange(len(seqs))[:, None]
    cols = np.arange(L)[None, :]
    valid = idx < 4
    oh[rows, idx, cols] = valid.astype(np.float32)
    return oh

# Reverse-complement onehot — flip channels (A<->T, C<->G) and reverse position axis
def rc_onehot(oh):
    return oh[:, [3, 2, 1, 0], ::-1].copy()

# ------------------ load JASPAR motifs ------------------
def parse_jaspar(path):
    motifs = []
    rows = {}
    header = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header and len(rows) == 4:
                    counts = np.stack([rows[b] for b in 'ACGT']).astype(float) + 0.01
                    pwm = counts / counts.sum(axis=0, keepdims=True)
                    motifs.append((header[0], header[1], pwm))
                parts = line[1:].split('\t')
                header = (parts[0], parts[1] if len(parts) > 1 else parts[0])
                rows = {}
            else:
                m = re.match(r'([ACGT])\s*\[(.*)\]', line)
                if not m:
                    continue
                rows[m.group(1)] = np.array([float(x) for x in m.group(2).split()])
        if header and len(rows) == 4:
            counts = np.stack([rows[b] for b in 'ACGT']).astype(float) + 0.01
            pwm = counts / counts.sum(axis=0, keepdims=True)
            motifs.append((header[0], header[1], pwm))
    return motifs

all_motifs = parse_jaspar(DATA / "JASPAR_vert.txt")
motifs = all_motifs[::MOTIF_STRIDE]
print(f"using {len(motifs)} motifs (strided from {len(all_motifs)})")

# Convert to log-odds matrices
log_odds_list = []
max_score_list = []
for _, _, pwm in motifs:
    lo = np.log2(pwm / 0.25).astype(np.float32)  # (4, Lm)
    log_odds_list.append(lo)
    max_score_list.append(float(lo.max(axis=0).sum()))

# ------------------ scan sequences with motifs ------------------
print("scanning motifs (this may take a minute)...")
oh = onehot(seqs)              # (N, 4, L)
oh_rc = rc_onehot(oh)          # (N, 4, L) RC strand
N_seq = oh.shape[0]
hits = np.zeros(N_seq, dtype=np.int32)   # count of motifs with strong hit

BATCH = 5000
for mi, lo in enumerate(log_odds_list):
    Lm = lo.shape[1]
    thr = THRESH_FRAC * max_score_list[mi]
    n_pos = L - Lm + 1
    found = np.zeros(N_seq, dtype=bool)
    for b0 in range(0, N_seq, BATCH):
        b1 = min(b0 + BATCH, N_seq)
        sub_oh = oh[b0:b1]
        sub_rc = oh_rc[b0:b1]
        # sliding windows: (B, 4, n_pos, Lm)
        win_f = np.lib.stride_tricks.sliding_window_view(sub_oh, window_shape=Lm, axis=2)
        win_r = np.lib.stride_tricks.sliding_window_view(sub_rc, window_shape=Lm, axis=2)
        # score per position: einsum over channels (b) and motif positions (j)
        s_f = np.einsum('nbpj,bj->np', win_f, lo)
        s_r = np.einsum('nbpj,bj->np', win_r, lo)
        max_f = s_f.max(axis=1)
        max_r = s_r.max(axis=1)
        best = np.maximum(max_f, max_r)
        found[b0:b1] = best >= thr
    hits += found.astype(np.int32)
    if (mi + 1) % 10 == 0:
        print(f"  {mi+1}/{len(motifs)} motifs, mean hits so far = {hits.mean():.2f}")

print(f"final motif-hit distribution:")
print(f"  mean = {hits.mean():.2f}, max = {hits.max()}, min = {hits.min()}")
print(f"  pct >= 5: {(hits>=5).mean()*100:.1f}%")
print(f"  pct >= 10: {(hits>=10).mean()*100:.1f}%")

# ------------------ select top by motif diversity ------------------
order = np.argsort(-hits)
top_idx = order[:N]
print(f"selected top {N}, hit range [{hits[top_idx[-1]]}, {hits[top_idx[0]]}]")
print(f"  mean hits in selected: {hits[top_idx].mean():.2f}")

selected = [seqs[i] for i in top_idx]

with OUT.open("w") as f:
    for s in selected:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
