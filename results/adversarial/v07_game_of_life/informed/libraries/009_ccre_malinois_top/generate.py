"""Experiment 009 — cCREs filtered by GC and re-ranked by Malinois activity.

All ~1M ENCODE cCREs → extract 200bp center windows → filter GC ∈ [0.45, 0.55]
→ score with Malinois → take top 50k by max(predicted log2FC).

Hypothesis: combining biological sequence statistics (real enhancers) with
oracle activity selection should beat either alone. Random+Malinois saturates
at 0.397; pure cCREs at 0.392. Their combination might break 0.40.

Generalization justification: real enhancer sequences contain natural TF
motif clusters and the spacing evolution has selected for. A library of real
high-activity enhancers exposes the model to biologically valid motif syntax,
which should transfer to held-out cell types better than random sequences
with the same predicted activity.
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np
import torch

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)
sys.path.insert(0, "/data/users/arao/.private/mpra_exp/boda2")

from utils.seqlib import extract, write_sequences, SEQLEN
from boda.model.basset import BassetBranched
from boda.common.constants import MPRA_UPSTREAM, MPRA_DOWNSTREAM

CCRE_BED = os.path.join(REPO, "data", "ENCODE_cCREs_v3.bed")
OUT = os.path.join(THIS, "sequences_0.txt")
SCORES_OUT = os.path.join(THIS, "scores.npz")
SEED = 1
N_TARGET = 50_000
GC_MIN, GC_MAX = 0.45, 0.55
BATCH = 1024
ART_DIR = "/data/users/arao/.private/mpra_exp/data/malinois_artifacts/artifacts"


def load_malinois(device):
    ckpt = torch.load(os.path.join(ART_DIR, "torch_checkpoint.pt"),
                       map_location="cpu", weights_only=False)
    hp = vars(ckpt["model_hparams"])
    model = BassetBranched(**hp)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    return model.to(device)


_UP_IDX, _DN_IDX = None, None
def get_flanks():
    global _UP_IDX, _DN_IDX
    if _UP_IDX is None:
        m = {"A": 0, "C": 1, "G": 2, "T": 3}
        _UP_IDX = np.array([m[c] for c in MPRA_UPSTREAM[-200:]], dtype=np.int8)
        _DN_IDX = np.array([m[c] for c in MPRA_DOWNSTREAM[:200]], dtype=np.int8)
    return _UP_IDX, _DN_IDX


M_LOOKUP = np.full(256, -1, dtype=np.int8)
for i, c in enumerate("ACGT"):
    M_LOOKUP[ord(c)] = i


def seqs_to_idx(seqs):
    """Convert list of ACGT strings to (N, 200) int8 array."""
    N = len(seqs)
    out = np.empty((N, SEQLEN), dtype=np.int8)
    for i, s in enumerate(seqs):
        for j, c in enumerate(s):
            out[i, j] = M_LOOKUP[ord(c)]
    return out


def onehot_batch_from_idx(core_idx, device):
    up, dn = get_flanks()
    B = core_idx.shape[0]
    full = np.empty((B, 600), dtype=np.int64)
    full[:, :200] = up[None, :]
    full[:, 200:400] = core_idx
    full[:, 400:600] = dn[None, :]
    x = torch.zeros(B, 4, 600, device=device)
    x.scatter_(1, torch.from_numpy(full).to(device).unsqueeze(1), 1.0)
    return x


def score_all(model, all_idx, device):
    N = all_idx.shape[0]
    preds = np.zeros((N, 3), dtype=np.float32)
    t0 = time.time()
    for i in range(0, N, BATCH):
        batch_idx = all_idx[i:i + BATCH]
        x = onehot_batch_from_idx(batch_idx, device)
        with torch.no_grad():
            y = model(x)
        preds[i:i + BATCH] = y.cpu().numpy()
        if (i // BATCH) % 50 == 0:
            dt = time.time() - t0
            done = i + batch_idx.shape[0]
            rate = done / max(dt, 1e-6)
            eta = (N - done) / max(rate, 1e-6)
            print(f"  scored {done:,}/{N:,} | {rate:.0f}/s | ETA {eta:.0f}s", flush=True)
    return preds


def load_ccres_windows(rng):
    """Load all cCREs, extract 200bp center windows (+ a few jittered),
    keep only ACGT-clean sequences with GC in [0.45, 0.55]."""
    print("Loading cCREs...")
    seqs = []
    rejected_n = 0
    rejected_gc = 0
    n_ccres = 0
    with open(CCRE_BED) as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) < 3:
                continue
            chrom = p[0]
            try:
                s = int(p[1]); e = int(p[2])
            except ValueError:
                continue
            n_ccres += 1
            # Center window
            center = (s + e) // 2
            # Try center + 2 jitters
            for off in (0, -50, 50, -100, 100):
                seq = extract(chrom, center + off)
                if seq is None:
                    rejected_n += 1
                    continue
                gc = (seq.count("G") + seq.count("C")) / SEQLEN
                if gc < GC_MIN or gc > GC_MAX:
                    rejected_gc += 1
                    continue
                seqs.append(seq)
            if n_ccres % 100000 == 0:
                print(f"  processed {n_ccres:,} cCREs, kept {len(seqs):,}")
    print(f"  Total cCREs: {n_ccres:,}; kept {len(seqs):,} sequences "
          f"(N rejected {rejected_n:,}, GC rejected {rejected_gc:,})")
    # De-duplicate
    seen = set()
    uniq = []
    for s in seqs:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    print(f"  After dedup: {len(uniq):,}")
    rng.shuffle(uniq)
    return uniq


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    print("Loading Malinois...")
    model = load_malinois(device)

    cand_seqs = load_ccres_windows(rng)
    if len(cand_seqs) < N_TARGET:
        raise RuntimeError(f"Only {len(cand_seqs)} candidates after GC filter — need >= {N_TARGET}")

    print(f"Converting {len(cand_seqs):,} candidates to index arrays...")
    all_idx = seqs_to_idx(cand_seqs)
    print(f"  shape {all_idx.shape}")

    print(f"Scoring with Malinois (batch={BATCH})...", flush=True)
    preds = score_all(model, all_idx, device)
    score = preds.max(axis=1)
    print(f"max-cell score: min={score.min():.2f}, median={np.median(score):.2f}, max={score.max():.2f}")

    sel = np.argpartition(-score, N_TARGET)[:N_TARGET]
    print(f"  Selected {N_TARGET}; min score = {score[sel].min():.3f}, max = {score[sel].max():.3f}")
    print(f"Selected preds:")
    sp = preds[sel]
    print(f"  K562  : mean {sp[:,0].mean():.2f}, std {sp[:,0].std():.2f}")
    print(f"  HepG2 : mean {sp[:,1].mean():.2f}, std {sp[:,1].std():.2f}")
    print(f"  SKNSH : mean {sp[:,2].mean():.2f}, std {sp[:,2].std():.2f}")

    sel_idx = all_idx[sel]
    gc = ((sel_idx == 1) | (sel_idx == 2)).mean(axis=1)
    print(f"  GC: mean={gc.mean():.3f}, std={gc.std():.3f}")

    rng.shuffle(sel)
    selected_seqs = [cand_seqs[i] for i in sel]
    np.savez(SCORES_OUT,
              selected_idx=sel.astype(np.int32),
              selected_preds=preds[sel].astype(np.float16))
    write_sequences(selected_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
