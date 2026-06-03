"""Experiment 014 — Exact replicate of 007 with SEED=2.

Purpose: measure the noise floor of the evaluation pipeline.
007 used SEED=1 → mean_r=0.3861. If 014 (same recipe, SEED=2) lands within
±0.001, the 0.397 plateau is robust. If within ±0.003, the plateau spans
the noise.

Identical recipe: 500k GC-50 random ACGT, score with Malinois,
top 50k by max(K562, HepG2, SKNSH).
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

from utils.seqlib import write_sequences, SEQLEN
from boda.model.basset import BassetBranched
from boda.common.constants import MPRA_UPSTREAM, MPRA_DOWNSTREAM

OUT = os.path.join(THIS, "sequences_0.txt")
SCORES_OUT = os.path.join(THIS, "scores.npz")
SEED = 2                          # <-- only change vs 007
N_CAND = 500_000
N_TARGET = 50_000
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


def idx_to_seqs(idx_arr):
    lookup = np.array(list("ACGT"))
    chars = lookup[idx_arr]
    return ["".join(row) for row in chars]


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    print("Loading Malinois...")
    model = load_malinois(device)

    print(f"Generating {N_CAND:,} GC-50 random candidates...", flush=True)
    all_idx = rng.integers(0, 4, size=(N_CAND, SEQLEN), dtype=np.int8)
    print(f"  Generated shape {all_idx.shape} in {time.time()-t0:.1f}s", flush=True)

    print(f"Scoring with Malinois (batch={BATCH})...", flush=True)
    preds = score_all(model, all_idx, device)
    print(f"All preds: K562 [{preds[:,0].min():.2f}, {preds[:,0].max():.2f}] mean {preds[:,0].mean():.3f}")
    print(f"           HepG2 [{preds[:,1].min():.2f}, {preds[:,1].max():.2f}] mean {preds[:,1].mean():.3f}")
    print(f"           SKNSH [{preds[:,2].min():.2f}, {preds[:,2].max():.2f}] mean {preds[:,2].mean():.3f}")

    score = preds.max(axis=1)
    sel = np.argpartition(-score, N_TARGET)[:N_TARGET]
    sel = sel[np.argsort(-score[sel])]
    print(f"  Selected {N_TARGET} sequences; min selected score = {score[sel].min():.3f}, max = {score[sel].max():.3f}")

    selected_preds = preds[sel]
    print(f"Selected pred ranges: K562 mean {selected_preds[:,0].mean():.3f}")
    print(f"                      HepG2 mean {selected_preds[:,1].mean():.3f}")
    print(f"                      SKNSH mean {selected_preds[:,2].mean():.3f}")
    selected_idx = all_idx[sel]
    gc = ((selected_idx == 1) | (selected_idx == 2)).mean(axis=1)
    print(f"Selected GC: mean={gc.mean():.3f}, std={gc.std():.3f}")

    rng.shuffle(sel)
    selected_seqs = idx_to_seqs(all_idx[sel])
    np.savez(SCORES_OUT,
              selected_idx=sel.astype(np.int32),
              selected_preds=preds[sel].astype(np.float16))
    write_sequences(selected_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
