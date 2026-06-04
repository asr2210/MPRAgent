"""Experiment 022 — Anti-oracle (bottom 50k by Malinois max-activity).

Mirror of 007 with sign flipped: take bottom 50k by max cell activity.
Sequences where ALL THREE cells predict LOW activity — i.e. inactive
across the panel.

Confirms oracle directionality. Should give mean_r < 0.378 (worse
than random) if oracle direction matters.
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
SEED = 1
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


_UP, _DN = None, None
def get_flanks():
    global _UP, _DN
    if _UP is None:
        m = {"A": 0, "C": 1, "G": 2, "T": 3}
        _UP = np.array([m[c] for c in MPRA_UPSTREAM[-200:]], dtype=np.int8)
        _DN = np.array([m[c] for c in MPRA_DOWNSTREAM[:200]], dtype=np.int8)
    return _UP, _DN


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
    for i in range(0, N, BATCH):
        b = all_idx[i:i + BATCH]
        x = onehot_batch_from_idx(b, device)
        with torch.no_grad():
            y = model(x)
        preds[i:i + BATCH] = y.cpu().numpy()
    return preds


def idx_to_seqs(idx_arr):
    lookup = np.array(list("ACGT"))
    return ["".join(row) for row in lookup[idx_arr]]


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")
    print("Loading Malinois...")
    model = load_malinois(device)

    print(f"Generating {N_CAND:,} GC-50 random candidates...")
    pool = rng.integers(0, 4, size=(N_CAND, SEQLEN), dtype=np.int8)
    print("Scoring with Malinois...")
    preds = score_all(model, pool, device)

    # ANTI-ORACLE: bottom by MAX cell activity
    score = preds.max(axis=1)
    sel = np.argpartition(score, N_TARGET)[:N_TARGET]   # smallest N
    print(f"Selected anti-oracle stats:")
    print(f"  Max-cell score: min={score[sel].min():.3f}, max={score[sel].max():.3f}, "
          f"mean={score[sel].mean():.3f}")
    sp = preds[sel]
    print(f"  K562 mean {sp[:,0].mean():.3f}, HepG2 {sp[:,1].mean():.3f}, "
          f"SKNSH {sp[:,2].mean():.3f}")
    gc = ((pool[sel] == 1) | (pool[sel] == 2)).mean(axis=1)
    print(f"  Selected GC: mean={gc.mean():.3f}, std={gc.std():.3f}")

    rng.shuffle(sel)
    seqs = idx_to_seqs(pool[sel])
    np.savez(SCORES_OUT, selected_preds=preds[sel].astype(np.float16))
    write_sequences(seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
