"""Experiment 011 — Hyperactive sequence design via simulated annealing.

Run 50k parallel SA trajectories. Each starts from a random GC-50 200bp
scaffold; at each step, propose a single-base substitution and accept it
with probability  exp(delta / T) where delta = new_Malinois_score -
old_Malinois_score, score = max(K562, HepG2, SKNSH).

This produces sequences with predicted activity well beyond what random
sampling can find (~10× more extreme tail).

Hypothesis: hyperactive designed sequences contain dense TF binding-site
clusters in proper spatial arrangements. Training on these should teach
the model the GRAMMAR of activation, which should transfer to held-out
cell types since the TF binding code is universal.

Risk: if Malinois has cell-type-specific biases or saturates, designed
sequences may overfit Malinois failure modes, hurting downstream.
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
N = 50_000
N_STEPS = 200
BATCH = 4096                 # parallel SA trajectories per batch
TEMP_START = 1.0
TEMP_END = 0.05
ART_DIR = "/data/users/arao/.private/mpra_exp/data/malinois_artifacts/artifacts"


def load_malinois(device):
    ckpt = torch.load(os.path.join(ART_DIR, "torch_checkpoint.pt"),
                       map_location="cpu", weights_only=False)
    hp = vars(ckpt["model_hparams"])
    model = BassetBranched(**hp)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    return model.to(device)


def _flank_indices():
    m = {"A": 0, "C": 1, "G": 2, "T": 3}
    up = torch.tensor([m[c] for c in MPRA_UPSTREAM[-200:]], dtype=torch.long)
    dn = torch.tensor([m[c] for c in MPRA_DOWNSTREAM[:200]], dtype=torch.long)
    return up, dn


def core_to_onehot(core_idx, up, dn):
    """core_idx: (B, 200) long. up, dn: (200,) long. Returns (B, 4, 600)."""
    B = core_idx.shape[0]
    device = core_idx.device
    full = torch.empty(B, 600, dtype=torch.long, device=device)
    full[:, :200] = up
    full[:, 200:400] = core_idx
    full[:, 400:600] = dn
    x = torch.zeros(B, 4, 600, device=device)
    x.scatter_(1, full.unsqueeze(1), 1.0)
    return x


def score_batch(model, core_idx, up, dn):
    x = core_to_onehot(core_idx, up, dn)
    with torch.no_grad():
        y = model(x)              # (B, 3)
    return y.max(dim=1).values    # (B,)


def sa_batch(model, B, n_steps, rng, device, up, dn):
    """Run B parallel SA trajectories. Returns final core_idx (B, 200) and scores (B,)."""
    core = torch.from_numpy(rng.integers(0, 4, size=(B, SEQLEN), dtype=np.int8)).long().to(device)
    cur_score = score_batch(model, core, up, dn)
    # Cooling schedule
    temps = torch.linspace(TEMP_START, TEMP_END, n_steps, device=device)
    for step in range(n_steps):
        pos = torch.randint(0, SEQLEN, (B,), device=device)
        new_base = torch.randint(0, 4, (B,), device=device)
        proposed = core.clone()
        proposed[torch.arange(B, device=device), pos] = new_base
        # Score proposed
        new_score = score_batch(model, proposed, up, dn)
        delta = new_score - cur_score
        T = temps[step]
        # Accept always if delta > 0; else with prob exp(delta/T)
        accept_prob = torch.where(delta > 0,
                                    torch.ones_like(delta),
                                    torch.exp(delta / T))
        accept = torch.rand(B, device=device) < accept_prob
        core = torch.where(accept.unsqueeze(1), proposed, core)
        cur_score = torch.where(accept, new_score, cur_score)
    return core.cpu().numpy().astype(np.int8), cur_score.cpu().numpy()


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
    up, dn = _flank_indices()
    up, dn = up.to(device), dn.to(device)

    all_core = np.empty((N, SEQLEN), dtype=np.int8)
    all_scores = np.empty(N, dtype=np.float32)
    n_done = 0
    while n_done < N:
        bs = min(BATCH, N - n_done)
        print(f"  SA batch starting at {n_done:,} (bs={bs})...", flush=True)
        t_b = time.time()
        core_np, scores = sa_batch(model, bs, N_STEPS, rng, device, up, dn)
        all_core[n_done:n_done + bs] = core_np
        all_scores[n_done:n_done + bs] = scores
        n_done += bs
        print(f"    done in {time.time()-t_b:.1f}s, mean final score={scores.mean():.2f}, max={scores.max():.2f}")

    print(f"All SA scores: min={all_scores.min():.2f}, mean={all_scores.mean():.2f}, max={all_scores.max():.2f}")

    # GC stats
    gc = ((all_core == 1) | (all_core == 2)).mean(axis=1)
    print(f"GC: mean={gc.mean():.3f}, std={gc.std():.3f}")

    # Filter GC ∈ [0.40, 0.60]
    keep = (gc >= 0.40) & (gc <= 0.60)
    print(f"GC filter keeps {keep.sum():,} of {N:,}")
    kept_core = all_core[keep]
    kept_scores = all_scores[keep]

    # If we lost too many, top off with extra rounds
    if len(kept_core) < N:
        deficit = N - len(kept_core)
        print(f"Need {deficit:,} more sequences — running additional SA rounds.")
        extra_core = []
        extra_scores = []
        while sum(len(x) for x in extra_core) < deficit:
            bs = min(BATCH, deficit - sum(len(x) for x in extra_core))
            print(f"  extra SA batch (bs={bs})...", flush=True)
            core_np, scores = sa_batch(model, bs, N_STEPS, rng, device, up, dn)
            gc_e = ((core_np == 1) | (core_np == 2)).mean(axis=1)
            keep_e = (gc_e >= 0.40) & (gc_e <= 0.60)
            extra_core.append(core_np[keep_e])
            extra_scores.append(scores[keep_e])
        extra_core = np.concatenate(extra_core)[:deficit]
        extra_scores = np.concatenate(extra_scores)[:deficit]
        kept_core = np.concatenate([kept_core, extra_core])
        kept_scores = np.concatenate([kept_scores, extra_scores])

    kept_core = kept_core[:N]
    kept_scores = kept_scores[:N]
    print(f"Final selected scores: min={kept_scores.min():.2f}, mean={kept_scores.mean():.2f}, max={kept_scores.max():.2f}")
    gc = ((kept_core == 1) | (kept_core == 2)).mean(axis=1)
    print(f"Final GC: mean={gc.mean():.3f}, std={gc.std():.3f}")

    # Shuffle and write
    order = rng.permutation(len(kept_core))
    kept_core = kept_core[order]
    kept_scores = kept_scores[order]
    seqs = idx_to_seqs(kept_core)
    np.savez(SCORES_OUT, scores=kept_scores.astype(np.float16))
    write_sequences(seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
