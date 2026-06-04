"""Experiment 010 — neural-lineage TF motif implants on GC-50 scaffolds.

50k random GC-50 200bp scaffolds with 3-5 implanted JASPAR motifs drawn
from a curated neural-lineage TF set (proneural, neural-stem, neuro-
specification factors). Targets lifting SK-N-SH (stuck at 0.13-0.15
across all 9 prior libraries).

Hypothesis: cell-type-specific motif augmentation provides cell-type-
specific training signal that lifts the lagging cell type. If SK-N-SH
rises and others stay flat, motif targeting is a new lever.

Generalization justification: neural TFs are universal in neural lineage.
A library teaching neural motif → SK-N-SH activity will transfer to any
other neural cell type held out.
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)
from utils.seqlib import write_sequences, SEQLEN
from utils.jaspar import parse_jaspar, sample_motif_instance

OUT = os.path.join(THIS, "sequences_0.txt")
SEED = 1
N = 50000

# Curated neural-lineage TF names (matched against JASPAR motif names).
# Mix of proneural, neural-stem, neuro-specification, neural-crest, glial,
# neural-specific repressors. Match by name prefix (e.g., "NEUROG" matches
# NEUROG1/2/3, "POU3F" matches POU3F1/2/3/4).
NEURAL_TF_PATTERNS = [
    "NEUROG", "NEUROD", "ASCL", "OLIG",   # proneural / oligodendro
    "SOX1", "SOX2", "SOX10", "SOX3",       # neural-stem / NCC
    "POU3F", "POU4F", "POU2F",             # POU-domain neural
    "LHX2", "LHX3", "LHX4", "LHX9",        # LIM-homeo neural patterning
    "ISL1", "ISL2",                         # motor neuron
    "ZIC1", "ZIC2", "ZIC3", "ZIC4",        # cerebellar / neural-spec
    "PAX3", "PAX6", "PAX7",                # neural-spec
    "OTX1", "OTX2",                         # forebrain
    "EN1", "EN2",                           # mid/hindbrain
    "GBX2",                                 # hindbrain
    "FOXA", "FOXG1", "FOXP1", "FOXP2", "FOXP4",  # neural fox
    "NHLH1", "NHLH2",                       # neuronal HLH
    "DLX",                                  # interneuron / fore
    "REST", "NRSF",                         # neural repressor
    "RFX",                                  # neural-spec ciliogenesis
    "TFAP2",                                # neural crest
    "HAND2", "PHOX2", "GATA3",              # SK-N-SH neuroblastoma TFs
    "MEIS",                                 # general but expressed in CNS
    "TBR1", "TBR2", "EOMES",                # cortical
    "NFIA", "NFIB", "NFIC",                 # glial / neural
    "ATOH",                                 # cerebellar / inner ear
]


def select_neural_motifs(all_motifs):
    selected = []
    for mid, name, pwm in all_motifs:
        nameU = name.upper()
        for pat in NEURAL_TF_PATTERNS:
            if pat in nameU:
                selected.append((mid, name, pwm))
                break
    return selected


def gen_gc50_scaffold(rng):
    return rng.choice(list("ACGT"), size=SEQLEN)


def implant_motifs(scaffold, n_motifs, motifs, rng):
    """Implant n_motifs non-overlapping motif instances into scaffold."""
    L = SEQLEN
    used = []
    for _ in range(n_motifs):
        mid, name, pwm = motifs[rng.integers(0, len(motifs))]
        mL = pwm.shape[1]
        if mL >= L:
            continue
        # Try up to 10 random positions for non-overlap
        for _ in range(10):
            pos = int(rng.integers(0, L - mL + 1))
            if all(not (pos < u_end and pos + mL > u_start) for (u_start, u_end) in used):
                instance = sample_motif_instance(pwm, rng)
                scaffold[pos:pos + mL] = list(instance)
                used.append((pos, pos + mL))
                break
    return scaffold


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    all_motifs = parse_jaspar()
    print(f"JASPAR motifs total: {len(all_motifs)}")
    neural = select_neural_motifs(all_motifs)
    print(f"Selected neural motifs: {len(neural)}")
    print("First 20:")
    for mid, name, pwm in neural[:20]:
        print(f"  {mid}  {name}  len={pwm.shape[1]}")

    seqs = []
    while len(seqs) < N:
        scaffold = gen_gc50_scaffold(rng)
        n_motifs = int(rng.integers(3, 6))   # 3-5 motifs
        seq_arr = implant_motifs(scaffold, n_motifs, neural, rng)
        seq = "".join(seq_arr.tolist())
        if len(seq) == SEQLEN and all(c in "ACGT" for c in seq):
            seqs.append(seq)
        if len(seqs) % 10000 == 0 and len(seqs) > 0:
            elapsed = time.time() - t0
            print(f"  {len(seqs):,}/{N:,} ({elapsed:.0f}s)")

    rng.shuffle(seqs)
    gcs = [(s.count("G") + s.count("C")) / SEQLEN for s in seqs]
    print(f"Final GC: mean={np.mean(gcs):.3f}, std={np.std(gcs):.3f}")
    write_sequences(seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
