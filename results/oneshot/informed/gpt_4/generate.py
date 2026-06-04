#!/usr/bin/env python3
"""Generate a deterministic 50k x 200 bp MPRA training library."""

from __future__ import annotations

import gzip
import hashlib
import math
import os
import random
import struct
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "library" / "sequences.txt"
N_TOTAL = 50_000
LENGTH = 200
SEED = 20260527
CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
BASES = "ACGT"
TWOBIT_BASES = "TCAG"


class TwoBitGenome:
    def __init__(self, path: Path):
        self.path = path
        self.handle = path.open("rb")
        sig = self.handle.read(4)
        if sig == bytes.fromhex("1A412743"):
            self.endian = ">"
        elif sig == bytes.fromhex("4327411A"):
            self.endian = "<"
        else:
            raise RuntimeError(f"{path} is not a 2bit file")
        version, seq_count, _reserved = self._unpack("III", self.handle.read(12))
        if version != 0:
            raise RuntimeError(f"Unsupported 2bit version: {version}")
        self.index: dict[str, int] = {}
        for _ in range(seq_count):
            name_len = self.handle.read(1)[0]
            name = self.handle.read(name_len).decode()
            (offset,) = self._unpack("I", self.handle.read(4))
            self.index[name] = offset
        self.records = {name: self._read_record(offset) for name, offset in self.index.items()}

    def _unpack(self, fmt: str, data: bytes):
        return struct.unpack(self.endian + fmt, data)

    def _read_uint_array(self, count: int) -> list[int]:
        if count == 0:
            return []
        data = self.handle.read(4 * count)
        return list(self._unpack("I" * count, data))

    def _read_record(self, offset: int) -> dict:
        self.handle.seek(offset)
        (dna_size,) = self._unpack("I", self.handle.read(4))
        (n_count,) = self._unpack("I", self.handle.read(4))
        n_starts = self._read_uint_array(n_count)
        n_sizes = self._read_uint_array(n_count)
        (mask_count,) = self._unpack("I", self.handle.read(4))
        self.handle.seek(8 * mask_count, os.SEEK_CUR)
        self.handle.seek(4, os.SEEK_CUR)
        packed_offset = self.handle.tell()
        return {
            "size": dna_size,
            "n_blocks": list(zip(n_starts, n_sizes)),
            "packed_offset": packed_offset,
        }

    def chroms(self) -> dict[str, int]:
        return {name: rec["size"] for name, rec in self.records.items()}

    def sequence(self, chrom: str, start: int, end: int) -> str:
        rec = self.records[chrom]
        byte_start = start // 4
        byte_end = (end + 3) // 4
        self.handle.seek(rec["packed_offset"] + byte_start)
        packed = self.handle.read(byte_end - byte_start)
        decoded: list[str] = []
        base_pos = byte_start * 4
        for byte in packed:
            for shift in (6, 4, 2, 0):
                if start <= base_pos < end:
                    decoded.append(TWOBIT_BASES[(byte >> shift) & 3])
                base_pos += 1
        seq = decoded[: end - start]
        for n_start, n_size in rec["n_blocks"]:
            n_end = n_start + n_size
            if n_end <= start:
                continue
            if n_start >= end:
                break
            a = max(start, n_start) - start
            b = min(end, n_end) - start
            seq[a:b] = "N" * (b - a)
        return "".join(seq)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def stable_u01(text: str) -> float:
    val = int.from_bytes(hashlib.blake2b(text.encode(), digest_size=8).digest(), "little")
    return (val + 0.5) / 2**64


def read_dhs_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    with gzip.open(path, "rt") as handle:
        header = handle.readline().strip().split("\t")
        col = {name: i for i, name in enumerate(header)}
        required = ["seqname", "start", "end", "identifier", "mean_signal", "numsamples", "summit", "component"]
        missing = [name for name in required if name not in col]
        if missing:
            raise RuntimeError(f"DHS index missing columns: {missing}")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[col["seqname"]]
            if chrom not in CHROMS:
                continue
            try:
                rows.append(
                    {
                        "chrom": chrom,
                        "start": int(parts[col["start"]]),
                        "end": int(parts[col["end"]]),
                        "summit": int(parts[col["summit"]]),
                        "identifier": parts[col["identifier"]],
                        "mean_signal": float(parts[col["mean_signal"]]),
                        "numsamples": int(parts[col["numsamples"]]),
                        "component": parts[col["component"]],
                    }
                )
            except ValueError:
                continue
    return rows


def read_mixture(path: Path, n_rows: int) -> np.ndarray | None:
    if not path.exists():
        return None
    with gzip.open(path, "rb") as handle:
        mix = np.load(handle)
    if mix.shape[0] == 16 and mix.shape[1] == n_rows:
        mix = mix.T
    if mix.shape[0] != n_rows or mix.shape[1] != 16:
        return None
    return np.asarray(mix, dtype=np.float32)


def weighted_without_replacement(rng: np.random.Generator, weights: np.ndarray, count: int) -> np.ndarray:
    weights = np.asarray(weights, dtype=np.float64)
    weights[~np.isfinite(weights)] = 0.0
    weights = np.maximum(weights, 0.0)
    if weights.sum() <= 0:
        weights = np.ones_like(weights)
    # Efraimidis-Spirakis keys avoid materializing a probability vector for
    # repeated calls and give deterministic weighted sampling without replacement.
    u = rng.random(len(weights))
    keys = np.log(u) / weights
    return np.argpartition(keys, -count)[-count:]


def valid_window(tb: TwoBitGenome, chrom: str, center: int, offset: int = 0) -> str | None:
    chrom_len = tb.chroms().get(chrom)
    if chrom_len is None:
        return None
    mid = center + offset
    start = mid - LENGTH // 2
    end = start + LENGTH
    if start < 0 or end > chrom_len:
        return None
    seq = tb.sequence(chrom, start, end).upper()
    if len(seq) != LENGTH or any(base not in BASES for base in seq):
        return None
    return seq


def add_dhs_sequences(
    seqs: list[str],
    seen: set[str],
    tb: TwoBitGenome,
    rows: list[dict],
    indices: list[int] | np.ndarray,
    offsets: list[int],
    quota: int,
) -> None:
    for idx in indices:
        row = rows[int(idx)]
        for offset in offsets:
            seq = valid_window(tb, row["chrom"], row["summit"], offset)
            if seq is None:
                continue
            # Canonicalize strand for uniqueness; orientation is not annotated
            # for DHSs and this avoids spending capacity on reverse complements.
            key = min(seq, revcomp(seq))
            if key in seen:
                continue
            seen.add(key)
            seqs.append(seq)
            if len(seqs) >= quota:
                return


def dinuc_shuffle(seq: str, rng: random.Random) -> str:
    # Lightweight Markov-preserving shuffle: shuffle following bases within
    # each predecessor bucket, preserving many local composition statistics.
    trans: dict[str, list[str]] = {b: [] for b in BASES}
    for a, b in zip(seq, seq[1:]):
        trans[a].append(b)
    for values in trans.values():
        rng.shuffle(values)
    current = seq[0]
    out = [current]
    for _ in range(len(seq) - 1):
        if trans[current]:
            current = trans[current].pop()
        else:
            current = rng.choice(BASES)
        out.append(current)
    return "".join(out)


IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


MOTIFS = [
    "TGASTCA",        # AP-1 / bZIP
    "CACGTG",         # E-box
    "GGGCGG",         # GC box
    "GGAAGT",         # ETS-like
    "CCAAT",          # NF-Y
    "TATAWA",         # TATA-like
    "WGATAR",         # GATA
    "TRTTTAC",        # FOX-like
    "AATGG",          # SOX-like core
    "TGACCTTTG",      # nuclear receptor-like half-sites
    "GCGCATGCGC",     # NRF1-like GC-rich core
    "CCACCAGGGGGCGCTA",  # CTCF core fragment
]


def instantiate_iupac(pattern: str, rng: random.Random) -> str:
    return "".join(rng.choice(IUPAC[ch]) for ch in pattern)


def random_background(rng: random.Random, gc: float) -> str:
    p_gc = gc / 2.0
    weights = [0.5 - p_gc, p_gc, p_gc, 0.5 - p_gc]
    return "".join(rng.choices(BASES, weights=weights, k=LENGTH))


def add_synthetic(seqs: list[str], seen: set[str], dhs_for_shuffle: list[str], quota: int) -> None:
    rng = random.Random(SEED + 17)
    while len(seqs) < quota and dhs_for_shuffle:
        src = dhs_for_shuffle[rng.randrange(len(dhs_for_shuffle))]
        seq = dinuc_shuffle(src, rng)
        key = min(seq, revcomp(seq))
        if key not in seen and all(base in BASES for base in seq):
            seen.add(key)
            seqs.append(seq)

    while len(seqs) < quota:
        gc = min(0.75, max(0.25, rng.betavariate(5, 5)))
        chars = list(random_background(rng, gc))
        motif_count = rng.choice([1, 2, 2, 3, 4])
        used: list[tuple[int, int]] = []
        for _ in range(motif_count):
            motif = instantiate_iupac(rng.choice(MOTIFS), rng)
            if rng.random() < 0.5:
                motif = revcomp(motif)
            for _attempt in range(50):
                pos = rng.randrange(8, LENGTH - len(motif) - 8)
                span = (pos, pos + len(motif))
                if all(span[1] <= a or span[0] >= b for a, b in used):
                    chars[pos : pos + len(motif)] = motif
                    used.append(span)
                    break
        seq = "".join(chars)
        key = min(seq, revcomp(seq))
        if key not in seen:
            seen.add(key)
            seqs.append(seq)


def main() -> None:
    dhs_path = DATA / "dhs_index_hg38.txt.gz"
    mix_path = DATA / "dhs_mixture.npy.gz"
    twobit_path = DATA / "hg38.2bit"
    for path in [dhs_path, twobit_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    rows = read_dhs_rows(dhs_path)
    if not rows:
        raise RuntimeError("No DHS rows loaded")
    mix = read_mixture(mix_path, len(rows))
    tb = TwoBitGenome(twobit_path)
    rng = np.random.default_rng(SEED)

    signal = np.array([r["mean_signal"] for r in rows], dtype=np.float64)
    breadth = np.array([r["numsamples"] for r in rows], dtype=np.float64)
    base_weight = np.sqrt(np.maximum(signal, 0.001)) * np.log1p(breadth)

    seqs: list[str] = []
    seen: set[str] = set()

    if mix is not None:
        topic_strength = np.maximum(mix, 0).max(axis=1)
        topic_total = np.maximum(mix, 0).sum(axis=1)
        weights = base_weight * (1.0 + topic_strength + 0.25 * topic_total)
        primary_idx = weighted_without_replacement(rng, weights, 39_000)
        add_dhs_sequences(seqs, seen, tb, rows, primary_idx, [0], 37_500)

        for topic in range(16):
            topic_w = np.maximum(mix[:, topic], 0) * np.sqrt(np.maximum(signal, 0.001))
            idx = weighted_without_replacement(rng, topic_w, 900)
            add_dhs_sequences(seqs, seen, tb, rows, idx, [0], 42_500)
            if len(seqs) >= 42_500:
                break
    else:
        weights = base_weight * (1.0 + np.array([stable_u01(r["component"]) for r in rows]))
        primary_idx = weighted_without_replacement(rng, weights, 44_000)
        add_dhs_sequences(seqs, seen, tb, rows, primary_idx, [0], 42_500)

    # Add nearby windows from high-confidence DHSs. These keep local genomic
    # context but reduce overconcentration on exact summits.
    flank_weights = base_weight * np.sqrt(breadth + 1.0)
    flank_idx = weighted_without_replacement(rng, flank_weights, 12_000)
    add_dhs_sequences(seqs, seen, tb, rows, flank_idx, [-75, 75, -125, 125], 47_000)

    # Component balancing backstop if any topics/components were under-sampled.
    by_component: dict[str, list[int]] = defaultdict(list)
    for i, row in enumerate(rows):
        by_component[row["component"]].append(i)
    for component in sorted(by_component, key=lambda c: stable_u01(c)):
        idxs = by_component[component]
        rng.shuffle(idxs)
        add_dhs_sequences(seqs, seen, tb, rows, idxs[:500], [0], 48_000)
        if len(seqs) >= 48_000:
            break

    dhs_for_shuffle = seqs[:]
    add_synthetic(seqs, seen, dhs_for_shuffle, N_TOTAL)

    if len(seqs) != N_TOTAL:
        raise RuntimeError(f"Expected {N_TOTAL} sequences, got {len(seqs)}")
    bad = [i for i, seq in enumerate(seqs) if len(seq) != LENGTH or any(base not in BASES for base in seq)]
    if bad:
        raise RuntimeError(f"Invalid sequences at indices {bad[:5]}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
