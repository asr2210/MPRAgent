#!/usr/bin/env python3
import bisect
import os
import random
import struct
from collections import defaultdict


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "library", "sequences.txt")
CCRE_BED = os.path.join(DATA, "GRCh38-cCREs.bed")
HG38_2BIT = os.path.join(DATA, "hg38.2bit")
SEED = 271828
LENGTH = 200
TARGET = 50_000
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


CCRE_COUNTS = {
    "dELS": 13_000,
    "pELS": 7_500,
    "PLS": 6_000,
    "CA-CTCF": 4_500,
    "CA": 4_000,
    "CA-H3K4me3": 3_000,
    "TF": 2_500,
    "CA-TF": 1_500,
}
BACKGROUND_COUNT = 3_000
SHUFFLE_COUNT = 3_000
SYNTH_COUNT = TARGET - sum(CCRE_COUNTS.values()) - BACKGROUND_COUNT - SHUFFLE_COUNT


class TwoBit:
    BASES = "TCAG"

    def __init__(self, path):
        self.handle = open(path, "rb")
        sig = self.handle.read(4)
        if sig == struct.pack(">I", 0x1A412743):
            self.endian = ">"
        elif sig == struct.pack("<I", 0x1A412743):
            self.endian = "<"
        else:
            raise ValueError("Not a 2bit file")
        version, count, _reserved = self._unpack("III")
        if version != 0:
            raise ValueError(f"Unsupported 2bit version {version}")
        self.index = {}
        self._records = {}
        for _ in range(count):
            name_len = self.handle.read(1)[0]
            name = self.handle.read(name_len).decode("ascii")
            offset = self._unpack("I")[0]
            self.index[name] = offset

    def _unpack(self, fmt):
        size = struct.calcsize(self.endian + fmt)
        return struct.unpack(self.endian + fmt, self.handle.read(size))

    def chroms(self, chrom=None):
        if chrom is not None:
            return self._record(chrom)["size"]
        return {name: self._record(name)["size"] for name in self.index}

    def _record(self, chrom):
        if chrom in self._records:
            return self._records[chrom]
        self.handle.seek(self.index[chrom])
        dna_size = self._unpack("I")[0]
        n_count = self._unpack("I")[0]
        n_starts = list(self._unpack(f"{n_count}I")) if n_count else []
        n_sizes = list(self._unpack(f"{n_count}I")) if n_count else []
        mask_count = self._unpack("I")[0]
        if mask_count:
            self.handle.seek(8 * mask_count, os.SEEK_CUR)
        self.handle.seek(4, os.SEEK_CUR)
        rec = {
            "size": dna_size,
            "n_blocks": list(zip(n_starts, n_sizes)),
            "data_offset": self.handle.tell(),
        }
        self._records[chrom] = rec
        return rec

    def sequence(self, chrom, start, end):
        rec = self._record(chrom)
        if start < 0 or end > rec["size"] or end < start:
            raise RuntimeError("Invalid interval")
        first_byte = start // 4
        last_byte = (end + 3) // 4
        self.handle.seek(rec["data_offset"] + first_byte)
        packed = self.handle.read(last_byte - first_byte)
        chars = []
        for genome_pos in range(first_byte * 4, last_byte * 4):
            if start <= genome_pos < end:
                byte = packed[(genome_pos // 4) - first_byte]
                shift = (3 - (genome_pos % 4)) * 2
                chars.append(self.BASES[(byte >> shift) & 3])
        for n_start, n_size in rec["n_blocks"]:
            n_end = n_start + n_size
            if n_end <= start or n_start >= end:
                continue
            for pos in range(max(start, n_start), min(end, n_end)):
                chars[pos - start] = "N"
        return "".join(chars)

    def close(self):
        self.handle.close()


MOTIFS = {
    "AP1": ["TGACTCA", "TGAGTCA"],
    "CREB": ["TGACGTCA"],
    "ETS": ["GGAA", "GGAT"],
    "GATA": ["GATAA", "GATAG"],
    "EBOX": ["CACGTG", "CAGCTG"],
    "SP1": ["GGGCGG", "CCGCCC"],
    "NRF1": ["GCGCATGCGC"],
    "NFY": ["CCAAT"],
    "NFKB": ["GGGRNNYYCC".replace("R", "A").replace("N", "A").replace("Y", "C")],
    "CTCF": ["CCGCGNGGNGGCAG".replace("N", "A")],
    "FOX": ["TGTTTAC"],
    "SOX": ["AACAAT"],
    "IRF": ["GAAA"],
    "RUNX": ["TGTGGT"],
    "MEF2": ["CTAATTTTAG"],
    "P53": ["RRRCWWGYYY".replace("R", "A").replace("W", "A").replace("Y", "C")],
}


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def clean(seq):
    seq = seq.upper()
    if len(seq) != LENGTH or any(c not in "ACGT" for c in seq):
        return None
    return seq


def read_ccres():
    by_class = defaultdict(list)
    intervals = defaultdict(list)
    with open(CCRE_BED) as handle:
        for line in handle:
            chrom, start, end, *_rest, cls = line.rstrip("\n").split("\t")[:6]
            if chrom not in CANONICAL:
                continue
            start, end = int(start), int(end)
            if end <= start:
                continue
            by_class[cls].append((chrom, start, end))
            intervals[chrom].append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return by_class, intervals


def overlaps(intervals, chrom, start, end):
    rows = intervals.get(chrom, [])
    i = bisect.bisect_left(rows, (start, -1))
    if i and rows[i - 1][1] > start:
        return True
    return i < len(rows) and rows[i][0] < end


def fetch_centered(tb, rng, record, jitter=40):
    chrom, start, end = record
    mid = (start + end) // 2 + rng.randint(-jitter, jitter)
    s = max(0, mid - LENGTH // 2)
    e = s + LENGTH
    try:
        seq = tb.sequence(chrom, s, e)
    except RuntimeError:
        return None
    seq = clean(seq)
    if seq and rng.random() < 0.5:
        seq = revcomp(seq)
    return seq


def sample_ccres(tb, rng, by_class):
    seqs = []
    seen = set()
    for cls, n in CCRE_COUNTS.items():
        pool = by_class[cls][:]
        rng.shuffle(pool)
        made = 0
        attempts = 0
        while made < n and attempts < n * 30:
            rec = pool[attempts % len(pool)]
            seq = fetch_centered(tb, rng, rec)
            attempts += 1
            if seq and seq not in seen:
                seen.add(seq)
                seqs.append(seq)
                made += 1
        if made != n:
            raise RuntimeError(f"Only made {made}/{n} sequences for {cls}")
    return seqs, seen


def sample_background(tb, rng, intervals, seen, n):
    chrom_sizes = tb.chroms()
    sizes = {c: chrom_sizes[c] for c in CANONICAL if c in chrom_sizes}
    chroms = sorted(sizes)
    weights = [sizes[c] for c in chroms]
    seqs = []
    attempts = 0
    while len(seqs) < n and attempts < n * 200:
        chrom = rng.choices(chroms, weights=weights, k=1)[0]
        start = rng.randint(0, sizes[chrom] - LENGTH)
        attempts += 1
        if overlaps(intervals, chrom, start, start + LENGTH):
            continue
        seq = clean(tb.sequence(chrom, start, start + LENGTH))
        if not seq:
            continue
        if rng.random() < 0.5:
            seq = revcomp(seq)
        if seq not in seen:
            seen.add(seq)
            seqs.append(seq)
    if len(seqs) != n:
        raise RuntimeError(f"Only made {len(seqs)}/{n} background sequences")
    return seqs


def dinuc_shuffle(seq, rng):
    # Eulerian-style shuffle: preserve exact dinucleotide counts when possible.
    edges = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        edges[a].append(b)
    for vals in edges.values():
        rng.shuffle(vals)
    out = [seq[0]]
    for _ in range(len(seq) - 1):
        cur = out[-1]
        if edges[cur]:
            out.append(edges[cur].pop())
        else:
            out.append(rng.choice("ACGT"))
    return "".join(out)


def sample_shuffled(rng, source, seen, n):
    seqs = []
    attempts = 0
    while len(seqs) < n and attempts < n * 20:
        base = rng.choice(source)
        seq = dinuc_shuffle(base, rng)
        attempts += 1
        if seq not in seen:
            seen.add(seq)
            seqs.append(seq)
    if len(seqs) != n:
        raise RuntimeError(f"Only made {len(seqs)}/{n} shuffled sequences")
    return seqs


def random_background(rng, gc):
    p_gc = gc / 2
    p_at = (1 - gc) / 2
    bases = []
    for _ in range(LENGTH):
        x = rng.random()
        if x < p_at:
            bases.append("A")
        elif x < 2 * p_at:
            bases.append("T")
        elif x < 2 * p_at + p_gc:
            bases.append("C")
        else:
            bases.append("G")
    return bases


def mutate_motif(rng, motif, rate=0.12):
    chars = list(motif)
    for i, base in enumerate(chars):
        if rng.random() < rate:
            chars[i] = rng.choice([b for b in "ACGT" if b != base])
    return "".join(chars)


def sample_synthetic(rng, seen, n):
    programs = [
        ["AP1", "ETS", "GATA"],
        ["SP1", "NRF1", "EBOX"],
        ["CTCF", "CTCF"],
        ["NFY", "SP1", "ETS"],
        ["FOX", "SOX", "RUNX"],
        ["CREB", "AP1", "IRF"],
        ["P53", "EBOX", "SP1"],
        ["MEF2", "GATA", "AP1"],
    ]
    seqs = []
    attempts = 0
    while len(seqs) < n and attempts < n * 50:
        attempts += 1
        gc = min(0.78, max(0.22, rng.betavariate(2.2, 2.2)))
        chars = random_background(rng, gc)
        program = rng.choice(programs)
        copies = rng.randint(3, 9)
        occupied = []
        for _ in range(copies):
            name = rng.choice(program)
            motif = rng.choice(MOTIFS[name])
            motif = mutate_motif(rng, motif)
            if rng.random() < 0.5:
                motif = revcomp(motif)
            for _try in range(30):
                pos = rng.randint(5, LENGTH - len(motif) - 5)
                if all(abs(pos - old) > 6 for old in occupied):
                    chars[pos : pos + len(motif)] = list(motif)
                    occupied.append(pos)
                    break
        seq = "".join(chars)
        if rng.random() < 0.5:
            seq = revcomp(seq)
        seq = clean(seq)
        if seq and seq not in seen:
            seen.add(seq)
            seqs.append(seq)
    if len(seqs) != n:
        raise RuntimeError(f"Only made {len(seqs)}/{n} synthetic sequences")
    return seqs


def main():
    rng = random.Random(SEED)
    by_class, intervals = read_ccres()
    tb = TwoBit(HG38_2BIT)
    seqs, seen = sample_ccres(tb, rng, by_class)
    source_for_shuffle = seqs[:]
    seqs.extend(sample_background(tb, rng, intervals, seen, BACKGROUND_COUNT))
    seqs.extend(sample_shuffled(rng, source_for_shuffle, seen, SHUFFLE_COUNT))
    seqs.extend(sample_synthetic(rng, seen, SYNTH_COUNT))
    rng.shuffle(seqs)
    if len(seqs) != TARGET:
        raise RuntimeError(len(seqs))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as handle:
        for seq in seqs:
            handle.write(seq + "\n")
    tb.close()


if __name__ == "__main__":
    main()
