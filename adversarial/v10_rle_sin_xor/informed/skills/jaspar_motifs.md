# JASPAR motif loading and consensus sampling

## File
`data/JASPAR2024_CORE_vertebrates.meme` — 879 vertebrate TF binding profiles
downloaded from `https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt`.

## Format
MEME format. Each motif:
```
MOTIF MA0004.1 Arnt
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
 0.200000  0.800000  0.000000  0.000000
 ... (one row per position, p(A) p(C) p(G) p(T))
URL ...
```

## Parsing snippet
```python
def load_jaspar(path):
    motifs = {}  # id -> (name, np.array [w,4])
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('MOTIF '):
            parts = line.split()
            mid, name = parts[1], parts[2]
            # find probability matrix
            i += 1
            while i < len(lines) and not lines[i].startswith('letter-probability'):
                i += 1
            hdr = lines[i].split()
            w = int(hdr[hdr.index('w=') + 1])
            rows = []
            i += 1
            for _ in range(w):
                rows.append([float(x) for x in lines[i].split()])
                i += 1
            motifs[mid] = (name, np.array(rows))
        else:
            i += 1
    return motifs
```

## Sampling sequences from a PWM
```python
def sample_from_pwm(pwm, rng, with_revcomp=False):
    # pwm: [w, 4] each row sums to ~1
    bases = np.array(list("ACGT"))
    seq = "".join(bases[rng.choice(4, p=p) for p in pwm])
    if with_revcomp and rng.random() < 0.5:
        comp = str.maketrans("ACGT", "TGCA")
        seq = seq.translate(comp)[::-1]
    return seq
```

## Notable HepG2-relevant TFs (hepatocyte)
HNF4A, HNF1A, FOXA1, FOXA2, CEBPA, CEBPB, ATF4, NR1I2, NR1I3

## Notable K562-relevant TFs (erythroid/myeloid)
GATA1, GATA2, KLF1, TAL1, RUNX1, SPI1, MYB, NFE2

## Notable SK-N-SH-relevant TFs (neuroblastoma)
ASCL1, NEUROD1, REST (silencer), MEF2C, MYCN, PHOX2B, TFAP2A
