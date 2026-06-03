# Skill: Parse JASPAR PFMs and sample motif instances

**Input**: `data/JASPAR_vert.txt` — JASPAR 2024 CORE vertebrate non-redundant PFMs (879 motifs).

**Format**: each motif is 5 lines —
```
>MA0004.1	Arnt
A  [     4     19      0      0      0      0 ]
C  [    16      0     20      0      0      0 ]
G  [     0      1      0     20      0     20 ]
T  [     0      0      0      0     20      0 ]
```
Columns are positions; rows are counts. Some have many more positions.

**Parser** (use this exact code):
```python
import re
import numpy as np
from pathlib import Path

def parse_jaspar(path):
    """Return list of (id, name, pwm) where pwm is shape (4, L) probabilities."""
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
                    counts = np.stack([rows[b] for b in 'ACGT']).astype(float)
                    counts += 0.01  # pseudocount
                    pwm = counts / counts.sum(axis=0, keepdims=True)
                    motifs.append((header[0], header[1], pwm))
                parts = line[1:].split('\t')
                header = (parts[0], parts[1] if len(parts) > 1 else parts[0])
                rows = {}
            else:
                m = re.match(r'([ACGT])\s*\[(.*)\]', line)
                if not m:
                    continue
                base = m.group(1)
                vals = [float(x) for x in m.group(2).split()]
                rows[base] = np.array(vals)
        if header and len(rows) == 4:
            counts = np.stack([rows[b] for b in 'ACGT']).astype(float)
            counts += 0.01
            pwm = counts / counts.sum(axis=0, keepdims=True)
            motifs.append((header[0], header[1], pwm))
    return motifs
```

**Sample a motif instance**: given a PWM of shape (4, L), draw one base per column:
```python
def sample_from_pwm(pwm, rng):
    L = pwm.shape[1]
    out = np.empty(L, dtype='U1')
    for j in range(L):
        out[j] = "ACGT"[rng.choice(4, p=pwm[:, j])]
    return "".join(out)
```

**Vectorised version** (faster for many draws):
```python
def sample_many(pwm, n, rng):
    """Return list of n sampled motif instances."""
    L = pwm.shape[1]
    # Per-column draw using inverse CDF
    cdf = pwm.cumsum(axis=0)  # (4, L)
    u = rng.random((n, L))
    # broadcast: for each (i, j), pick smallest k where cdf[k, j] >= u[i, j]
    bases_idx = (u[:, None, :] >= cdf[None, :, :]).sum(axis=1)  # 0..3
    chars = np.array(list("ACGT"))
    return ["".join(chars[bases_idx[i]]) for i in range(n)]
```

**Tips**:
- Pseudocount 0.01 avoids log(0) and zero-probability draws.
- JASPAR 2024 vertebrate non-redundant has 879 motifs, length 4–30bp (most 6–15).
- PWM lengths vary; check `pwm.shape[1]` before placement so motif fits within 200bp window.
- TF coverage: covers ~600 unique TFs across vertebrate lineages — broader than what's expressed in any one cell line. Good for generalization-by-grammar libraries.
