# MPRA Library Design Notebook

## Initial framing

The objective is not to maximize activity in K562, HepG2, or SK-N-SH directly.
The library should train a model that learns transferable regulatory sequence
grammar. I interpret that as needing many examples from real regulatory DNA,
spread across many regulatory programs, while retaining enough negative and
edge-case sequence diversity that the model does not only learn "open
chromatin-like" composition.

The prior runs strongly favor DHS-derived sequence at 50k. The best 50k method
was DHS sampled by topic loading, and the next-best methods all retained a
large DHS component. Fully random sequence helped eval_08 but usually hurt the
other test sets, so synthetic sequence should be used as a minority component
rather than a large replacement for genomic regulatory DNA.

## Candidate sources considered

- Meuleman et al. DHS index/vocabulary: include as the primary source. It is
  broad across many biosamples, has a 16-topic NMF representation, and matches
  the strongest prior strategy conceptually.
- ENCODE/S creen cCREs or other chromatin-state catalogs: useful for enhancer,
  promoter, and CTCF class diversity, but the prior results show that replacing
  too much DHS with SEI-like state catalogs can reduce 50k performance. I will
  use DHS itself first, and only add extra non-DHS material if easy to obtain
  without increasing risk.
- Existing MPRA sequence sets: exclude as a main source. The prior MPRA oracle
  result was substantially worse at 50k, likely because its distribution is
  narrower than genome-wide regulatory grammar.
- Pure random synthetic sequence: include only as a small component. It covers
  sequence space and helps one eval set, but large fractions dilute natural
  motif syntax.

## Current design hypothesis

Use a mostly DHS-topic library with four parts:

1. Topic-weighted DHS centers, reproducing the strongest baseline signal.
2. Topic-stratified DHS centers, ensuring rare regulatory programs are not
   lost when sampling only by maximum loading or total loading.
3. Boundary/offset DHS windows, so the model sees motif grammar in nearby
   weak flanks instead of only summit-centered high-confidence windows.
4. A small synthetic tail: dinucleotide/genomic-background random sequences
   and motif-spiked sequences with varied motif combinations, intended to
   add controlled coverage of composition and motif interaction space.

The library size is fixed at 50,000. My initial allocation target is about
80-90% DHS-derived sequence and 10-20% synthetic/perturbed sequence. I will
avoid exact duplicate sequences and filter non-ACGT bases. The final generator
will be deterministic.

## Data acquisition started

I started downloading:

- `DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz`
- `2018-06-08NC16_NNDSVD_Mixture.npy.gz`
- `hg38.2bit`

The plan is to sample DHS rows using the NMF mixture matrix, extract 200 bp
windows from hg38, and then fill any remaining quota with deterministic
synthetic designs.

## Generator implementation

`py2bit` could not be installed in this environment because the Python
development headers are unavailable, so I implemented a small read-only 2bit
parser directly in `library/generate.py`. This keeps the generator
self-contained and avoids converting the whole genome to FASTA.

The downloaded DHS table has the expected columns:

- chromosome, start, end, identifier
- mean DNase signal
- number of biosamples
- summit and core coordinates
- vocabulary component

The NMF mixture file is a 16 x 3,591,898 matrix. The generator transposes this
to DHS rows x 16 topics when the row count matches the DHS table.

Final allocation implemented:

- Up to 37,500 summit-centered DHS windows sampled by a weight combining
  mean signal, biosample breadth, maximum topic strength, and total topic
  loading.
- Additional summit-centered DHS windows sampled separately from each of the
  16 topics, filling to 42,500. This is meant to keep the strongest
  DHS-topic behavior while making rare topics less likely to disappear.
- Offset windows from high-confidence/broad DHSs at +/-75 and +/-125 bp,
  filling to 47,000. These are still biological DHS sequence, but not exact
  summit replicas.
- Component-balanced DHS backfill to 48,000.
- 2,000 synthetic/perturbed sequences: dinucleotide-like shuffles of selected
  DHS windows first, then GC-varied motif-spiked random sequences. The motif
  list includes generic AP-1, E-box, GC-box, ETS, CCAAT/NF-Y, TATA, GATA,
  FOX/SOX, nuclear receptor, NRF1-like, and CTCF-core patterns. This is a
  deliberately small tail because the prior 50% synthetic baseline helped
  eval_08 but slightly reduced overall mean performance.

Validation after generation:

- 50,000 lines
- all lines length 200
- all bases in {A,C,G,T}
- 50,000 unique sequences
- mean GC fraction about 0.483

I will now make the single allowed `prepare.py` call.

## Final evaluation

I ran the one allowed command:

`python3 prepare.py library/sequences.txt`

The run took about 1293 seconds and wrote `library/result.json`.

Final mean_r by anonymous eval set:

| eval | mean_r |
|---|---:|
| eval_01 | 0.7505 |
| eval_02 | 0.8459 |
| eval_03 | 0.8297 |
| eval_04 | 0.8076 |
| eval_05 | 0.7502 |
| eval_06 | 0.8458 |
| eval_07 | 0.7894 |
| eval_08 | 0.7304 |
| eval_09 | 0.8799 |
| eval_10 | 0.8234 |
| eval_11 | 0.7375 |
| eval_12 | 0.7121 |
| eval_13 | 0.7824 |
| eval_14 | 0.8463 |

Average of the 14 mean_r values: 0.7951.

Compared with the prior 50k table, this design improved on the listed best
baseline for every eval set. The largest gains appear on eval_08 and the
harder lower-correlation evals, which supports the decision to add offset
DHS windows plus a small perturbation/synthetic tail rather than exactly
recreating summit-only DHS-topic sampling.

## What I would try next

If another shot were available, I would treat this result as evidence that
the main opportunity is still within DHS sequence selection, not wholesale
replacement with other catalogs. I would next test:

- A controlled allocation sweep around the synthetic/perturbation tail
  between 2% and 10%.
- More principled DHS offset sampling based on DHS core width and summit
  uncertainty instead of fixed +/-75 and +/-125 bp offsets.
- Component-aware constraints during the first sampling stage rather than
  only as backfill.
- Additional natural non-DHS negatives from matched genomic flanks, avoiding
  repetitive or ambiguous sequence.
