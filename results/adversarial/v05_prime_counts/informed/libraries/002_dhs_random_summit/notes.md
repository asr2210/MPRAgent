# 002_dhs_random_summit

## Design
50k 200bp windows from hg38, each centered on a randomly-sampled DHS summit
(Meuleman 2020 index, 3.59M canonical-chrom elements). Discard windows
containing N or that fall off chromosome ends.

## Hypothesis
Real genomic context should beat synthetic. instructions.md baseline claims
dhs_random eval_01 = 0.7089. So I expect 0.6-0.7 here.

## Result
**eval_01 = 0.0410** — same magnitude as exp 001 (motif cocktail) and
strategies.md `random_uniform` (0.0399). DHS provides NO advantage over
random in this environment.

## Interpretation
The instructions.md baseline table is misleading for this prepare.py.
strategies.md is the reliable reference: all standard libraries hit a
~0.04 ceiling. The eval must have unusual structure that DHS doesn't
capture.

Notable: eval_08 = 0.045 here. random_uniform got 0.121 on eval_08.
So DHS is *worse* than random on eval_08 — likely because DHS sequences
have learned biases (motif content, GC isochores, repeat elements) that
the random-trained model lacks but the eval_08 eval set distribution
shares with random sequences.

## What to try next
Probe with promoters (Exp 003), ENCODE cCREs, and high-redundancy
libraries to triangulate what eval distribution looks like.
