# Lab Notebook

## Starting objective and assumptions
The task is to design a 50,000-sequence, 200 bp MPRA library for training a model of regulatory activity. The evaluation is a black box and must be treated like a wet-lab assay; I will not inspect `prepare.py`.

My working theory is that the best one-shot library should combine: (1) real human regulatory DNA sampled broadly across regulatory classes, because natural sequence context captures motif syntax, spacing, CpG composition, repeats, local grammar, and weak/background motifs; (2) genomic/random-like negatives and weak elements, because a model trained only on functional elements will learn an unrealistic prior; and (3) synthetic motif perturbation sequences, because a compact library should deliberately cover TF motif identity, copy number, strand, spacing, GC background, motif combinations, and promoter/enhancer-like architectures more evenly than random genomic sampling.

The target is not K562/HepG2/SK-N-SH specificity. I will therefore avoid over-weighting those lines and instead use pan-cell-type cCRE classes plus general vertebrate TF motifs. Evaluation includes three measured cell types, but a robust model should benefit from examples spanning promoters, enhancers, CTCF/insulators, accessible-but-unannotated regions, neutral genomic-like backgrounds, and synthetic grammar probes.

## Initial data plan
I will use ENCODE SCREEN / Registry cCRE class BED files for GRCh38 as broad regulatory coordinate sources: promoter-like, enhancer-like, CTCF-associated, accessible-only, accessible-with-TF, TF-only, and similar classes. These are not specific to one tissue and cover many cell and tissue types. If I can extract sequence context from hg38 reliably, I will sample 200 bp windows centered or jittered around those elements. If genome download becomes a bottleneck, I will fall back to a synthetic-only library with realistic motif/background variation.

I will also generate synthetic sequences from a curated motif vocabulary representing common regulatory families: promoter core and CpG features, CTCF, ETS, AP-1, GATA, FOX, nuclear receptors, bHLH, RUNX, NF-kB, STAT, SMAD, TEAD, SOX, OCT, CREB/ATF, CEBP, KLF/SP, YY1, REST/NRSE, p53, MYC/MAX, MEF2, RFX, IRF, HNF, and others. The synthetic design should include motif dropout/mutation controls and multiple background GC bins.

## Concrete library architecture
I downloaded the ENCODE cCRE V4 GRCh38 BED file from the Moore/Weng lab supplementary data mirror and am downloading the UCSC hg38 2bit genome. The BED file has approximately 2.35 million cCREs with class counts dominated by distal enhancer-like signatures (dELS), but also promoter-like signatures (PLS), proximal enhancer-like signatures (pELS), CTCF-associated accessible elements, accessible-only elements, and TF-only classes. This breadth is useful because the goal is not a model for one tissue or one active set; it should see many sequence regimes.

I will allocate most of the library to 200 bp genomic windows sampled from cCREs with jitter around the element center. I am intentionally not using only the highest-confidence active classes. Accessible-only, TF-only, and CTCF classes are valuable because they teach the model insulator/boundary grammar and partially active or context-dependent regulatory DNA. I will add random genomic windows as background controls; they are important for calibrating the model against ordinary genomic sequence, including realistic GC distribution and low motif density.

The synthetic component will cover motif grammars that random genomic sampling may under-cover in only 50,000 probes: motif copy number, co-occurrence, spacing, strand, GC background, promoter-like core architectures, CTCF/insulator-like motifs, and mutated controls. The motifs are consensus/IUPAC abstractions, not precise PWMs, because the objective is broad coverage rather than optimizing for one motif database. This should still expose the model to key sequence words and combinations.

Planned final allocation:
- 36,000 cCRE-derived genomic sequences, stratified across cCRE classes.
- 6,000 random genomic background windows.
- 8,000 synthetic motif-grammar sequences.

I expect this mixed design to have a better training-performance-to-size ratio than any single source alone: natural cCREs provide real syntax and local sequence context, random genomic windows provide negatives and calibration, and synthetic designs deliberately fill the combinatorial motif space.

## Generation and validation
I implemented `library/generate.py` as a deterministic generator with seed 20260522. It reads `data/GRCh38-cCREs-V4.bed.gz` and `data/hg38.2bit` when present, samples 36,000 cCRE-centered/jittered genomic windows stratified by class, samples 6,000 random genomic windows, and adds 8,000 synthetic motif grammar probes. The script deduplicates and fills any shortfall with additional synthetic enhancer-like probes.

The first generated `library/sequences.txt` passed the hard validation checks:
- 50,000 lines.
- 50,000 unique sequences.
- Every sequence is exactly 200 bp.
- Alphabet is limited to A/C/G/T.
- Overall GC fraction is approximately 0.473.

The GC value is a reasonable compromise: it is not too promoter/CpG-heavy, and it still includes higher-GC synthetic promoter-like probes and cCRE classes. This should help prevent the model from treating high GC as the only regulatory signal while still giving it enough CpG/promoter grammar.

## Evaluation result
I ran `python3 prepare.py library/sequences.txt` once, after the hard validation checks passed. The evaluator completed in 544 seconds and wrote `library/result.json`.

Per-test mean correlations:
- eval_01: 0.6981
- eval_02: 0.7859
- eval_03: 0.7639
- eval_04: 0.7518
- eval_05: 0.6984
- eval_06: 0.7860
- eval_07: 0.7188
- eval_08: 0.6321
- eval_09: 0.8129
- eval_10: 0.7323
- eval_11: 0.6857
- eval_12: 0.6569
- eval_13: 0.7103
- eval_14: 0.7858

The unweighted average of the 14 `mean_r` values is approximately 0.7299. The weakest sets were eval_08, eval_12, and eval_11. Since I cannot iterate after the committed result, I will keep this design. If I had another shot, I would try to diagnose whether those lower-scoring evals represent a missing class such as repetitive/low-complexity sequence, specific promoter architecture, long-range motif spacing not captured in 200 bp, or a domain where purely natural cCREs are underperforming relative to more synthetic perturbation coverage. I would also consider a small family of alternate allocations, especially increasing random genomic and synthetic negative controls at the expense of dELS-heavy cCRE sampling.
