# Experiment 003 — Gosai MPRA, train chromosomes

## Hypothesis
v08's eval distribution comes from the Gosai 2024 MPRA dataset (K562/HepG2/SK-N-SH,
200bp — the only public dataset matching all v08 cell types). 50K random Gosai
sequences from chromosomes 1-6, 8, 10-12, 14-20, 22 (excluding likely held-out
chroms 7, 9, 13, 21, X, Y) should yield a substantially higher eval_01 than random
DHS/cCRE.

## Method
- Downloaded Table_S2__MPRA_dataset.txt (798K Gosai sequences)
- Filtered to 200bp ACGT sequences from non-held-out chromosomes (647K)
- Sampled 50K uniformly at random

## Result
**eval_01 = 0.0000.** Most eval sets still near zero, BUT:
- eval_04 = 0.0131 (up from -0.0018 random baseline)
- eval_07 = 0.0091
- eval_09 = 0.0131 (eval_04 ≡ eval_09 pattern)
- eval_08 = 0.0030

So a real but tiny signal emerged for some eval sets. The Gosai distribution
partially aligns with eval — but only on certain eval indices.

## Interpretation
Direction is right but magnitude is small. Possible reasons:
1. The eval composition is more specific (e.g., only CRE-class, or specific
   activity range)
2. Random Gosai sampling dilutes signal with low-info variant sequences
3. Need quality filtering (low measurement SE) and activity diversity

## Theory update
The eval likely uses Gosai-like sequences but with composition different from
random. To maximize signal, we should preferentially sample:
- Sequences with low measurement error (lfcSE < 1.0, matching Malinois training)
- Sequences with diverse / extreme activity (informative training data)
- Possibly weight toward the CRE-class subset (14K curated)
