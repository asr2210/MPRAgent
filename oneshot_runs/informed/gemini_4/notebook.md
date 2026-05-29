# Lab Notebook — MPRA Library Design for Gene Regulatory Grammar

## 1. Theory of What Makes a Good MPRA Training Library and Why

Massively Parallel Reporter Assays (MPRAs) measure how DNA sequences drive gene expression by looking at the transcriptional output of hundreds of thousands of candidate regulatory elements. To train a high-performing deep learning model that decodes gene regulatory grammar across *all* human tissues and cell lines, the library's sequence composition must satisfy several orthogonal objectives:

1.  **Capturing Specificity and Active Grammar (Active Space):** The library must contain active regulatory sequences (promoters, enhancers, enhancers with cell-type-specific activity) that represent the functional TF binding sites, spacing, and combinations found in nature. This teaches the model the "positive" regulatory grammar.
2.  **Covering Broad Genomic Context (General Space):** The model must generalize beyond specific cell-type-specific programs. It needs a diverse set of genomic accessible regions representing housekeeping or universally open elements across different chromatin states.
3.  **Establishing a Robust Baseline & Negative Space (Background Noise):** Deep learning models are prone to learning spurious features (e.g., GC-content biases, assembly artifacts) if trained only on active genomic regions. A good training library must include inactive or fully random sequences. This negative control teaches the model what *does not* drive activity, preventing false positives and grounding the model's predictions.

## 2. Analysis of Baselines and Trade-offs

We began by parsing and analyzing the 14 systematic baseline runs provided in the instructions (Table 1):
*   **`dhs_topic` (Mean r = 0.7630):** The top overall baseline, performing exceptionally well on 9 out of 14 evaluation sets. It samples DHS regions with probability proportional to 16 NMF topic loadings, heavily upweighting elements with strong cell-type-specific activity. This indicates that tissue-specific chromatin signatures are highly informative for general regulatory grammar.
*   **`dhs_synth` (Mean r = 0.7591):** Consists of 50% `dhs_topic` and 50% fully random synthetic sequences. It performs nearly identically to `dhs_topic` on almost all biological sets, but scores a massive **0.7523** on `eval_08` (compared to `dhs_topic`'s **0.7011**). This demonstrates that introducing synthetic sequences provides an immense boost to synthetic-sensitive test sets (or general background learning) with an incredibly small penalty on biological sets.
*   **`dhs_random` (Mean r = 0.7503):** Performs best on `eval_13` (**0.7639** vs `dhs_topic`'s **0.7271**). This suggests that some test sets are composed of highly diverse or non-specific genomic elements, where the diversity of the unweighted DHS pool is much more useful than narrow cell-type-specific elements.

These observations led to our key hypothesis: **A hybrid library combining topic-specific active elements, diverse genomic accessible background, and fully random synthetic sequences will outperform any single-component strategy.**

## 3. Exploratory Data Analysis & Empirical Discoveries

We downloaded and analyzed the official Meuleman et al. 2020 DNase I Hypersensitive Sites (DHS) Index and the NMF topic loadings (Mixture matrix). 

### Analysis 1: Component Distribution and Signal Strengths
Running our `analyze_components.py` script on the 3.59 million DHS elements revealed the exact distributions across the 16 NMF components:
*   **Constitutive Components:** `Stromal A` (Mean Signal: 2.08) and `Tissue invariant` (Mean Signal: 1.49) have extremely high signals and are active in over 117+ biosamples, acting as housekeeping/universal regulatory regions.
*   **Specific Components:** Components like `Neural` and `Primitive / embryonic` have lower signals and are active in fewer biosamples, representing cell-type-specific active chromatin programs.

### Analysis 2: Correlation of NMF Loadings and Signals
We analyzed the correlation between continuous NMF loadings, `numsamples` (number of active biosamples), and `mean_signal` (signal magnitude) for the first 100,000 DHS elements:
*   **Sum of NMF Loadings vs. `numsamples`:** Pearson $r = \mathbf{0.9766}$ (almost perfect correlation). The sum of NMF topic loadings is a direct proxy for how many tissues the DHS is open in (tissue ubiquity).
*   **Sum of NMF Loadings vs. `mean_signal`:** Pearson $r = 0.2804$.
*   **Max of NMF Loadings vs. `mean_signal`:** Pearson $r = 0.3131$.

**Key Insight:** The `mean_signal` column provides highly independent, orthogonal information representing the *confidence and magnitude* of accessibility. A tissue-specific element (high $L_{\text{topic}}$) with a high `mean_signal` is an extremely strong, functional enhancer. A tissue-specific element with a weak `mean_signal` might represent weak or noisy accessibility.

## 4. Specific Design Decisions & Implementation

Based on our EDA, we developed a **Multi-objective Stratified Diversity-Maximized Hybrid Library** of exactly 50,000 sequences of 200bp length:

### A. Library Composition (70% / 15% / 15%)
1.  **70% (35,000 sequences) Topic-Specific High-Signal DHS ("Cream of the Crop"):**
    *   We divide the 35,000 sequences equally across the 16 NMF components (~2,187 or 2,188 per component).
    *   For each component $k$, we define the sampling weight for each DHS $j$ as:
        $$W_{k, j} = L_{k, j} \times \text{mean\_signal}_j$$
        This multiplication ensures we select elements that are both highly specific to the given biological program and have strong, high-confidence regulatory activity.
2.  **15% (7,500 sequences) Broad Genomic DHS (Chromatin Diversity):**
    *   We sample 7,500 sequences uniformly from the entire DHS index pool. This covers the general genomic open chromatin space, preserving performance on diverse test sets like `eval_13`.
3.  **15% (7,500 sequences) Synthetic Sequences (Sequence Space Diversity):**
    *   We generate 7,500 fully random sequences (i.i.d. uniform drawing {A, C, G, T}). This provides the ultimate negative/inactive background, helping with sets like `eval_08` and preventing false positives.

### B. High-Fidelity Sequence Extraction
*   **Summit-Centering:** Instead of taking arbitrary coordinates, we use the `summit` column of the DHS index to extract exactly 200bp centered at the peak of accessibility ($[\text{summit} - 100, \text{summit} + 100]$). This places the core transcription factor binding sites and regulatory motifs directly in the middle of our 200bp window, maximizing motif density and training signal.
*   **Assembly Gap ('N') Filtering:** Any extracted sequence containing any assembly gap characters (`N` or `n`) is immediately discarded and replaced with another sampled element.
*   **Strict Uniqueness:** We enforce 100% uniqueness of sequences across the entire library. Any duplicate sequence is rejected, maximizing sequence-space representation.
*   **Deterministic Reproducibility:** We set fixed random seeds (`numpy` and `random` to `42`) to guarantee identical library generation.

## 5. What to Try Next if We Had Another Shot

If we had another iteration, we would explore the following high-impact directions:
1.  **In-Silico Model Screening (Active Learning):** Use a pre-trained genomic model (like Enformer, Sei, or DNABERT) to predict the regulatory activity of millions of sequences, then specifically curate a training set with high "information content" (e.g., highly predicted active sequences, high-entropy sequences, or those with high predicted cell-type variance).
2.  **Explicit Motif density Balancing:** Scan DHS sequences for known TF motifs using Position Weight Matrices (PWMs) and balance the library to ensure every human TF family (e.g., bHLH, Zinc Fingers, Homeobox) has adequate motif representation.
3.  **Controlled Synthetic Perturbations:** Instead of fully random synthetic sequences, design synthetic sequences containing specific known regulatory motifs embedded in random backgrounds at varying positions and spacings. This would directly probe the positional and combinatorial logic of regulatory grammar.
