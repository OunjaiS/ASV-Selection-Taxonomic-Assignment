# ASV Authentication Framework for Single-Specimen COI Barcoding

Analysis scripts for the manuscript:

> **Diversity and Phylogenetic Placement of Amplified Sequence Variants in Single-Specimen Barcoding of Tropical Insects**
> Sarawut Ounjai, Huaxi Liu, Zichen Zhou, Maria Pestana Correia, Thomas J. Creedy, Carmelo Andújar, Paula Arribas, Alfried P. Vogler
> *Running title: ASV Authentication in Tropical Beetles*

---

## Overview

This repository contains the analysis scripts used to authenticate COI Amplicon Sequence Variants (ASVs) from single-specimen barcoding of tropical Coleoptera. The framework integrates abundance filtering, phylogenetic placement, codon usage bias analysis, and machine learning to classify 175,954 ASV records from 18,533 specimens across 14 countries.

**Key results:**
- 86.52% authentication success (15,901 authenticated ASV records from 14,715 specimens)
- 4 failure categories: Technical Artefacts (58.02%), Environmental Contamination (14.20%), Intra-individual Variants/NUMTs (11.27%), Cross-sample Contamination (7.48%)
- NUMT validation via codon usage bias: 1,790 NUMTs validated across 4,496 specimens
- DADA2 benchmarking: 75,481 NUMT/Contamination records passed DADA2 filters but failed phylogenetic authentication

> **Raw data and sequences** are archived on Dryad: https://doi.org/10.5061/dryad.v41ns1s9q

---

## Repository Structure

```
script/
├── 1_sequences_analysis.ipynb          # Sequence QC and composition
├── 2_classification_analysis.ipynb     # Classification pipeline
├── 3_statistical_analysis.ipynb        # Statistical tests
├── DADA2_Benchmarking_Comparison.py    # DADA2 vs framework comparison
├── NUMT_Codon_Usage_Analysis.py        # NUMT codon usage validation
└── 5_dada2_chimera_detection.R         # DADA2 chimera detection (R)
```

---

## Scripts — Detailed Description

### 1. `1_sequences_analysis.ipynb`
**Manuscript section:** Methods 2.3, Results 3.1–3.2

Performs sequence quality control and composition analysis on the raw ASV sequences.

**What it does:**
- Translates each COI sequence in all 6 reading frames and validates the correct open reading frame
- Calculates nucleotide composition (GC content, AT skew, codon position GC%)
- Calculates amino acid composition (hydrophobicity, leucine content, molecular weight)
- Performs motif analysis for conserved COI domains (TGA/TGG codon validation, transmembrane helices)
- Assigns sequence quality grade (A+, A, B, C, Fail) based on multiple criteria
- Outputs `ASV_Complete_Analysis.csv` and `Analysis_Summary.csv`

**Input files (download from Dryad — https://doi.org/10.5061/dryad.v41ns1s9q):**
- `ASV_sequences_64544.fasta` — unique ASV sequences (FASTA format, 64,544 sequences)

---

### 2. `2_classification_analysis.ipynb`
**Manuscript section:** Methods 2.3–2.5, Results 3.3–3.4, Table 1, Table 2, Figure 4

Implements the multi-stage ASV authentication and classification pipeline.

**What it does:**
- Applies abundance thresholding (MRCT = 4 reads) to remove low-count ASVs
- Classifies each ASV into one of five categories: Authenticated, Technical Artefacts, Environmental Contamination, Intra-individual Variants (NUMTs), Cross-sample Contamination
- Implements rule-based classification using phylogenetic distance + read abundance criteria
- Trains and evaluates a Random Forest classifier for ASV-level feature importance
- Performs co-occurrence analysis to identify stable NUMTs present in ≥2 specimens
- Outputs `ASV_Final_Classification.csv`, `Classification_Statistics.csv`, `Feature_Importance.csv`

**Input files (download from Dryad):**
- `ASV_Authentication_Results_030925.csv` — Main dataset: 175,954 ASV records with phylogenetic distance, read counts, specimen ID, family assignment, and collection metadata

---

### 3. `3_statistical_analysis.ipynb`
**Manuscript section:** Results 3.5, Table S1–S5

Performs all statistical comparisons reported in the manuscript.

**What it does:**
- Tests differences in authentication success by country, collection method, sequencing batch, and taxonomic family
- Uses chi-squared tests with Cramér's V effect sizes for categorical comparisons
- Uses Spearman correlation for continuous predictors (sequencing depth, batch size)
- Generates Table S1 (country success), Table S2 (method success), Table S3 (batch success), Table S4 (family success), Table S5 (effect sizes for all comparisons)

**Input files (download from Dryad):**
- `ASV_Authentication_Results_030925.csv` — Main dataset with collection metadata (country, method, batch, family)

---

### 4. `NUMT_Codon_Usage_Analysis.py`
**Manuscript section:** Methods 2.5, Results 3.6, Table 3, Figure 7

Performs codon usage bias analysis to provide molecular evidence supporting NUMT classification.

**What it does:**
- Compares Relative Synonymous Codon Usage (RSCU) between Authenticated ASVs (n = 8,533 unique sequences) and putative NUMTs (n = 26,416 unique sequences)
- Calculates Effective Number of Codons (ENC) as a measure of codon bias
- Performs Mann-Whitney U tests for nucleotide composition and amino acid composition (Table 3A & B)
- Applies Benjamini-Hochberg correction across 62 sense codons (Supplementary Table S7)
- Performs PCA on per-ASV RSCU values to visualise group separation (Figure 7C)
- Generates Figure 7 (4-panel: GC violin, ENC violin, PCA, RSCU bar chart)

**Input files (download from Dryad):**
- `ASV_Authentication_Results_030925.csv` — Full classification with ASV sequences
- `ASV_sequences_64544.fasta` — ASV sequences for codon table generation

---

### 5. `DADA2_Benchmarking_Comparison.py`
**Manuscript section:** Methods 2.3, Results 3.4, Discussion 4.5, Supplementary Figure S2, Supplementary Table S6

Benchmarks the DADA2 denoising pipeline against the phylogenetic authentication framework.

**What it does:**
- Compares feature-by-feature: DADA2 (5/15 features) vs this framework (15/15 features)
- Runs MRCT sensitivity analysis showing the tradeoff between sensitivity and artefact removal
- Quantifies the "authentication gap": 75,481 records (NUMTs + Environmental Contamination) that pass DADA2 denoising but fail phylogenetic authentication
- Generates Supplementary Figure S2 (4-panel benchmarking figure)
- Generates Supplementary Table S6 (DADA2 analysis summary)

**Input files (download from Dryad):**
- `ASV_Authentication_Results_030925.csv` — Main classification results
- Requires `5_dada2_chimera_detection.R` output (`dada2_chimeras.csv`, `dada2_nonchimeras.csv`)

---

### 6. `5_dada2_chimera_detection.R`
**Manuscript section:** Methods 2.3, Results 3.4, Discussion 4.5

Runs DADA2 chimera detection (`isBimeraDenovo`) on the ASV sequences using the R `dada2` package.

**What it does:**
- Reads unique ASV sequences in FASTA format
- Runs DADA2's bimera (two-parent chimera) detection algorithm
- Outputs two CSV files: `dada2_chimeras.csv` and `dada2_nonchimeras.csv`
- These are used as input for `DADA2_Benchmarking_Comparison.py`

**Note:** This is the improved version (v2) that uses unique sequence deduplication. The original v1 used abundance-weighted mapping.

**Input files (download from Dryad):**
- `ASV_sequences_64544.fasta` — Unique ASV sequences

**R dependencies:**
```R
install.packages("BiocManager")
BiocManager::install("dada2")
```

---

## Installation

```bash
git clone https://github.com/OunjaiS/ASV-Selection-Taxonomic-Assignment
cd ASV-Selection-Taxonomic-Assignment

# Python dependencies
pip install pandas numpy scipy scikit-learn matplotlib seaborn biopython jupyter

# R dependencies
# In R: BiocManager::install("dada2")
```

---

## Running Order

```bash
jupyter notebook

# Suggested order:
# 1. script/1_sequences_analysis.ipynb        → generates sequence quality data
# 2. script/2_classification_analysis.ipynb   → generates classification results
# 3. script/3_statistical_analysis.ipynb      → generates Table S1–S5
# 4. script/5_dada2_chimera_detection.R       → generates DADA2 chimera results (requires R)
# 5. python script/DADA2_Benchmarking_Comparison.py
# 6. python script/NUMT_Codon_Usage_Analysis.py
```

---

## Data Availability

| Resource | Location |
|---|---|
| Raw sequences & authentication results | [Dryad](https://doi.org/10.5061/dryad.v41ns1s9q) |
| Analysis scripts | This repository |
| Code archive (Zenodo) | https://doi.org/10.5281/zenodo.XXXXXXX |
| Specimen images | [Flickr](https://www.flickr.com/photos/site-100/) |
| Taxonomic data | [iNaturalist](https://www.inaturalist.org/people/site_100) |

---

## Citation

> Ounjai S, Liu H, Zhou Z, Pestana Correia M, Creedy TJ, Andújar C, Arribas P, Vogler AP.
> *Diversity and Phylogenetic Placement of Amplified Sequence Variants in Single-Specimen Barcoding of Tropical Insects.*
> (in review, 2026)

---

**Keywords:** DNA barcoding, COI, ASV authentication, NUMTs, DADA2, phylogenetic placement, tropical Coleoptera, codon usage bias
