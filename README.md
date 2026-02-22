# ASV Authentication Framework for Single-Specimen COI Barcoding

Analysis scripts for the manuscript:

> **Diversity and Phylogenetic Placement of Amplified Sequence Variants in Single-Specimen Barcoding of Tropical Insects**
> Sarawut Ounjai, Huaxi Liu, Zichen Zhou, Maria Pestana Correia, Thomas J. Creedy, Carmelo Andújar, Paula Arribas, Alfried P. Vogler
> *Running title: ASV Authentication in Tropical Beetles*

---

## Overview

This repository contains the analysis scripts used to authenticate COI Amplicon Sequence Variants (ASVs) from single-specimen barcoding of tropical Coleoptera. The framework integrates abundance filtering, phylogenetic placement, codon usage bias analysis, and machine learning to classify 175,954 ASV records from 18,533 specimens across 14 countries.

**Key results:**
- Authentication success: 86.52% of quality-passing specimens (15,901 authenticated ASV records)
- 4 failure categories identified: Technical Artefacts (58.02%), Environmental Contamination (14.20%), Intra-individual Variants / NUMTs (11.27%), Cross-sample Contamination (7.48%)
- NUMT validation via codon usage bias: 1,790 NUMTs validated across 4,496 specimens
- DADA2 benchmarking: 75,481 records passed DADA2 filters but failed phylogenetic authentication

> **Raw data, sequences, and supplementary tables** are archived on Dryad:
> https://doi.org/10.5061/dryad.v41ns1s9q

---

## Repository Structure

```
script/
├── 1_sequences_analysis.ipynb        # ASV sequence QC, nucleotide/amino acid composition
├── 2_classification_analysis.ipynb   # Classification pipeline + ML feature importance
├── 3_statistical_analysis.ipynb      # Statistical tests (Mann-Whitney U, chi-squared, effect sizes)
└── 3_statistical_analysis.py         # Standalone Python version of notebook 3

data_analysis/
├── DADA2_Benchmarking_Comparison.py  # DADA2 vs framework comparison
├── NUMT_Codon_Usage_Analysis.py      # NUMT validation: RSCU, ENC, PCA analysis
└── dada2_benchmarking/
    ├── dada2_chimera_detection_v2.R  # R script for DADA2 chimera detection (latest)
    └── dada2_chimera_detection.R     # R script v1

figures/
└── Figure_1_Revised.py              # Workflow diagram (Figure 1) generation script
```

---

## Analysis Pipeline

Run notebooks in sequence to reproduce the analysis:

```bash
git clone https://github.com/OunjaiS/ASV-Selection-Taxonomic-Assignment
cd ASV-Selection-Taxonomic-Assignment
pip install pandas numpy scikit-learn matplotlib seaborn biopython scipy jupyter

jupyter notebook
```

| Step | Script | Description |
|------|--------|-------------|
| 1 | `1_sequences_analysis.ipynb` | Sequence QC, codon usage, GC content, motif analysis |
| 2 | `2_classification_analysis.ipynb` | Rule-based + ML classification, MRCT thresholding |
| 3 | `3_statistical_analysis.ipynb` | Success rate analysis, effect sizes, NUMT validation stats |
| 4 | `data_analysis/NUMT_Codon_Usage_Analysis.py` | RSCU, ENC, PCA for NUMT vs Authenticated ASVs |
| 5 | `data_analysis/DADA2_Benchmarking_Comparison.py` | DADA2 pipeline benchmarking |

---

## Classification Categories

| Category | Records | % | Description |
|---|---|---|---|
| Authenticated | 15,901 | 9.04% | Zero phylogenetic distance to reference; highest read abundance |
| Technical Artefacts | 102,083 | 58.02% | Low-abundance sequences; sequencing/PCR errors |
| Environmental Contamination | 24,977 | 14.20% | Non-target taxa (prey, parasites, environmental DNA) |
| Intra-individual Variants (NUMTs) | 19,820 | 11.27% | Co-occurring variants; validated by codon usage bias |
| Cross-sample Contamination | 13,163 | 7.48% | Sequences authenticated in other samples (78.4% lab-sourced) |

---

## NUMT Validation

The NUMT validation (Section 3.6) uses codon usage bias to provide molecular evidence for NUMT classification:
- **Codon bias:** 30 of 62 sense codons significantly different (BH-adjusted p < 0.05)
- **GC content:** Authenticated ASVs 34.72 ± 3.53% vs NUMTs 34.26 ± 3.62% (p = 3.67 × 10⁻²⁷)
- **ENC:** NUMTs show higher codon bias scores (9.68 ± 2.49 vs 9.36 ± 2.31; p = 4.10 × 10⁻²²)
- **PCA:** PC1 (14.6% variance) separates Authenticated and NUMT clusters

Script: `data_analysis/NUMT_Codon_Usage_Analysis.py`

---

## Software Dependencies

- Python 3.8+: `pandas`, `numpy`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`, `biopython`
- R 4.0+: `dada2`, `ShortRead` (for DADA2 chimera detection scripts)
- VSEARCH / UNOISE3 (external — see Supplementary Material 2 on Dryad)
- MAFFT (external — for sequence alignment)
- FastTree2 (external — for phylogenetic tree construction)

---

## Data Availability

| Resource | Location |
|---|---|
| Raw sequences & authentication results | [Dryad](https://doi.org/10.5061/dryad.v41ns1s9q) |
| Analysis scripts | This repository |
| Code archive (Zenodo DOI) | https://doi.org/10.5281/zenodo.XXXXXXX |
| Specimen images | [Flickr](https://www.flickr.com/photos/site-100/) |
| Taxonomic data | [iNaturalist](https://www.inaturalist.org/people/site_100) |

---

## Citation

> Ounjai S, Liu H, Zhou Z, Pestana Correia M, Creedy TJ, Andújar C, Arribas P, Vogler AP.
> *Diversity and Phylogenetic Placement of Amplified Sequence Variants in Single-Specimen Barcoding of Tropical Insects.*
> (in review, 2026)

---

**Keywords:** DNA barcoding, COI, ASV authentication, NUMTs, DADA2, phylogenetic placement, tropical Coleoptera, codon usage bias
