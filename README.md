# ASV Authentication for COI DNA Barcoding in Insects

This repository contains all data, analysis scripts, and supplementary materials for the manuscript on Amplicon Sequence Variant (ASV) authentication in COI DNA barcoding of insects.

## Overview

This project develops and validates methods to distinguish authentic COI sequences from technical artifacts, contamination, and pseudogenes in high-throughput DNA barcoding datasets. Using a combination of rule-based classification and machine learning approaches, we analyzed 175,954 ASV records from 18,533 insect specimens, achieving an 86.52% authentication success rate.

## Dataset Summary

- **Total ASV Records**: 175,954
- **Unique ASVs**: 64,544
- **Total Specimens**: 18,533
- **Authenticated Specimens**: 14,715 (86.52%)
- **Total Sequencing Reads**: 36,459,895
- **Primary Taxonomic Group**: Coleoptera (beetles)
- **Geographic Focus**: Primarily French Guiana (South America)

## Repository Structure

```
Manuscript_ASV_selection/
├── raw_data/                   # Original sequence and authentication data
├── script/                     # Analysis scripts (Jupyter notebooks)
├── data_analysis/             # Processed analysis results
│   ├── sequences_analysis/    # Sequence quality metrics and composition
│   └── classification_analysis/ # Classification results and ML model
├── suplementary_data/         # Supplementary tables and materials
├── figures/                   # Publication-ready figures (PNG and PDF)
├── paper/                     # Manuscript document
└── README.md                  # This file
```

## Data Files

### Raw Data (`raw_data/`)

| File | Description | Records/Lines |
|------|-------------|---------------|
| `ASV_Authentication_Results_030925.csv` | Main authentication results with metadata | 175,954 records |
| `ASV_sequences_64544.fasta` | Unique ASV sequences (FASTA format) | 64,544 sequences |
| `ASV_sequences_64544.txt` | Text version of sequences | - |

### Analysis Scripts (`script/`)

Four Jupyter notebooks implement the complete analysis pipeline:

1. **`1_sequences_analysis.ipynb`** - Sequence Quality and Composition Analysis
   - DNA/protein composition
   - Codon usage patterns
   - Motif discovery
   - GC content and phylogenetic features
   - Translation validation

2. **`2_classification_analysis.ipynb`** - Classification Pipeline
   - Rule-based classification
   - Machine learning (Random Forest) classification
   - Feature importance analysis
   - Classification validation

3. **`3_statistical_analysis.ipynb`** - Statistical Testing and Validation
   - Success rate analysis by geographic location
   - Success rate analysis by collection method
   - Batch-wise comparisons
   - Confidence calculations

4. **`4_manuscript_figures.ipynb`** - Publication Figure Generation
   - Dataset overview visualizations
   - Classification breakdown plots
   - Success factor analysis
   - Workflow diagrams

### Processed Data (`data_analysis/`)

#### Sequences Analysis (`data_analysis/sequences_analysis/`)

| File | Description |
|------|-------------|
| `ASV_Corrected_Sequences.fasta` | Corrected and validated sequences |
| `ASV_Complete_Analysis.csv` | Comprehensive sequence metrics |
| `Analysis_Summary.csv` | Summary statistics (100% QC pass rate) |
| `AA_Composition.csv` | Amino acid composition statistics |
| `Codon_Usage_Table.csv` | Codon usage patterns |
| `Motif_Analysis_Summary.csv` | Conserved motif analysis |
| `visualizations/Sequence_Characteristics.png` | Quality visualizations |
| `visualizations/Quality_Distribution.png` | Distribution plots |

#### Classification Analysis (`data_analysis/classification_analysis/`)

| File | Description |
|------|-------------|
| `ASV_Final_Classification.csv` | Final classification results (175,954 records) |
| `Classification_Statistics.csv` | Classification breakdown by category |
| `Classification_Report.txt` | Detailed methodology and statistics |
| `Classification_Model.pkl` | Trained Random Forest model |
| `ASV_Sequence_Summary.csv` | Per-sequence classification summary |
| `Feature_Importance.csv` | ML feature importance rankings |

### Supplementary Data (`suplementary_data/`)

| File | Description |
|------|-------------|
| `Table_S1_Classification_Breakdown.csv` | Classification categories and percentages |
| `Table_S2_Country_Success.csv` | Success rates by country |
| `Table_S3_Method_Success.csv` | Success rates by collection method |
| `Table_S4_Batch_Success.csv` | Batch-wise success metrics |
| `Table_S5_All_Families.csv` | Taxonomy family breakdown (121 families) |
| `Table_S6_Statistical_Tests.csv` | Statistical test results |
| `metadata.csv` | Sample metadata (175,954 records) |
| `Manuscript_Basic_Statistics.csv` | Key manuscript metrics |
| `Supplementary_Material_1_Primers.xlsx` | Primer sequences used |
| `Supplementary_Material_2_Bioinformatics_Pipeline.docx` | Detailed pipeline documentation |
| `Supplementary_Material_3_mitogenomes.xlsx` | Mitogenome reference data |
| `Supplementary_Material_4_Phylogenetic_Framework.docx` | Phylogenetic methodology |

### Figures (`figures/`)

All figures available in both PNG (screen) and PDF (print) formats:

- **Figure 1**: Dataset Overview
- **Figure 2**: MRCT (Most Recent Common Time) Analysis
- **Figure 3**: Classification Breakdown
- **Figure 4**: Success Factors
- **Figure 5**: Feature Importance (Machine Learning)
- **Figure 6**: Analysis Workflow

## Classification Categories

The pipeline classifies ASVs into five categories:

| Category | Percentage | Avg Reads | Description |
|----------|-----------|-----------|-------------|
| **Authenticated** | 9.04% | 1,903.8 | High-quality, validated sequences matching specimen identification |
| **Cross-Contamination** | 7.48% | 100.4 | Sequences authenticated in other samples but not in current sample |
| **Intra-Species Variants** | 11.27% | 59.8 | Co-occurring variants (potential NUMTs/heteroplasmy) |
| **Environmental Contamination** | 14.20% | 138.5 | Non-target sequences (prey, parasites, environmental DNA) |
| **Technical Artifacts** | 58.02% | 2.2 | Low-abundance sequences, sequencing errors |

## Requirements

### Software Dependencies

- Python 3.8+
- Jupyter Notebook or JupyterLab
- Required Python packages:
  - pandas
  - numpy
  - scikit-learn
  - matplotlib
  - seaborn
  - biopython

### Installation

```bash
# Clone the repository
git clone https://github.com/OunjaiS/ASV-Selection-Taxonomic-Assignment

# Create virtual environment (optional but recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install pandas numpy scikit-learn matplotlib seaborn biopython jupyter
```

## Usage

### Reproducing the Analysis

Run the Jupyter notebooks in sequence:

```bash
# Start Jupyter
jupyter notebook

# Then run notebooks in order:
# 1. script/1_sequences_analysis.ipynb
# 2. script/2_classification_analysis.ipynb
# 3. script/3_statistical_analysis.ipynb
# 4. script/4_manuscript_figures.ipynb
```

### Using Pre-computed Results

All processed data and figures are provided in the `data_analysis/` and `figures/` directories. You can directly use these results without re-running the analysis.

### Classification Model

The trained Random Forest model is available at `data_analysis/classification_analysis/Classification_Model.pkl` and can be loaded using:

```python
import pickle
with open('data_analysis/classification_analysis/Classification_Model.pkl', 'rb') as f:
    model = pickle.load(f)
```

## Quality Metrics

### Sequence Quality
- **QC Pass Rate**: 100% (175,112 sequences)
- **COI Confidence**: 98.85% high confidence
- **Quality Grades**: 82.41% A+, 11.75% A, 5.84% B or lower

### Authentication Thresholds
- **Sequence length**: 200-900 bp
- **Minimum reads**: 4
- **Phylogenetic distance** (same species): < 0.02
- **Phylogenetic distance** (same genus): < 0.05
- **Intra-species maximum**: 0.15

## Citation

If you use this data or methods in your research, please cite:

```
xxx
```

## Data Availability

- **GitHub Repository**: https://github.com/OunjaiS/ASV-Selection-Taxonomic-Assignment
- **Dryad Dataset**: https://doi.org/[Dryad DOI:xxx]

## License

You are free to:
- Share — copy and redistribute the material
- Adapt — remix, transform, and build upon the material

Under the following terms:
- Attribution — You must give appropriate credit and indicate if changes were made

## Contact

For questions or issues regarding this dataset, please:
- Open an issue on GitHub
- Contact: ounjai.sa@gmail.com


## Version History

- **v1.0** (2025-XX-XX): Initial release with complete dataset and analysis pipeline

---

**Keywords**: DNA barcoding, COI, Amplicon Sequence Variants, ASV, authentication, insects, Coleoptera, machine learning, quality control, phylogenetics
