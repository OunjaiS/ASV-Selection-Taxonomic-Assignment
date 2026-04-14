# Reproducibility Code
**Phylogenetic Authentication of Amplicon Sequence Variants in Single-Specimen Metabarcoding of Tropical Insects**  
*Molecular Ecology Resources* — in press (2026)

---

## Directory Structure

```
data_analysis/
├── config.py                    ← Python path config (auto-detects all paths)
├── config.R                     ← R path config (auto-detects all paths)
├── requirements.txt             ← Python dependencies (pip install -r requirements.txt)
│
├── input/
│   ├── raw_data/                ← Large input files (download from Dryad)
│   │   ├── ASV_Authentication_Results_030925.csv
│   │   ├── ASV_sequences_64544.fasta
│   │   └── ASV_sequences_66595.fasta
│   └── reference/               ← Small reference files (included in repo)
│       ├── Table_S_NUMT_Anchor_Pairs_n930.csv
│       ├── proven_numt_family_concordant_ids.csv   ← required by script 4
│       └── proven_numt_930_ids.csv                 ← required by script 4
│
├── scripts/                     ← All analysis scripts (run in numbered order)
│   ├── 1_sequence_qc.ipynb
│   ├── 2_asv_classification.ipynb
│   ├── 3_statistical_analysis.ipynb
│   ├── 4_numt_codon_usage.py
│   ├── 5_lmm_analyses.R
│   └── 6_dada2_benchmarking.py
│
└── output/                      ← Results (numbered to match scripts)
    ├── 1_sequences_analysis/    ← QC metrics, corrected sequences (large → Dryad)
    ├── 2_classification_analysis/ ← ASV classification results (large → Dryad)
    ├── 3_statistical_analysis/  ← Tables S1–S5 (small → in repo)
    ├── 4_numt_analysis/         ← NUMT validation tables, Figure 7 draft
    ├── 5_lmm_analysis/          ← LMM results, Table S8, phylogenetic signal
    └── 6_dada2_comparison/      ← DADA2 comparison, Table S6, Figure S2
```

---

## Setup

### 1. Download data from Dryad

Place the following files in `input/raw_data/` before running any scripts:

| File | Description | Size |
|------|-------------|------|
| `ASV_Authentication_Results_030925.csv` | Main ASV records with classification labels (175,954 rows) | ~127 MB |
| `ASV_sequences_64544.fasta` | FASTA sequences for 64,544 unique ASVs | ~25 MB |
| `ASV_sequences_66595.fasta` | Full FASTA sequences (used in script 1 only) | ~27 MB |

Place these files in `input/reference/` (also available from Dryad):

| File | Description |
|------|-------------|
| `Table_S_NUMT_Anchor_Pairs_n930.csv` | Pre-computed NUMT–anchor pair table (n=930 proven NUMTs) |
| `proven_numt_family_concordant_ids.csv` | NUMT IDs passing family-concordance filter (n=1,016) |
| `proven_numt_930_ids.csv` | Final 930 proven NUMT IDs (all four filter stages) |

> Note: `Table_S_NUMT_Anchor_Pairs_n930.csv` is already included in this repo. The two ID files are available from Dryad.

### 2. Python environment

**Option A — conda (recommended):**

```bash
conda create -n asv-auth python=3.11
conda activate asv-auth
pip install -r requirements.txt
```

**Option B — pip only:**

```bash
pip install -r requirements.txt
```

Key Python packages installed by `requirements.txt`:

| Package | Purpose |
|---------|---------|
| `pandas`, `numpy`, `scipy` | Data analysis (all scripts) |
| `scikit-learn` | ML classification (script 2) |
| `matplotlib`, `seaborn` | Figures (scripts 3, 4, 6) |
| `statsmodels` | Statistical tests (script 3) |
| `jupyter`, `nbconvert`, `ipykernel` | Running `.ipynb` notebooks |

### 3. R packages

```r
install.packages(c("lme4", "lmerTest", "ape", "phytools", "rotl"))
```

Tested versions: lme4 1.1.38, lmerTest 3.2.0, ape 5.8, phytools 2.5.2, rotl 3.1.1

> The phylogenetic signal analysis (script 5, Analysis 3) queries the Open Tree of Life API and requires an internet connection.

---

## Running the Scripts

Scripts must be run **in numbered order** — each script's outputs feed into the next.

### Script 1 — Sequence QC

**File:** `scripts/1_sequence_qc.ipynb`  
**Input:** `input/raw_data/ASV_sequences_66595.fasta`, `ASV_Authentication_Results_030925.csv`  
**Output:** `output/1_sequences_analysis/`

```bash
cd scripts
jupyter nbconvert --to notebook --execute 1_sequence_qc.ipynb
```

Key outputs:
- `ASV_Complete_Analysis.csv` — QC metrics for all ASVs
- `ASV_Corrected_Sequences.fasta` — frameshift-corrected sequences
- `Codon_Usage_Table.csv`, `AA_Composition.csv`, `Motif_Analysis_Summary.csv`

---

### Script 2 — ASV Classification

**File:** `scripts/2_asv_classification.ipynb`  
**Input:** `output/1_sequences_analysis/ASV_Complete_Analysis.csv`  
**Output:** `output/2_classification_analysis/`

```bash
jupyter nbconvert --to notebook --execute 2_asv_classification.ipynb
```

Key outputs:
- `ASV_Final_Classification.csv` — rule-based + ML classification for all ASVs
- `Classification_Statistics.csv`, `Feature_Importance.csv`

---

### Script 3 — Statistical Analysis

**File:** `scripts/3_statistical_analysis.ipynb`  
**Input:** `output/2_classification_analysis/ASV_Final_Classification.csv`  
**Output:** `output/3_statistical_analysis/`

```bash
jupyter nbconvert --to notebook --execute 3_statistical_analysis.ipynb
```

Key outputs (manuscript supplementary tables):
- `Table_S1_Country_Success.csv` — authentication success by country
- `Table_S2_Method_Success.csv` — authentication success by collection method
- `Table_S3_Batch_Success.csv` — authentication success by sequencing batch
- `Table_S4_All_Families.csv` — authentication success by beetle family (119 families)
- `Table_S5_Statistical_Tests.csv` — chi-square and Mann-Whitney test results
- `Classification_Summary.csv` — ASV count per classification category
- `Manuscript_Basic_Statistics_CORRECTED.csv` — key summary statistics cited in manuscript

---

### Script 4 — NUMT Codon Usage Validation

**File:** `scripts/4_numt_codon_usage.py`  
**Input:**
- `input/raw_data/ASV_Authentication_Results_030925.csv`
- `input/raw_data/ASV_sequences_64544.fasta`
- `input/reference/Table_S_NUMT_Anchor_Pairs_n930.csv`
- `input/reference/proven_numt_family_concordant_ids.csv`
- `input/reference/proven_numt_930_ids.csv`

**Output:** `output/4_numt_analysis/`

```bash
python3 scripts/4_numt_codon_usage.py
```

Key outputs:
- `proven_numt_nucleotide_composition.csv` — GC/AT, position-specific GC, ENC
- `proven_numt_rscu_per_codon.csv` — RSCU per codon (→ **Table S7**)
- `proven_numt_per_asv_rscu_matrix.csv` — per-ASV RSCU matrix for PCA
- `proven_numt_anchor_pdist.csv` — NUMT–anchor p-distances
- `proven_NUMT_Validation_Figure7.pdf` — 4-panel Figure 7 draft
- `lmm_numt_all_metrics_data.csv` ← input for script 5
- `lmm_numt_gc_data.csv` ← input for script 5

---

### Script 5 — Within-Specimen LMM Analyses

**File:** `scripts/5_lmm_analyses.R`  
**Input:**
- `output/4_numt_analysis/lmm_numt_all_metrics_data.csv`
- `output/4_numt_analysis/lmm_numt_gc_data.csv`
- `output/5_lmm_analysis/lmm_readcount_data.csv` *(from Dryad)*
- `output/5_lmm_analysis/family_gc_content.csv` *(from Dryad)*

**Output:** `output/5_lmm_analysis/`

```bash
Rscript scripts/5_lmm_analyses.R
```

Key outputs:
- `LMM1_ReadCount_results.txt` — LMM: read count vs authentication status
- `LMM2_NUMT_GC_results.txt` — LMM: GC content within specimen
- `LMM_WithinSpecimen_AllMetrics_Table.csv` — **Table S8**
- `PhyloSignal_GC_results.txt` — Blomberg's K and Pagel's λ

> Requires internet connection for Open Tree of Life API (Analysis 3).

---

### Script 6 — DADA2 Benchmarking

**File:** `scripts/6_dada2_benchmarking.py`  
**Input:** `input/raw_data/ASV_Authentication_Results_030925.csv`  
**Output:** `output/6_dada2_comparison/`

```bash
python3 scripts/6_dada2_benchmarking.py
```

Key outputs:
- `Table_S6_DADA2_analysis.csv` — pipeline feature comparison (→ **Table S6**)
- `authentication_gap_analysis.csv` — what DADA2 alone would miss
- `mrct_threshold_comparison.csv` — MRCT sensitivity analysis
- `Figure_S2_DADA2_Benchmarking_Comparison.pdf` — **Figure S2**

> Note: This script can be run independently at any time after script 1. Direct re-processing with DADA2 is not feasible because raw FASTQ files span multiple independent sequencing runs with heterogeneous parameters (Section 4.5 of manuscript).

---

## How Path Detection Works

All scripts use `config.py` (or `config.R`) located in `data_analysis/`. The config auto-detects the `data_analysis/` directory from its own file path — no manual path editing is needed.

For **Python scripts**, the config is loaded with:
```python
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CSV_PATH, NUMT_OUT, ...
```

For **Jupyter notebooks**, the config is found by searching upward from the CWD:
```python
ANALYSIS_DIR = next(p for p in [Path().resolve(), Path().resolve().parent, ...]
                    if (p / "config.py").exists())
sys.path.insert(0, str(ANALYSIS_DIR))
from config import ...
```

For the **R script**:
```r
source(file.path(dirname(normalizePath(sys.frames()[[1]]$ofile)), "..", "config.R"))
```

If you move the entire `data_analysis/` folder, paths update automatically.

---

## Key Definitions

**Proven NUMT** — A Nuclear_Pseudogene ASV that satisfies all four criteria:
1. Reads ≥ 4 (above MRCT) in at least one specimen
2. Co-occurs with the **same** authenticated anchor ASV in ≥ 2 independent specimens
3. NUMT and anchor morphologically/phylogenetically in the same family (family concordance)
4. Aligned p-distance ≤ 10% to anchor (NCBI Genetic Code Table 5, Needleman-Wunsch)

Final dataset: **n = 930 proven NUMTs** paired with **n = 330 unique anchor ASVs**.

**MRCT** — Minimum Read Count Threshold = 4 reads (equivalent to `--minsize 4` in VSEARCH UNOISE3).

**Genetic code** — All translations use **NCBI Table 5 (Invertebrate Mitochondrial)**: TGA → Trp, TAA/TAG → Stop, ATA → Met, AGA/AGG → Ser. This gives **62 sense codons** (not 61).

---

## Software Versions

| Software | Version |
|----------|---------|
| Python | 3.11 |
| R | ≥ 4.2 |
| pandas | ≥ 1.5 |
| numpy | ≥ 1.23 |
| scikit-learn | ≥ 1.1 |
| scipy | ≥ 1.9 |
| matplotlib | ≥ 3.6 |
| seaborn | ≥ 0.12 |
| statsmodels | ≥ 0.14 |
| lme4 | ≥ 1.1-31 |
| lmerTest | ≥ 3.1-3 |
| rotl | ≥ 3.0 |

---

## Citation

If you use these scripts, please cite:

> Ounjai, S., et al. Diversity and Phylogenetic Placement of Amplicon Sequence Variants in Single-Specimen Metabarcoding of Tropical Insects. *Molecular Ecology Resources* (in press, 2026).

Raw data and large intermediate files are deposited on Dryad. See the manuscript Data Availability section for accession details.
