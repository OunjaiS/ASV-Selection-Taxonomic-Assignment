# GitHub & Zenodo Upload Guide
**Project:** ASV Authentication Manuscript
**Date:** 2026-02-22
**GitHub repo URL:** https://github.com/OunjaiS/ASV-Selection-Taxonomic-Assignment

---

## สคริปต์ที่ต้องอัปขึ้น GitHub

### 📁 `script/` — Core Analysis Scripts (สำคัญที่สุด)

| ไฟล์ | คำอธิบาย |
|---|---|
| `1_sequences_analysis.ipynb` | ASV sequence quality analysis, nucleotide composition, codon usage |
| `2_classification_analysis.ipynb` | ASV classification pipeline, feature importance, authentication stats |
| `3_statistical_analysis.ipynb` | Statistical tests (Mann-Whitney U, effect sizes, chi-squared) + outputs |
| `3_statistical_analysis.py` | Python script version ของ notebook 3 |
| `4_manuscript_figures.ipynb` | Figure generation for all manuscript figures (F1-F7) |

### 📁 `data_analysis/` — Analysis Modules

| ไฟล์ | คำอธิบาย |
|---|---|
| `DADA2_Benchmarking_Comparison.py` | DADA2 vs framework benchmarking comparison |
| `NUMT_Codon_Usage_Analysis.py` | NUMT validation via codon usage bias (RSCU, ENC, PCA) |
| `dada2_benchmarking/dada2_chimera_detection.R` | R script for DADA2 chimera detection |
| `dada2_benchmarking/dada2_chimera_detection_v2.R` | R script v2 (latest version) |

### 📁 `figures/` — Figure Generation

| ไฟล์ | คำอธิบาย |
|---|---|
| `Figure_1_Revised.py` | Workflow diagram Figure 1 |
| `Figure_7_NUMT_Validation.png` | Final Figure 7 (NUMT codon bias) |
| `Figure_7_NUMT_Validation.pdf` | Figure 7 PDF version |

### 📁 `suplementary_data/` — Supplementary Tables

| ไฟล์ | คำอธิบาย |
|---|---|
| `Table_S6_DADA2_analysis.csv` | DADA2 benchmarking supplementary (Feature comparison, MRCT, Auth gap) |
| `Table_S7_RSCU_comparison.csv` | Full RSCU codon comparison (62 codons) Auth vs NUMT |

---

## ขั้นตอน: อัปโหลด GitHub

### ขั้นตอนที่ 1: Initialize Git Repository

```bash
cd /Users/sarawut/Desktop/Manuscript_ASV_selection
git init
git add .
git commit -m "Initial commit: ASV authentication analysis scripts and manuscript files"
```

### ขั้นตอนที่ 2: เชื่อมกับ Remote Repository

```bash
git remote add origin https://github.com/OunjaiS/ASV-Selection-Taxonomic-Assignment.git
git branch -M main
git push -u origin main
```

> ⚠️ ถ้า repo นั้นมีไฟล์อยู่แล้ว ให้ pull ก่อน:
> ```bash
> git pull origin main --allow-unrelated-histories
> ```

### ขั้นตอนที่ 3: สร้าง Release บน GitHub

1. ไปที่ GitHub repo → **Releases** → **"Create a new release"**
2. Tag version: `v1.0.0`
3. Title: `"ASV Authentication Framework v1.0 — Initial Release"`
4. Description: ระบุ version ที่ตรงกับ manuscript submission

---

## ขั้นตอน: รับ Zenodo DOI

### ขั้นตอนที่ 1: เชื่อม Zenodo กับ GitHub

1. ไปที่ [https://zenodo.org](https://zenodo.org) → Login ด้วย GitHub account
2. ไปที่ **GitHub** tab (มุมขวาบน)
3. เปิด (Flip ON) สำหรับ repo: **OunjaiS/ASV-Selection-Taxonomic-Assignment**

### ขั้นตอนที่ 2: สร้าง GitHub Release

1. กลับไปที่ GitHub → ทำ **New Release** (ถ้ายังไม่ได้ทำ)
2. Zenodo จะ **automatically** สร้าง DOI ให้ทันทีที่มี Release ใหม่

### ขั้นตอนที่ 3: รับ DOI

1. กลับไปที่ Zenodo → **Uploads** → จะเห็น record ใหม่
2. Copy DOI ในรูปแบบ: `https://doi.org/10.5281/zenodo.XXXXXXX`

### ขั้นตอนที่ 4: อัปเดต Manuscript

แก้ใน Word/PDF:
```
เก่า: https://doi.org/10.5281/zenodo.XXXXXXX
ใหม่: https://doi.org/10.5281/zenodo.[DOI จริง]
```

---

## ไฟล์ที่ **ไม่ควร** อัปขึ้น GitHub

| ไฟล์/โฟลเดอร์ | เหตุผล |
|---|---|
| `raw_data/` | ไฟล์ขนาดใหญ่ → อัปไว้ที่ Dryad แทน |
| `script/__pycache__/` | auto-generated ไม่จำเป็น |
| `.DS_Store` | macOS system file |
| `.gemini/`, `.claude/` | AI assistant working files |
| `.Rhistory` | R session history |

(สิ่งเหล่านี้อยู่ใน `.gitignore` แล้ว)

---

## Zenodo Metadata ที่ต้องกรอก

เมื่อสร้าง record บน Zenodo ให้กรอก:

| Field | ค่า |
|---|---|
| **Title** | ASV Authentication Framework for Single-Specimen Barcoding |
| **Authors** | Sarawut Ounjai, Huaxi Liu, Zichen Zhou, Maria Pestana Correia, Thomas J. Creedy, Carmelo Andújar, Paula Arribas, Alfried P. Vogler |
| **Description** | Authentication pipeline and analysis scripts for "Diversity and Phylogenetic Placement of Amplified Sequence Variants in Single-Specimen Barcoding of Tropical Insects" |
| **License** | MIT หรือ CC BY 4.0 (แล้วแต่ตกลงกับ co-authors) |
| **Related Identifiers** | DOI ของ Paper จาก Journal (ถ้ามีแล้ว) |
| **Keywords** | DNA barcoding, ASV authentication, COI, Coleoptera, tropical insects, NUMTs, DADA2 |
