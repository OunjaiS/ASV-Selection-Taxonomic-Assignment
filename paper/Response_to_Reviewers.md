# Response to Reviewers

## Manuscript: "Diversity and Phylogenetic Placement of Amplified Sequence Variants in Single-Specimen Barcoding of Tropical Insects"

We thank both reviewers for their constructive and insightful comments. Below we address each comment point-by-point. Reviewer comments are in **bold**, and our responses follow in regular text. Page and line numbers refer to the revised manuscript.

---

## REVIEWER 1

### Major Comments

**1. NUMT Validation (R1 Major): The classification of intra-individual variants as NUMTs relies on co-occurrence with authentic ASVs and phylogenetic proximity, but lacks direct molecular evidence (e.g., genomic mapping or PCR validation of nuclear localization).**

We appreciate this important comment. In the revised manuscript we have provided three additional lines of evidence supporting the NUMT classification:

(a) **Codon usage bias analysis:** We performed a comprehensive comparison of codon usage patterns between Authenticated ASVs (n = 8,533 unique) and putative NUMTs (n = 26,416 unique). The analysis revealed statistically significant differences across multiple metrics:
- **Nucleotide composition:** Authenticated ASVs showed significantly higher GC content (34.72 ± 3.53%) than NUMTs (34.26 ± 3.62%; Mann-Whitney U, p = 3.67 × 10⁻²⁷), with differences at the 1st (p = 1.54 × 10⁻³⁸) and 3rd (p = 7.54 × 10⁻¹⁷) codon positions, but not at the functionally constrained 2nd position (p = 0.159).
- **RSCU patterns:** 30 of 62 sense codons showed significantly different usage (Benjamini-Hochberg adjusted p < 0.05). PCA on per-ASV RSCU values separated Authenticated and NUMT clusters along PC1 (14.6% variance).
- **Codon bias score:** NUMTs showed higher codon bias scores (9.68 ± 2.49) than Authenticated ASVs (9.36 ± 2.31; p = 4.10 × 10⁻²²), consistent with relaxed selection pressure.
- **Amino acid composition:** Authenticated ASVs had higher hydrophobic amino acid content (65.54 ± 1.65%) than NUMTs (65.31 ± 2.45%; p = 1.13 × 10⁻⁵), consistent with functional constraints on transmembrane COI domains.

See new **Figure 7** and **Table 3** (now in main manuscript text, Results Section 3.6).

(b) **Case studies:** We present three representative NUMTs from different beetle families:
1. **Curculionidae (French Guiana):** uniq627 — 0.10% reads, phylogenetic distance 2.04, GC 43.9% vs Main ASV 33.1%, ENC 31.1 vs 25.1, detected in 4/16 specimens (validated in 3).
2. **Staphylinidae (French Guiana):** uniq627 — 0.08% reads, phylogenetic distance 1.91, GC 43.9% vs Main ASV 36.9%, ENC 31.1 vs 26.9, detected in 4/16 specimens (validated in 3).
3. **Cerambycidae (Ecuador):** MSL4763 — 5.17% reads, phylogenetic distance 1.68, GC 41.0% vs Main ASV 34.1%, ENC 29.1 vs 23.0, detected in 2/2 specimens.

Each case shows convergent evidence from abundance, phylogenetic placement, and codon usage supporting NUMT classification.

(c) **Limitation acknowledgement:** We have added text in Discussion 4.3 acknowledging that definitive confirmation of nuclear localisation requires genomic approaches (long-read sequencing, chromosome mapping, or PCR with nuclear-specific flanking primers), and noting that the convergence of multiple independent lines of evidence provides strong support for the classification.

**Status: Completed.** See tracked changes: Change #13 in REVISION_TRACKED_CHANGES.md.

---

**2. Benchmarking with Alternative Methods (R1 Major): The manuscript does not explicitly compare the proposed framework with alternative ASV filtering tools (e.g., DADA2, Deblur).**

We agree that a benchmarking comparison would strengthen the manuscript. We have conducted a comprehensive comparison between our authentication framework and the DADA2 pipeline:

**(a) Methodological comparison:** We provide a detailed feature-by-feature comparison table (Supplementary Table S6) showing that DADA2 addresses 5 of 15 pipeline features (denoising, error modelling, chimera detection, quality filtering, abundance filtering), while our framework addresses all 15, including specimen-level authentication, phylogenetic placement, NUMT detection, cross-contamination detection, and ML-based confidence scoring.

**(b) Abundance threshold analysis:** We compared the performance of different Minimum Read Count Thresholds (MRCT). At MRCT=4 (our pipeline), sensitivity for Authenticated ASVs is 100% (no authenticated sequences lost) while all 83,709 Technical Artifacts are removed. At MRCT=2 (DADA2's typical effective minimum), 27,690 Technical Artifacts would still be retained.

**(c) Authentication gap analysis:** The critical finding is that DADA2-like denoising alone (reads ≥ 2) would retain 119,092 records, of which 103,171 (86.6%) are non-authenticated sequences:
- 39,267 Environmental Contamination records
- 36,214 Nuclear Pseudogene (NUMT) records
- 27,690 Technical Artifact records (removable by a stricter abundance threshold)

Even with MRCT=4 (removing all Technical Artifacts), 75,481 problematic records (NUMTs + Environmental Contamination) — representing biological contaminants that pass all quality filters — would remain undetected without the specimen-level authentication framework.

**(d) Chimera detection:** The absence of chimeric sequences in the final ASV table confirms that our UCHIME3 de novo chimera detection step, applied during the bioinformatics pipeline (Section 2.3), effectively removed chimeric sequences prior to ASV table construction.

**(e) Framework complementarity:** We emphasise that DADA2 and our framework are **complementary**, not competing approaches. DADA2 excels at statistical error correction through its parametric error model, while our framework extends the analysis with specimen-level authentication, phylogenetic placement, and multi-category classification. The authentication framework is pipeline-agnostic and could be applied downstream of either DADA2 or VSEARCH/UNOISE denoising.

See new additions to Methods Section 2.3 and Results Section 3.4, along with the detailed breakdown in new Discussion Section 4.5, Supplementary Figure S2, and Supplementary Table S6.

**Status: Completed.** See tracked changes: Change #14 in REVISION_TRACKED_CHANGES.md.

---

**3. Mechanisms Underlying Environmental Effects (R1 Major): The biological mechanisms driving differences in authentication success across trap types and preservation fluids are not fully discussed.**

We have expanded Discussion Section 4.2 with specific mechanistic explanations:

- **DNA degradation mechanisms:** Added detail on hydrolysis, oxidative damage, and depurination that occur during prolonged trap exposure in tropical environments (revised manuscript, Section 4.2, paragraph 2).
- **Preservation fluid comparison:** Added mechanistic explanation of why ethanol outperforms SDS+EDTA (rapid tissue dehydration and nuclease inactivation vs. residual enzymatic activity at tropical temperatures).
- **Pitfall trap mechanisms:** Expanded discussion of humic acid PCR inhibition, UV radiation effects, and exogenous DNA from non-target organisms (Schrader et al. 2012).

See tracked changes: Changes #1–2 in REVISION_TRACKED_CHANGES.md.

---

### Minor Comments

**4. Cross-contamination Mitigation Strategies (R1 Minor): The discussion should briefly address strategies to mitigate cross-contamination in future studies.**

We have added a new paragraph in Discussion Section 4.3 covering:
- Field-level mitigation (minimising co-exposure time, individual preservation vials)
- Laboratory protocols (physical separation, UV decontamination)
- Unique dual indexing (UDI) to eliminate index hopping (Costello et al. 2018; MacConaill et al. 2018)
- Negative controls (extraction blanks and PCR no-template controls)
- Computational approaches including occupancy modelling (Ficetola et al. 2015)

See tracked changes: Change #3 in REVISION_TRACKED_CHANGES.md.

---

## REVIEWER 2

### Major Comments

**5. Illumina vs PacBio Platform Comparison (R2 Major): A broader comparison between Illumina versus PacBio HiFi sequencing is missing. With 58% assigned to technical artifacts, this indicates Illumina may be less suitable for megabarcoding.**

We have added a comprehensive platform comparison paragraph in Discussion Section 4.4 addressing:

- **Clarification of the 58% artifact rate:** We clarify that this high proportion reflects the deep sequencing strategy rather than Illumina platform limitations. 92.1% of these artefacts were low-abundance sequences (1–3 reads) effectively removed by the MRCT=4 threshold.
- **Platform comparison:** We compare Illumina (cost-effective ASV-level resolution), PacBio HiFi (reference-quality haplotypes, 94.9% Sanger concordance; Congrains et al. 2025), and ONT (maximum throughput, OTU-level resolution).
- **Framework applicability:** We note that the authentication framework is platform-agnostic but benefits most from high per-base accuracy platforms (Illumina, PacBio).

See tracked changes: Change #4 in REVISION_TRACKED_CHANGES.md.

---

**6. Primer Reference Issue (R2 Major): The forward primer does not match any primer in Arribas et al. 2016. Was it created for this study?**

The reviewer is correct that only the reverse primer matches the primers used in Arribas et al. (2016). We have corrected the primer references in Methods Section 2.2:

- The reverse primer (5'–TANACYTCNGGRTGNCCRAARAAYCA–3') is a degenerate variant of jgHCO2198 (Geller et al. 2013), with inosine bases replaced by fully degenerate positions.
- The forward primer (5'–CCNGAYATRGCNTTYCCNCG–3'): [Authors to verify from laboratory records — see two options provided in REVISION_TRACKED_CHANGES.md, Change #5]
- Added justification for selecting this primer pair for Coleoptera-focused study.

See tracked changes: Change #5 in REVISION_TRACKED_CHANGES.md.

---

### Minor Comments

**7. Manuscript Length (R2 Minor): Section 3.5 could be shortened.**

We have retained Section 3.5 in full, as it is within the journal's word limit and provides valuable biological context about factors affecting authentication success that would be difficult to convey through supplementary material alone.

---

**8. Figure 1 Formatting (R2 Minor): Text size, alignment, legend, and box numbering need improvement.**

We have completely redesigned Figure 1:
- Box numbers now correspond to Methods sections 2.1–2.6
- Increased text size (minimum 9pt)
- Compact horizontal legend
- Improved text alignment and arrow clarity
- Updated caption to reference Methods section numbers

See tracked changes: Changes #6–7 in REVISION_TRACKED_CHANGES.md.

---

**9. Broken Link (R2 Minor): Line 169 gives a 404 error.**

Fixed. The incomplete URL has been replaced with the full GitHub repository address: https://github.com/OunjaiS/ASV-Selection-Taxonomic-Assignment

See tracked changes: Change #8 in REVISION_TRACKED_CHANGES.md.

---

**10. Reproducibility — Script Versioning and Repository (R2 Minor): All scripts should be referenced with versions and made available in a single repository with DOI.**

We have:
- Added specific version numbers for all software tools in Methods sections 2.3, 2.4, and 2.6
- Added tool names and repository URLs for third-party scripts (tjcreedy/biotools, catfasta2phyml)
- Added access dates for all third-party resources
- Created a Zenodo DOI for the GitHub repository [DOI to be inserted before submission]

See tracked changes: Changes #9–12 in REVISION_TRACKED_CHANGES.md.

---

## Summary of Changes

| Comment | Type | Status | Action Taken |
|---------|------|--------|-------------|
| #1 NUMT Validation | R1 Major | **Completed** | Codon usage bias analysis, 3 case studies, limitation statement |
| #2 DADA2 Benchmarking | R1 Major | **Completed** | Framework comparison, gap analysis |
| #3 Environmental Mechanisms | R1 Major | **Completed** | Discussion 4.2 expanded |
| #4 Cross-contamination | R1 Minor | **Completed** | New paragraph in Discussion 4.3 |
| #5 Platform Comparison | R2 Major | **Completed** | New paragraph in Discussion 4.4 |
| #6 Primer Reference | R2 Major | **Partially** | Corrected; author verification needed |
| #7 Manuscript Length | R2 Minor | **Completed** | No changes (section retained) |
| #8 Figure 1 Formatting | R2 Minor | **Completed** | Figure redesigned |
| #9 Broken Link | R2 Minor | **Completed** | URL fixed |
| #10 Script Versioning | R2 Minor | **Completed** | Versions added, Zenodo DOI pending |
