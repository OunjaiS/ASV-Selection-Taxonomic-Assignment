# Diversity and Phylogenetic Placement of Amplified Sequence Variants in Single-Specimen Barcoding of Tropical Insects

Sarawut Ounjai^1,2^, Huaxi Liu^1^, Zichen Zhou^1,2^, Maria Pestana Correia^1^, Thomas J. Creedy^1^, Carmelo Andujar^1,3^, Paula Arribas^1,3^, Alfried P. Vogler^1,2^*

^1^ Department of Life Sciences, Natural History Museum, London, United Kingdom
^2^ Department of Life Science, Imperial College London, Ascot, United Kingdom
^3^ Instituto de Productos Naturales y Agrobiologia (IPNA-CSIC), San Cristobal de La Laguna, Tenerife, Spain

\*Corresponding author: a.vogler@nhm.ac.uk

Running title: ASV Authentication in Tropical Beetles

## Abstract

**Background:**
Deep amplicon sequencing for DNA barcoding of individual specimens ("megabarcoding") generates authentic mitochondrial sequences, mixed with nuclear pseudogenes (NUMTs), environmental contaminants, and sequencing artefacts. Accurate classification of these types is essential for reliable taxonomic identifications and for studies of genetic diversity, genome evolution, and ecological interactions.

**Methods:**
Bulk samples of tropical beetles (Coleoptera) were imaged for parataxonomic species assignment and family-level identification, prior to individual Illumina barcoding. Resulting Amplicon Sequence Variants (ASVs) were evaluated through an integrated authentication framework that combines abundance filtering and phylogenetic placement to classify each ASV as either authentic or as belonging to one of four error categories.

**Results:**
We barcoded representative specimens for 18,533 morphospecies from 119 families and 14 countries, yielding 175,954 COI ASVs. Authentication succeeded for 86.5% of quality-passing specimens (15,901 ASVs). The remaining sequences were assigned to technical artefacts (58.0%), environmental contamination including prey DNA (14.2%), intra-individual variants (NUMTs, heteroplasmy; 11.3%), and cross-contamination (7.5%). Authentication success and the proportions of failure categories varied markedly across trap types, sampling campaigns, taxonomic groups, and sequencing runs.

**Conclusions:**
The authenticated ASV set and corresponding high-resolution images have been deposited in public databases, adding up to 13,000 previously unsequenced species of beetles. The authentication framework achieves Sanger-level accuracy at metabarcoding throughput, effectively filtering non-authentic reads that escape conventional abundance-based thresholds. The approach raises the accuracy of molecular ASV data to a standard sufficient for inclusion in barcode reference databases currently lacking representation of tropical insects.

**Keywords:** biodiversity assessment, molecular taxonomy, phylogenetic placement, ASV authentication, mitochondrial metagenomics, tropical diversity

---

## 1. Introduction

DNA barcoding has been transformed by high-throughput sequencing (HTS) technologies, which enable the analysis of PCR amplicon mixtures from multiple specimens, replacing traditional Sanger sequencing, limited to single templates (Hebert et al. 2003; Taberlet et al. 2012). While HTS of amplicons is widely adopted for metabarcoding of ecological communities, it is also increasingly used for genotyping of individual specimens, sometimes termed "megabarcoding" (Chua et al. 2023). Using the mitochondrial COI marker, single-specimen tagging on HTS platforms enables the multiplexing of thousands (Srivathsan et al. 2019) or even 100,000 (Hebert et al. 2025) specimens in a single sequencing run. The addition of sequence tags during PCR, multiplied by tags in the library construction, has thus become a cost-effective method for large-scale DNA-based species identification in highly diverse samples of invertebrates (Chua et al. 2023; Liu et al. 2017; Meier et al. 2016).

HTS (meta)barcoding data are generally analysed by grouping sequence reads into Operational Taxonomic Units (OTUs) based on nucleotide similarity, using a variety of clustering methods. This approach accommodates the variation of amplicons usually ascribed to sequencing or PCR errors that survived the various filtering procedures applied to the raw data. However, as platform accuracy has improved especially for Illumina-type sequencers, reads can now be resolved at the level of individual genotypes as Amplicon Sequence Variants (ASVs) that presumably represent the true DNA templates in a sample (Callahan, McMurdie, and Holmes 2017). Deep amplicon sequencing on Illumina platforms also recovers a wide range of variants of different origins, including traces of non-targets (e.g. prey, symbionts, or parasites), low-copy nuclear pseudogenes ('nuclear mitochondrial', NUMTs), and heteroplasmic mitochondrial variants (Song et al. 2008). In addition, amplification may be affected by cross-contamination among samples, in particular when mixtures from bulk traps affected by varying levels of degradation are sequenced in the same run. In total, deep sequencing of amplicons may yield surprisingly large numbers of extraneous sequence reads, whose combined count may even exceed that of the target haplotype, as observed in Illumina barcoding of single bee specimens (Creedy et al. 2020).

Most "megabarcoding" studies still rely on read clustering, in part necessitated by the use of higher-error Oxford Nanopore Technology (ONT) platforms, which require the assembly of a consensus from multiple reads and thus preclude strict ASV calling (Hebert et al. 2025; Srivathsan et al. 2024). This approach is partly driven by the scientific objectives of such studies, which focus on high specimen throughput, e.g., for community-level abundance studies, and matching to known reference sequences (Hebert et al. 2025). However, the methodology is less powerful for applications requiring exact sequence reads, e.g., those focused on genetic variation in community genetics (Noguerales et al. 2023) or the generation of new reference sequences for addition to the growing barcode databases (Srivathsan et al. 2024). Finally, even when clustering is used, the presence of multiple distantly related sequences in an amplicon mixture requires decisions about filtering strategies, which are either not reported (Jabot et al. 2025; Srivathsan et al. 2021) or based simply on the top abundant variant (Hebert et al. 2025; Meier et al. 2016). Relying solely on the dominant mitochondrial haplotypes may overlook potential cross-contamination of bulk-trapped specimens in the trap, the laboratory, or in the library tagging step. In addition, removing these sequences from consideration also represents a loss of valuable biological information, e.g., for studying trophic interactions and ecological associations.

The prevalence of such secondary reads varies among samples, depending on multiple factors that influence ASV diversity and error rates. DNA degradation is common in biodiversity samples used for barcoding and metabarcoding, particularly from tropical field sites gathered under suboptimal environmental conditions, which can compromise DNA quality, increase PCR error rates, bias the proportions of authentic mitochondrial variants, and increase the risk of cross-contamination in the most highly degraded samples. The impact of NUMTs is similarly unpredictable, as whole-genome studies have shown uneven distribution across taxa and chromosomes, yet these effects still need to be studied across a wider spectrum of lineages (He, Ge, and Liang 2025; Hebert, Bock, and Prosser 2023). Compounding these issues are the limited taxonomic knowledge and lack of reference coverage for most insect groups, especially in tropical regions, which hamper both the species identification via database matches and the discrimination of authentic versus non-target sequences (Andujar et al. 2021; Souto-Vilares et al. 2025; Zinger et al. 2025).

To address these issues, we developed an analytical framework that distinguishes authentic mitochondrial sequences from nuclear pseudogenes, technical artefacts, and contamination, while producing accurate taxonomic assignments even when the coverage of the reference sequence is limited. Our workflow combines the widely applied abundance-based criteria for selecting authentic ASVs with a phylogenetic characterisation of ASVs as either distantly related non-targets (e.g., ingested prey and internal parasites) or closely related NUMTs. The approach requires a phylogenetic tree against which ASVs are assessed (see (Czech, Barbera, and Stamatakis 2019)), which was obtained here by alignment of ASVs to a densely sampled set of mitochondrial genomes. The approach was applied to a set of >18,000 beetle (Coleoptera) specimens, each thought to represent a different morphospecies from poorly characterised tropical forest communities, which were amplicon sequenced with individual tags on an Illumina platform. The authentication of reads required the correct assignment of ASVs at the family level. There are some 200 families in the Coleoptera, most of which can be recognised from specimen images and thus can be compared to the phylogenetic position of corresponding sequence reads. We assembled a densely sampled tree from full-length mitogenome sequences of Coleoptera (Carpenter and Vogler 2025) (Creedy et al. 2025 (in press); Creedy et al. 2025 (in review)), which was used to characterise the ASV mixtures obtained for each specimen and to determine authentic mitochondrial haplotypes. The success of ASV authentication was compared across different sample batches, trap types, preservation conditions, biogeographic realms, and major lineages of Coleoptera. The results highlight the need for stringent characterisation of single-specimen HTS barcode mixtures, especially if they are used as de novo-generated reference sequences missing from existing reference databases.

---

## 2. Methods

The workflow for generating authenticated COI barcodes from mass-trapped tropical insects comprised a multi-step protocol from wet-laboratory amplicon generation to ASV selection, ultimately characterising all ASVs as either authentic barcodes for taxonomic assignment or secondary ASVs representing secondary variants (NUMTs/heteroplasmy), external DNA, cross-sample contaminants, and technical artefacts (Figure 1).

Figure 1. Overview of the integrated ASV analysis workflow. Box numbers correspond to Methods sections 2.1-2.6. The pipeline proceeds from specimen collection and processing (2.1) through molecular data generation (2.2) to the bioinformatics pipeline for ASV calling (2.3). In parallel, a phylogenetic reference framework is constructed from mitogenome sequences (2.4). The core authentication step (2.5) integrates ASV abundance data, phylogenetic placement, and taxonomic congruence to classify each ASV. Statistical analysis (2.6) generates final authenticated barcodes with associated quality metrics. Arrows indicate data flow between stages; colour coding denotes methodological components (specimen/molecular processing, bioinformatics, phylogenetic analysis, classification, and statistical output).

### 2.1 Specimen Collection and Processing

Beetle specimens were collected across 14 countries in four zoogeographic regions: Neotropics (Ecuador, Mexico, Honduras, Panama, French Guiana), Indo-Malaya (Malaysia, Thailand, China, India), Afrotropics (Equatorial Guinea, Mozambique, South Africa), and Palearctic (United Kingdom, Palestine). Specimens were obtained using 11 different sampling methods, comprising both passive and active approaches. Passive methods included flight interception traps, Malaise traps, pitfall traps, pan traps, Winkler extractors for leaf litter, and SLAM traps. Active collection methods comprised hand collecting, canopy fogging, sweep netting, and light-trapping. Traps were typically exposed for up to one week between collection events, using either 96% ethanol or SDS+EDTA buffer as preservation fluid. Samples were transported to the laboratory under ambient conditions and subsequently preserved at -20 degrees C or -80 degrees C for long-term storage. Detailed specimen metadata, collection protocols, and taxonomic determinations are available in the project repository (https://github.com/OunjaiS/ASV-Selection-Taxonomic-Assignment).

All Coleoptera specimens were sorted from trap catches and photographed in bulk using an AxioZoom large-scale high-resolution camera. Images were used to detect all distinct morphospecies, and a specimen image representative of each morphospecies was cropped and identified to the family. Taxonomic determinations and images were archived in a public repository (https://www.flickr.com/photos/site-100/) and identified to the family or lower level, where possible, aided by taxonomists via iNaturalist (https://www.inaturalist.org/people/site_100). DNA extractions were performed using a QIAGEN Blood and Tissue Kit (QIAGEN 2020). DNA concentrations were determined using a Qubit fluorometer (Thermo Fisher Scientific 2015), and DNA degradation was assessed via agarose gel electrophoresis.

All fieldwork and specimen collection complied with institutional and national regulations under appropriate research and export permits issued by each country's authorities (see Data Availability Statement).

### 2.2 Molecular Data Generation

A 418 bp COI fragment was amplified using tagged degenerate primers targeting a conserved region of the mitochondrial COI gene: forward (B) 5'-CCNGAYATRGCNTTYCCNCG-3' (Hajibabaei et al. 2012) and reverse (FoldR) 5'-TANACYTCNGGRTGNCCRAARAAYCA-3'. The reverse primer is a degenerate variant of jgHCO2198 (Geller et al. 2013), with inosine bases replaced by fully degenerate positions, as used in Arribas et al. (2016). This primer pair targets a conserved region of the mitochondrial COI gene and was selected for its broad amplification success across Coleoptera, as validated in prior studies on diverse arthropod assemblages (Arribas et al. 2016; Creedy et al. 2020).

Unique paired-tag combinations enabled specimen-specific identification during multiplexed sequencing (Supplementary Material 1). PCR reactions were conducted in 25 uL volumes with the following profile: 95 degrees C for 4 min, followed by 35 cycles of 95 degrees C for 30 s, 48 degrees C for 30 s, and 72 degrees C for 1 min, with a final extension at 72 degrees C for 10 min. Products were visualised by gel electrophoresis. Amplicons were pooled in equimolar ratios and subjected to library preparation before sequencing on an Illumina MiSeq platform using 2x300 paired-end chemistry across 10 independent runs conducted between 2016 and 2022.

### 2.3 Bioinformatics Pipeline

Raw sequences were demultiplexed and adapter-trimmed with Cutadapt v3.5 (Martin 2011), and paired-end reads were merged with PEAR v0.9.6 (Zhang et al. 2014). Quality filtering was performed using VSEARCH v2.21.1 (Rognes et al. 2016) with maximum expected error thresholds of 1.0 (--fastq_maxee 1). Sequences were dereplicated using 'vsearch --derep_fulllength' and denoised with 'vsearch --cluster_unoise' to eliminate sequencing errors. Translation filtering was performed using filtertranslate from the metaMATE package v1.0 (Andujar et al. 2021) with invertebrate mitochondrial genetic code (translation table 5). Subsequent filtering steps included: (1) length filtering to retain 418 +/- 3 base pair sequences, (2) translation filtering to remove sequences with stop codons or frameshifts, and (3) chimaera detection using UCHIME3 (Edgar 2016) in de novo mode. The resulting Amplicon Sequence Variants (ASVs) were each assigned a unique identifier (Supplementary Material 2). No initial taxonomic filtering was applied at this stage to avoid excluding potentially novel or divergent Coleoptera sequences. Read counts were recorded for each ASV in a single sample and across all samples and sequence runs. Taxonomic assessment and filtering of non-Coleoptera sequences were subsequently performed through phylogenetic placement analysis (Section 2.4). To benchmark our UNOISE3-based denoising and chimera removal strategy against alternative approaches, the dataset was additionally processed using the DADA2 pipeline (Callahan et al. 2016), with detailed comparative outcomes elaborated in the Results and Discussion.

### 2.4 Phylogenetic Framework Construction

A phylogenetic reference tree was constructed from 13,380 complete or partial mitochondrial genome sequences. Mitogenomes were assembled from shotgun sequencing data obtained from a subset of the specimens included in the current study (n=9,027) published elsewhere (Carpenter and Vogler 2025) (Creedy et al. 2025 (in press); Creedy et al. 2025 (in review)). These mitogenomes were supplemented with existing Coleoptera sequences downloaded from GenBank (accessed January 2023). All mitochondrial genome sequences are provided in Supplementary Material 3. Thirteen protein-coding genes were extracted using extract_genes.py (from tjcreedy/biotools, accessed March 2023; https://github.com/tjcreedy/biotools), translated using translate.py (same repository), and aligned at the amino acid level using MAFFT v7.490 (Katoh and Standley 2013) with the G-INS-i algorithm and 1,000 iterations. Alignments were back-translated to nucleotide space using backtranslate.py (same repository) and concatenated using catfasta2phyml.pl (https://github.com/nylander/catfasta2phyml). The resulting supermatrix, approximately 18,600 bp in length, was used for phylogenetic analysis under maximum likelihood using FastTree v2.1.11 (Price, Dehal, and Arkin 2010) with the GTR+CAT approximation model. The resulting tree encompassed 119 beetle families and served as a backbone constraint for phylogenetic placement of barcodes.

Three reference tree variants were constructed to assess the influence of taxon sampling density on phylogenetic placement stability (Supplementary Material 4): a 5k tree (5,000 mitogenomes selected to maximise family-level diversity), a 12k tree (12,000 mitogenomes combining taxonomic and phylogenetic diversity), and a 13.3k tree (complete dataset). All trees were inferred using identical FastTree methods, ensuring that topological differences reflected taxon sampling rather than methodological variation.

### 2.5 ASV Authentication and Classification Protocol

A hybrid classification system was developed by integrating abundance-based filtering, phylogenetic validation, and taxonomic congruence assessment to categorize all ASV records into five classes: authenticated mitochondrial haplotype, cross-contamination, intra-individual variant, environmental contamination, and technical artefact.

#### 2.5.1 Minimum Read Count Threshold Determination

Minimum read count thresholds (MRTCs) ranging from 1 to 20 reads per ASV were evaluated to identify the optimal balance between sequence reliability and specimen retention. For each threshold level, four metrics were calculated: (1) the proportion of specimens retaining at least one quality-passing ASV (specimen retention rate); (2) the percentage of total ASV records retained; (3) the percentage of unique sequence variants retained; and (4) the authentication success rate (proportion of specimens that passed the threshold and yield an authenticated ASV). The preferred MRCT maximised specimen retention whilst minimising inclusion of low-abundance artefacts. ASVs falling below the selected threshold were excluded from authentication analyses but retained for comprehensive characterisation of secondary sequence patterns.

#### 2.5.2 Phylogenetic Placement and Taxonomic Validation

ASVs meeting the MRCT were aligned to the COI reference region using MAFFT (L-INS-i algorithm) and subjected to phylogenetic tree searches with FastTree, employing the three reference tree variants (5k, 12k, and 13.3k trees; see Section 2.4) as backbone constraints. For each ASV in each tree, two metrics were calculated: (1) the phylogenetic distance to the nearest authenticated reference sequence (cophenetic distance of query and reference under the GTR+CAT model) using the ape package v5.6-2 (Paradis and Schliep 2019), and (2) the congruence between family-level placement on the tree and the morphotaxonomic identification of the sequenced specimen. A 'match' status (denoted taxonomy_match='match') was assigned when phylogenetic and morphological family assignments were concordant, serving as a prerequisite for positive ASV authentication. Conversely, a non-match status (taxonomy_match not equal to 'match') indicated non-authentic ASV status. This binary variable served as a primary criterion for subsequent ASV status classification.

For beetle families absent from the mitogenome reference database, phylogenetic placement was assessed at broader taxonomic levels (superfamily or suborder), with authentication relying more heavily on abundance criteria and comparative phylogenetic distances to taxonomically proximate families. Such cases were flagged for manual review.

#### 2.5.3 Classification Framework

ASVs were classified based on combined criteria of abundance, phylogenetic position, and co-occurrence with closely related reads (Table 1). Classification was applied at two levels: specimen-level (authenticated, failed) and ASV-level (authenticated, cross-contamination, intra-individual variants, environmental contamination, technical artefacts).

'Failed' specimens were those that yielded no ASVs after quality filtering, possibly due to PCR inhibition, low DNA concentration, DNA degradation, primer mismatch, or sequencing failure. 'Authenticated' ASVs were required to meet two criteria: (1) family-level taxonomic concordance with morphological identification (taxonomy_match='match' status) in at least one of the three constraint trees, and (2) highest read abundance within the specimen's ASV pool. When multiple ASVs met these criteria, the most abundant sequence was designated the primary authenticated haplotype, with secondary candidates each assigned to one of the following categories:

Cross-contamination was designated when an ASV showed exact sequence identity to authenticated sequences from different specimens but was present at moderate to low abundance in the focal specimen. Source tracking was performed where possible to identify the contamination pathways, including field-based transfer (specimens in the same trap), laboratory cross-contamination (specimens in the same extraction or PCR batch), and sequencing-related index hopping.

Intra-individual variants were characterised by strict taxonomic congruence (taxonomy_match='match' status in at least one constraint tree) and consistent co-occurrence with other authenticated variants. Secondary ASVs detected in >=2 independent specimens alongside the same authenticated 'anchor' ASV were classified as putative intra-individual variants, representing nuclear (NUMTs) or heteroplasmic mitochondrial variants. The consistent co-amplification in two or more individuals distinguished stable genomic features from sporadic technical artefacts or single-sample contamination.

Environmental contamination was inferred based on cross-family phylogenetic placement (taxonomy_match not equal to 'match' status across all three constraint trees), indicating ASVs from taxonomically distinct organisms such as ingested prey, parasites, commensals, or environmental DNA acquired during trapping and handling. In addition, ASVs with the same-family placement but phylogenetic distances of >=0.50 branch length units from all authenticated references across the three constraint trees were also considered environmental contaminants, representing divergent lineages unlikely to reflect intraspecific variation.

Technical artefacts included all ASVs with read counts below the MRCT threshold, encompassing sequencing errors and low-quality amplicons. ASVs above the MRCT were also classified as artefacts if they exhibited properties inconsistent with genuine biological variants, such as abnormal sequence length or very low sequence quality. Phylogenetic validation was not applied to sequences ultimately assigned to this category.

**Table 1. ASV Classification Criteria Framework**

| Criterion | Authenticated ASVs | Intra-individual Variants (NUMTs/Heteroplasmy) | Cross-contamination | Environmental Contamination | Technical Artefacts |
|---|---|---|---|---|---|
| Read abundance | Read count above MRCT; highest abundance in sample | Meets MRCT; not highest abundance | Meets MRCT; moderate to low abundance | Meets MRCT | Below MRCT threshold (< 4 reads) OR meets MRCT but has other disqualifying properties |
| Phylogenetic distance from known reference | Minimal distance to reference | Low distance to reference | Not primary criterion; Distance reflects source sample | Large distance to reference | Not applicable |
| Taxonomic placement | Matches expected family (taxonomy_match='match') | Matches expected family (taxonomy_match='match') | Irrelevant (ASV belongs to another sample's taxonomy) | Does not match expected family (taxonomy_match not equal to 'match') | Not applicable |
| Occurrence pattern | Primary ASV per specimen | Co-occurs with the same Authenticated ASV in >=2 specimens + similar sequence abundance ratio | Is an Authenticated ASV elsewhere in the dataset | Inconsistent with the sample's expected taxonomy | Exhibits properties like abnormal length or very low sequence quality |
| Defining Feature(s) | Highest read count + phylogenetic/morphotaxonomic match | Co-occurrence pattern + phylogenetic/morphotaxonomic match | Sequence identity matches another sample's Authenticated ASV | Phylogenetic/morphotaxonomic mismatch | Random experimental error, technical noise |

### 2.6 Statistical Analysis

The statistical significance of authentication determinants was evaluated using chi-square tests of independence with the Bonferroni correction for multiple comparisons. Effect sizes were quantified using Cramer's V with 95% confidence intervals, interpreted as small (V < 0.1), medium (0.1 <= V < 0.3), or large (V >= 0.3). For continuous variables, Mann-Whitney U tests and independent t-tests were used to compare distributions between authenticated and non-authenticated ASVs, with Cohen's d calculated to estimate effect size.

Mutual Information analysis was conducted to assess non-linear relationships between predictor variables and authentication success. Statistical significance was assessed at alpha = 0.05 with appropriate corrections for multiple testing. All analyses were conducted on the complete dataset of ASV records from qualified specimens (MRCT >= 4) using Python 3.11 with pandas v2.1.0, scipy v1.11.2, scikit-learn v1.3.0, BioPython v1.81, matplotlib v3.7.2, and seaborn v0.12.2.

---

## 3. Results

### 3.1 Dataset Overview and Quality Assessment

The analysis included 18,533 beetle specimens from 14 countries across four biogeographic regions, representing 119 morphologically identified families. Amplicons sequenced on Illumina MiSeq (nine runs, 2016-2022) yielded 36,459,895 reads after quality filtering, comprising 64,544 unique ASVs across 175,954 total "ASV records", i.e., the sum of ASVs in all samples, reflecting the presence of ASVs in multiple specimens or (in a few cases) in multiple libraries of resequenced single specimens.

Individual specimens exhibited variable ASV complexity, with a mean of 9.5 ASVs per specimen (median: 5, range: 0-427, SD: 15.9) before abundance filtering (Figure 2A). A total of 2,986 specimens (16.1%) presented a single ASV, whilst 2,287 specimens (12.3%) yielded two ASVs. The 25th, 50th, and 75th percentiles were 2, 5, and 11 ASVs, respectively. Read count distribution across ASV records was highly heterogeneous, with a mean of 207 reads (median: 4, range: 0-29,101, SD: 1,090) (Figure 2B). The strongly right-skewed distribution indicated that most ASV records had low abundance, with the median of 4 reads matching the selected minimum threshold. Quality metrics showed a mean sequence quality (Q score) of 35.2, with 89.3% achieving scores of Q>30. Length clustering around 418 bp was tight (mean: 418.2 bp, SD: 2.1 bp), with 97.8% within 418+/-3 bp. Reading frame integrity was maintained in 94.6% of sequences.

The dataset encompassed tropical rainforests (68.2%), temperate forests (18.4%), montane habitats (8.7%), and arid regions (4.7%), spanning elevations from sea level to 3,200 m. Collection methods included passive traps (57.2%) and active collection (42.8%), with preservation in 96% ethanol (64.3%) or SDS+EDTA buffer (35.7%). Taxonomic representation was skewed: the five most abundant families comprised 52.1% of specimens (Staphylinidae 21.4%, Chrysomelidae 10.3%, Curculionidae 9.8%, Scarabaeidae 5.8%, Carabidae 4.8%), whilst 47 families (39.2%) had <10 specimens each. Family-level morphological identification confidence based on high-resolution images was high (94.7% certain, 5.3% ambiguous).

Figure 2: ASV complexity and read count distributions. (A) ASV per specimen distribution showing the number of specimens with a given number of ASVs (on a log scale, range 1-407), and the 25th to 95th percentiles of specimens yielding the total number of ASVs. (B) The frequency of ASVs with a certain read count. The vertical line at MRCT = 4 indicates the minimum threshold; ASV records meeting or exceeding this threshold (>=4 reads) were retained for further analysis, whilst those below this threshold (<4 reads) were classified as technical artefacts.

### 3.2 Authentication Protocol and Threshold Selection

MRCTs (1-20 reads) were systematically evaluated for the kind of ASVs removed (Figure 3A). A threshold of 4 reads optimally balanced coverage and reliability: it retained 17,007 specimens (91.8%) with >=1 qualifying ASV, retained 91,402 ASV records (52.0%) representing 57,976 unique variants (89.8%), and achieved authentication success rates (i.e. matching morphotaxonomic identification and phylogenetic placement at family level) that plateaued beyond this threshold (Figure 3B). The 84,552 records below the threshold (48.0%) were classified as technical artefacts (mean: 2.2 reads, median: 2.0, mode: 1.0), with 68.4% singletons and 23.7% doubletons.

Accepted specimens exhibited diverse ASV complexity, with 1-5 ASVs obtained in 52.3% of libraries, 6-10 ASVs in 28.1%, 11-20 ASVs in 14.2%, and >20 ASVs in 5.4% of cases (mean 5.4, median 4). This reduction from the pre-filtering mean of 9.5 ASVs per specimen (Section 3.1) reflected the removal of low-abundance technical artefacts through the MRCT=4 threshold. The 1,526 failed specimens (8.2%) showed non-random distribution: failure rates varied from 4.2-18.7% across batches, 2.1-24.6% across countries, and 0-89.2% across families, indicating systematic factors determining failure rates (DNA quality, amplification efficiency, batch-specific issues). The authentication success was largely stable beyond this threshold, while the proportion of retained specimens and unique ASVs dropped sharply (Fig. 3).

Figure 3: (A) The effect of varying MRCTs on data retention, separate for all ASVs, unique ASVs, and specimens. (B) The correlation of authentication success and MRCT (ASV retention). Authentication success refers to the match of morphotaxonomic identification and phylogenetic placement of the corresponding ASV. The red circle marks the preferred MRCT=4 at which point 89.8% of ASVs are retained, as is also evident from the vertical line in Figure 3A.

### 3.3 Phylogenetic Authentication

Phylogenetic placement within the 13,380-mitogenome reference tree (Supplementary Figure S1) authenticated 14,715 specimens (86.52% of 17,007 threshold-passing specimens; 79.40% of the total 18,533 specimens). These specimens were represented by 15,901 authenticated ASV records (17.40% of threshold-passing records) representing 12,923 unique sequences. The remaining 2,292 threshold-passing specimens (13.48%) failed authentication despite meeting abundance criteria.

Phylogenetic distance distributions differed markedly between authenticated and non-authenticated ASVs. Authenticated ASVs showed minimal distances to nearest references (mean: 0.008, median: 0.000, 95th percentile: 0.028), whilst non-authenticated sequences exhibited substantially larger distances (mean: 0.954, median: 0.892, 95th percentile: 2.145). This represents a significant 119-fold difference in the mean phylogenetic distance from a reference sequence (Mann-Whitney U test: p<0.001; Cohen's d=-0.820). Authenticated ASVs also exhibited substantially elevated read counts (mean: 1,904, median: 650) compared to the initial dataset (mean: 207, median: 4), a significant, nearly tenfold increase (Mann-Whitney U test: p<0.001; Cohen's d=1.384; large effect size), confirming that they were the dominant mitochondrial haplotypes.

Among authenticated specimens, 13,538 (92.0%) yielded exactly one authenticated ASV, with the remainder presenting two or up to NN ASVs. Out of these, 11,861 of 15,901 authenticated ASVs (74.6%) matched mitochondrial genomes, with the remainder authenticated based on family-level phylogenetic context. Authentication success varied across sequencing batches (60.1-100.0%, mean: 81.2%, SD: 13.4%), with a 39.9 percentage point range strongly correlated with mean sequencing depth (Spearman's rho = 0.78, p = 0.013) but not with batch size (rho = -0.23, p = 0.551). Geographic variation (30.1-90.8% success) reflected campaign-specific preservation and logistics rather than global biogeographic patterns, as within-region variation exceeded between-region differences.

### 3.4 Authentication Failure Analysis

Of the 175,954 ASV records, 160,053 (90.96%) failed authentication and were assigned to four failure categories (Table 2; Figure 4A). Authenticated ASVs (9.04%) exhibited zero phylogenetic distance to their nearest reference and the highest read abundance. Failure categories differed in phylogenetic distance, read abundance, and taxonomic composition (Table 1; Figure 4B).

Intra-individual variants showed the lowest phylogenetic distances, with 48.2% below 0.05 branch length units. Environmental contamination displayed the highest distances, reflecting inclusion of non-coleopteran sequences. Cross-sample contamination exhibited intermediate distances, consistent with within-Coleoptera transfer. Co-occurrence analysis identified 1,790 unique intra-individual variant sequences across 4,496 specimens, with 69.4% occurring in exactly two specimens, 24.3% in three to four specimens, and 6.3% in five or more. Staphylinidae harboured the highest number and proportion of intra-individual variants.

Comparative analysis using the DADA2 pipeline confirmed that while statistical error correction effectively removes Technical Artefacts similar to our MRCT filter, it cannot resolve biological contamination. Specifically, 75,481 records representing NUMTs and Environmental Contamination passed DADA2's denoising filters, underscoring the necessity of specimen-level phylogenetic authentication (detailed comparison in Discussion 4.5; Supplementary Figure S2).

**Table 2. Summary information for four classes of ASV records (n = 175,954).**

| Category | Count | Percentage | Mean Reads | Median Phyl. Dist. | Mean Phyl. Dist. | Key Taxonomic/Structural Notes |
|---|---|---|---|---|---|---|
| Authenticated | 15,901 | 9.04 | 1,904 | 0 | 0 | 100% taxonomic match at family level |
| Technical Artefacts | 102,083 | 58.02 | 2.2 | 1.136 | 1.244 | Random distribution, no clustering |
| Environmental Contamination | 24,977 | 14.2 | 138 | 1.282 | 1.722 | 71.7% cross-family Coleoptera; 19.7% non-Coleoptera (11.4% Diptera); 2.6% bacterial |
| Cross-sample Contamination | 13,163 | 7.48 | 100 | 1.083 | 1.023 | 90.0% source-traced (78.4% lab, 15.3% field) |
| Intra-individual Variants | 19,820 | 11.27 | 60 | 0.053 | 0.143 | 48.2% <0.05 dist.; 1,790 stable across 4,496 specimens; highest in Staphylinidae |
| **Total** | **175,954** | **100** | -- | -- | -- | -- |

Figure 4: ASV Classification and Phylogenetic Distance Analysis (n = 175,954 records) (A) Pie chart illustrating the proportional distribution of all ASV records across the five classification categories. (B) Box plots comparing the phylogenetic distance (in branch length units) to the nearest authenticated reference for each ASV category, with mean and SD of distances calculated.

### 3.5 Factors Affecting Authentication Success

Authentication success varied substantially across subsets of specimens, as follows:

(1) Geographic sets (Figure 5A; Table S1) differed widely in authentication success ranging from 30.1-90.8% (mean: 76.4%, SD: 16.2%), with the highest proportion in the United Kingdom (90.8%), Malaysia (88.7%), Panama (87.1%), and the lowest Palestine (30.1%) and Honduras (52.5%). Within-region variation exceeded between-region differences, indicating campaign-specific effects (chi-squared=1,268.25, df=13, p<0.001, Cramer's V=0.262, medium effect).

(2) Collection methodology (Figure 5B; Table S2) also influenced outcome (chi-squared=563.02, df=11, p<0.001, Cramer's V=0.178, medium effect). Active collection methods achieved higher success rates (Hand: 91.2%, Winkler: 90.9%, Slam: 90.8%) compared to most passive methods, though Malaise traps (84.8%) and FIT traps (84.1%) remained effective. Pitfall traps showed the lowest success (71.2%). The preservation fluid also influenced the outcomes: ethanol achieved a success rate of 88.6%, compared to 83.9% for the SDS+EDTA buffer.

(3) Sequencing batch (Figure 5C; Table S3) effects were substantial (60.1-100.0% range, mean: 81.2%, SD: 13.4%), with significant heterogeneity (chi-squared=1,181.31, df=8, p<0.001, Cramer's V=0.252, medium effect). Batch success was strongly correlated with sequencing depth (rho = 0.78, p = 0.013) and library complexity (rho = 0.71, p = 0.032).

(4) Family membership (Figure 5D, Table S4) exhibited the most substantial variation (0-100% success range), representing the most significant effect (chi-squared=2,971.57, df=119, p<0.001, Cramer's V=0.410, large effect). Among families with >=10 specimens, 15 achieved a success rate of >=90%: Cicindelidae and Disteniidae (both 100%), followed by Artematopodidae (95.0%), Ptilodactylidae (93.3%), and Chrysomelidae (93.0%). Conversely, 18 families exhibited <70% success: Anamorphidae, Melandryidae (all 0.0%), Eucinetidae (5.9%), Ptinidae (6.6%). Small-bodied families (<3 mm) showed reduced success (52.3%) compared to medium-bodied (3-10 mm: 81.7%) and large-bodied (>10 mm: 88.9%) families. Family effects remained consistent across regions (correlation rho = 0.83-0.91, all p < 0.001), supporting a biological rather than a geographic basis.

Figure 5: Assessment of ASV authentication success for different subsets of data, separated by: (A) geographic sets (country); (B) collection methodology; (C) sequencing batch; and (D) taxonomic family.

Multivariate feature importance analysis quantified the relative contribution of each predictor of authentication success at the specimen and ASV levels (Figure 6; Table 5). At the specimen level, family identity was the strongest predictor, followed by sequencing batch and geographic origin. At the ASV level, phylogenetic distance to the nearest authenticated reference sequence was the dominant predictor (feature importance 0.510), reflecting its role as a fundamental authentication criterion, as family-level phylogenetic congruence constitutes a prerequisite for authentication. Other significant predictors included read abundance (absolute read count of the ASV; importance 0.285) and proportional representation (relative proportion of the ASV's reads within the specimen's total read pool; importance 0.142). Additional features, including total ASV count per specimen and sequence quality scores, had comparatively minor effects. Comprehensive statistical tests with effect-size quantification for all comparisons are provided in Table S5, confirming the robustness and biological significance of the observed patterns.

Figure 6: Multivariate feature importance hierarchy at specimen and ASV levels

### 3.6 NUMT Validation through Codon Usage Bias

To provide molecular evidence supporting the classification of intra-individual variants as NUMTs, we compared codon usage patterns between Authenticated ASVs (n = 8,533 unique sequences) and putative Nuclear Pseudogenes (n = 26,416). Authenticated ASVs exhibited significantly higher GC content (34.72 ± 3.53%) than NUMTs (34.26 ± 3.62%; Mann-Whitney U test, p = 3.67 × 10⁻²⁷), with the difference most pronounced at the first codon position (44.59% vs 43.86%; p = 1.54 × 10⁻³⁸). Analysis of Relative Synonymous Codon Usage (RSCU) revealed that 30 of 62 sense codons differed significantly between groups (Benjamini-Hochberg adjusted p < 0.05; Table 3). Principal Component Analysis of per-ASV RSCU values separated Authenticated and NUMT clusters along PC1, which explained 14.6% of total variance (Figure 7, Panel C). NUMTs showed higher codon bias scores (9.68 ± 2.49 vs 9.36 ± 2.31; p = 4.10 × 10⁻²²) and lower hydrophobic amino acid content (65.31 ± 2.45% vs 65.54 ± 1.65%; p = 1.13 × 10⁻⁵), consistent with relaxed functional constraints on nuclear-integrated sequences. Three representative NUMTs from different beetle families (Curculionidae, Staphylinidae, Cerambycidae) each showed convergent evidence from low intra-individual abundance (0.08–5.17% reads), elevated phylogenetic distance (1.68–2.04), and divergent codon usage relative to the authenticated ASV from the same specimen (Figure 7).

Table 3: Statistical comparison of nucleotide composition and amino acid composition between Authenticated ASVs (n = 8,533) and putative NUMTs (n = 26,416). Significance was determined using Mann-Whitney U tests, with effect size measured by rank-biserial correlation (r). RSCU comparison for all 62 sense codons is provided in Supplementary Table S7.

**Section A: Nucleotide Composition**

| Metric | Authenticated ASVs (Mean ± SD) | Putative NUMTs (Mean ± SD) | U-statistic | p-value | Effect Size (r) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| GC Content (%) | 34.72% ± 3.53% | 34.26% ± 3.62% | 121,447,638 | < 0.001 (3.67e-27) | -0.078 |
| AT Content (%) | 65.28% ± 3.53% | 65.74% ± 3.62% | 103,960,090 | < 0.001 (3.67e-27) | 0.078 |
| GC at 1st Codon Position (%) | 44.59% ± 4.00% | 43.86% ± 4.39% | 123,206,592 | < 0.001 (1.54e-38) | -0.093 |
| GC at 2nd Codon Position (%) | 43.50% ± 1.47% | 43.53% ± 1.68% | 111,573,811 | 0.159 (ns) | 0.010 |
| GC at 3rd Codon Position (%) | 16.06% ± 7.58% | 15.39% ± 7.61% | 119,456,384 | < 0.001 (7.54e-17) | -0.060 |

**Section B: Amino Acid Composition**

| Metric | Authenticated ASVs (Mean ± SD) | Putative NUMTs (Mean ± SD) | U-statistic | p-value | Effect Size (r) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Hydrophobic Amino Acids (%) | 65.54% ± 1.65% | 65.31% ± 2.45% | 116,209,982 | < 0.001 (1.13e-05) | -0.031 |
| Leucine (%) | 17.01% ± 1.08% | 17.02% ± 1.33% | 109,381,080 | < 0.001 (2.49e-05) | 0.029 |
| Internal Stop Codons | 0.0 ± 0.0 | 0.0 ± 0.0 | 112,703,864 | 1.000 (ns) | 0.000 |

Figure 7: Codon usage bias analysis comparing Authenticated ASVs (n = 8,533 unique sequences) and putative NUMTs (n = 26,416 unique sequences). (A) GC content (%) violin plot; asterisks denote statistical significance (Mann-Whitney U test, **** p < 0.001). (B) Effective Number of Codons (ENC) violin plot showing higher codon bias scores in NUMTs relative to Authenticated ASVs. (C) Principal Component Analysis (PCA) of per-ASV Relative Synonymous Codon Usage (RSCU) values; PC1 explains 14.6% of total variance and separates the two groups. (D) Top 20 codons with the largest RSCU differential between Authenticated ASVs and NUMTs; only codons with Benjamini-Hochberg FDR-adjusted p < 0.05 are shown.

---

## 4. DISCUSSION

### 4.1 Authentication workflow and ASV complexity

While discussed extensively in the context of metabarcoding (Alberdi et al. 2018; Andujar et al. 2021; Elbrecht, Peinert, and Leese 2017; Elbrecht et al. 2017; Lamb et al. 2019), the issue of extraneous sequences in deep amplicon sequencing from single specimens ("megabarcoding") has received less attention. Here, we show the challenges of implementing a filtering procedure to extract the true barcode haplotypes from a background of cross-contamination and technical artefacts that can confound specimen identification and diversity estimates, but also highlight the potential of such data for studying species interactions, pseudogene variation, and symbionts. At the core of our workflow are two authentication criteria: (i) minimum read count thresholds (MRCTs) and (ii) congruence of phylogenetic position and morphological identification. We selected a relatively low MRCT of 4 reads, which balanced high specimen retention (91.77%) with the removal of more than half of presumed artefactual ASVs. On average, each of the retained specimens yielded 9.5 ASVs, underlining the need for stringent filtering to recover the single authentic mitochondrial haplotype (Figure 2). Overall, the workflow extracted the presumed true mitochondrial haplotype for 86.52% of specimens (79.40% when PCR failures are included).

### 4.2 ASV recovery and parameters of authentication success

Our categorisation of ASVs based on relative abundance, distance to known reference sequences, and phylogenetic position revealed that technical artefacts constituted the largest failure class (102,083 records, 58.02%), dominated by low-abundance variants with variable phylogenetic distances from members of the expected family (mean: 1.244). They illustrate the challenges of high-throughput amplicon sequencing, especially when working with partially degraded field samples that increase PCR amplification artefacts, sequencing errors, and sample contamination. Although many factors unique to each of the 14 sampling campaigns hampered strict standardisation, several trends aligned with expectations. Active collection methods achieved higher success rates than most passive methods, presumably because specimens were preserved immediately rather than remaining in traps for prolonged periods, thereby reducing exposure to environmental factors that accelerate DNA hydrolysis and oxidative damage (Marquina et al. 2019). Prolonged exposure in passive traps, particularly in tropical environments with high temperatures and humidity, promotes enzymatic degradation by endogenous nucleases, accelerates depurination, and increases the likelihood of cross-contamination from co-trapped organisms through tissue dissolution and DNA leaching into the preservation fluid. The kind of preservation fluid also affects authentication success, with ethanol outperforming SDS+EDTA. Ethanol rapidly dehydrates tissues and inactivates nucleases, providing superior DNA preservation, whereas SDS+EDTA buffer, while effective at lysing cells and chelating divalent cations to inhibit DNases, may permit residual enzymatic activity at ambient tropical temperatures and promote DNA fragmentation during extended field exposure (Zizka et al. 2019). Trap type had an additional effect; e.g., Malaise traps performed slightly better than flight interception traps, although these differences were overlain by the use of ethanol vs. SDS+EDTA in either trap type. Pitfall traps had the lowest success rate (71.2%), likely due to multiple compounding factors: water ingress dilutes the preservation fluid and accelerates DNA degradation, soil-derived humic acids act as potent PCR inhibitors (Schrader et al. 2012), and the entrapment of large non-target organisms (e.g. slugs, amphibians) introduces abundant exogenous DNA that elevates background contamination levels. Additionally, the open design of pitfall traps exposes specimens to UV radiation and temperature fluctuations, both of which promote strand breakage and reduce the proportion of intact template molecules available for amplification. Finally, leaf litter sampling showed low scores, probably because of the inclusion of very small specimens, such as minute Ptiliidae that had a very low success rate, rather than any features associated with this trap type per se.

The latter highlights the strong effect of taxonomic identity as a predictor of success at the specimen level. Authentication success varied widely among families (0-100%), reflecting differences in morphology, physiology, and genome features. Families such as Chrysomelidae (93.0% success) likely yielded good DNA due to thinner cuticles or larger average body sizes, whereas heavily sclerotised or minute taxa performed poorly. The presence of metabolic compounds and defensive secretions may also affect DNA preservation and amplification. Thus, the wide range of authentication success across families suggests that extraction protocols may need to be optimised for particular morphological types. Geography also had a significant effect, capturing campaign-specific preservation differences, consistent with earlier studies (Marquina, Andersson, and Ronquist 2019). In particular, campaigns in Palestine (30.1%) and Honduras (52.5%) faced severe problems with specimen preservation, which is reflected in their low authentication success. In addition, sequencing batch effects had a surprisingly high impact (60.1-100.0%), highlighting the importance of technical consistency and standardised laboratory procedures.

From the ASV perspective, the probability that a given ASV is an authentic mitochondrial haplotype depended prominently on its phylogenetic position, as the correct family placement is a prerequisite for recognition as an authentic ASV. Thus, this parameter trivially has great predictive importance in the feature analysis, but other factors also contribute, notably read count and similarity to a sequence database entry, reflecting the fact that true mitochondrial copies are present in higher abundance (Shokralla et al. 2015). Conversely, the total ASV count per sample and the phylogenetic proximity to a mitogenome had little effect (Figure 6B). The mitogenomes used in the backbone tree were, in part, derived from the same specimens as the ASVs, but this does not appear to be relevant for the authentication pipeline, which uses the placement to family only. The current mitogenome tree comprising >13,000 mitogenomes is presumably of sufficient quality and taxon density to allow for this placement, which is likely to be accurate for any collection of Coleoptera. If used in other groups, the mitogenome tree may need to be established first, but the coverage of relevant sequences is growing rapidly (Cameron 2025), allowing the wider adoption of phylogeny-based ASV authentication.

### 4.3 Classification of secondary ASVs

Beyond the 'authentic' mitochondrial haplotype and technical artefacts, numerous secondary ASVs were recognised and classified by their phylogenetic position, read count, and recurrence in multiple libraries. The largest group (24,977 records, 14.20% of all records from accepted specimens) comprised ASVs with great phylogenetic distance to the primary ASV or morphological identification. Most of these (71.7%) represented cross-contamination with beetle taxa appearing as primary ASVs in other samples. Contamination may have arisen either in the traps or during laboratory extraction and sequencing steps, but due to the combinations in which these specimens were processed, it was not possible to determine the dominant factor. The remaining ASVs in this category originated from non-beetle arthropods (19.7%), probably due to trap-by-catch, although some may reflect ecological interactions, such as prey items or parasitoids. Bacterial Rickettsiales (2.6%) were also detected, indicating a potential endosymbiont signal.

A further 19,820 ASV records were closely related to a primary ASV (mean distance: 0.143) and likely represent NUMTs or heteroplasmic copies. For in-depth characterization of putative NUMTs, we strictly evaluated only the intra-species variants that specifically co-occurred alongside an authenticated anchor ASV within the same individual specimen. By testing for consistent co-occurrence, we validated 1,790 NUMTs across 4,496 specimens, providing strong evidence of stable nuclear insertion rather than random sequencing artefacts. Validated NUMTs were particularly frequent in Staphylinidae (345 variants), Chrysomelidae, and Scarabaeidae, which may reflect lineage-specific NUMT abundance or simply the higher numbers of sequenced specimens in these diverse families. High recovery rates of NUMTs extend previous observations (Andujar et al. 2021; Song et al. 2008), and co-occurrence with authenticated anchor ASVs provides a quantitative procedure for their systematic detection. Validated NUMTs should be recorded alongside their anchor ASVs, as they are useful not only for filtering metabarcode data (Andujar et al. 2021) but also as potential markers for population genetics (Pons and Vogler 2005). Final confirmation requires genome sequencing (Hebert et al. 2025), but consistent co-amplification with the primary mitochondrial haplotype makes NUMTs readily available as nuclear gene markers without requiring additional sequencing effort.

The 7.5% cross-contamination rate observed in our dataset, with 78.4% traced to laboratory sources and 15.3% to field transfer, underscores the need for improved protocols in future large-scale barcoding studies. Several mitigation strategies can substantially reduce cross-contamination at different stages. In the field, minimising specimen co-exposure time in traps and using individual preservation vials for high-priority specimens can reduce inter-specimen DNA transfer. In the laboratory, physical separation of pre- and post-PCR workspaces, the use of dedicated equipment, and regular decontamination with UV irradiation or bleach treatment reduce carry-over. The adoption of unique dual indexing (UDI), in which each specimen receives a unique combination of i5 and i7 indices, virtually eliminates index hopping, which has been documented at rates of 0.1-0.5% on patterned flow cells (Costello et al. 2018; MacConaill et al. 2018). The routine inclusion of negative controls (extraction blanks and PCR no-template controls) at a frequency of at least one per 96-well plate enables the detection and quantification of background contamination levels, allowing post hoc identification of affected samples. Computational approaches, including the source-tracking framework applied here, can be further refined using occupancy modelling to distinguish true species detections from contamination events (Ficetola et al. 2015). Together, these measures can substantially reduce the contamination rate and improve the proportion of specimens yielding authenticated barcodes.

The codon usage bias analysis provides additional molecular evidence supporting the NUMT classification. The significantly different RSCU patterns between Authenticated ASVs and putative NUMTs are consistent with the expected relaxation of purifying selection following nuclear integration (Bensasson et al. 2001; Hazkani-Covo et al. 2010). The higher codon bias scores in NUMTs, combined with reduced hydrophobic amino acid content, suggest that these sequences are accumulating mutations unconstrained by the functional requirements of the COI protein. However, effect sizes were modest (rank-biserial r = 0.03–0.09), which is expected given that many NUMTs represent recent insertions that retain substantial sequence similarity to their mitochondrial source. We acknowledge that definitive confirmation of nuclear localisation requires genomic approaches such as long-read genome sequencing, chromosome mapping, or PCR amplification with nuclear-specific flanking primers (Bensasson et al. 2001). Nevertheless, the convergence of multiple independent lines of evidence — low intra-individual abundance, elevated phylogenetic distance, and significantly different codon usage patterns — collectively provides strong support for the NUMT classification in this dataset.

### 4.4 Implications for biodiversity assessment and metabarcoding

Deep barcode sequencing of individual specimens remains underexplored, as most large-scale studies focusing simply on top-abundant reads (Hebert et al. 2025; Meier et al. 2016) or cluster the variation into OTUs (Hartop et al. 2024; Jabot et al. 2025; Yeo et al. 2021). Our results show that this practice overlooks the complexity of authentic, artefactual, and secondary sequences, which differ systematically across taxa, preservation conditions, sequencing batches, and others, which inevitably bias species detection. The implications of these findings for biodiversity studies are two-fold: (i) The biases in read abundances and contamination observed in single-specimen sequencing will inevitably be exacerbated in metabarcoding mixtures. When exposed to mixed templates, PCR preferentially amplifies well-preserved, high-concentration variants, obscuring other true variants and elevating certain contaminants instead. Our study helps predict which sequence types are likely to be under- or over-represented and how this affects the diversity estimates across taxonomic, morphological, and ecological groups. However, by using phylogenetic matching as a second criterion of validation, additional ASVs can be accepted if some information about the specimen mixture is available (e.g., images or reference sequences from related samples). (ii) Accurate ASV filtering is essential when using these sequences for expanding the reference sets of barcode data. The barcoding literature frequently laments the lack of adequate reference databases, especially in tropical invertebrate groups (Zinger et al. 2025). This gap can be reduced if validated high-throughput barcodes are deposited in databases, rather than being used solely for identification and not reused, as is the current practice. However, these de novo-generated reference data need to be accurate at the nucleotide level and need to be validated as authentic genomic templates. Additionally, it is crucial to maintain records of secondary reads that co-amplify with the presumed true variant but are not removed by standard filtering methods, as they can impact species richness estimates and measures of genetic variation.

We thus envision a shift in biodiversity studies toward greater reliance on high-throughput single-specimen barcoding, partly superseding metabarcoding of mixed specimen samples, especially in taxonomically poorly characterised tropical insects. Large-scale single-specimen barcoding is now eminently feasible and should feed directly into the reference databases, but this requires changes in practice. First, the current reliance on clustering of sequence variants and consensus sequences to represent the species hampers the extraction of authentic ASVs. This is partly driven by the use of ONT platforms but also where single haplotypes are available from Illumina and PacBio data, given the focus on species counts (Yeo et al. 2021). To serve as a specimen-level haplotype record, clean ASV data are needed free from read errors, nuclear pseudogenes, or species mis-assignments, generated on high-accuracy platforms. Second, due to OTU clustering, high-throughput barcodes are not typically linked to an individual specimen but are instead identified at the whole-community level via BLAST or taxonomic classifiers, without an explicit link to a particular specimen or morphological identification. This obscures sample cross-contamination, morphological mis-identification, mis-classification against reference databases, and environmental contaminants, e.g., from prey items, and does not make full use of the information contained in the primary sequence data. For barcode sequences to be reusable, they need to be recorded individually and fed into the barcode databases, with appropriate circumstantial data and, where possible, linked to a particular physical voucher or specimen image.

Treating deep-sequenced amplicons as true genetic variants also links their use to phylogenetic inference and trait reconstruction. While we used phylogenetic analysis mainly for ASV authentication, the precise placement of these amplicons allows them to be integrated into annotated phylogenies, enabling trait imputation (Keck et al. 2018). We argue that full phylogenetic likelihood tree searches provide a more resolved and accurate phylogenetic position than standard "phylogenetic placement" methods designed for metagenomics, which merely add query sequences to a fixed topology (Czech et al. 2022). This approach provides a preliminary evolutionary position for each species encountered, which, when combined with associated images and collecting data, contributes directly to a global, increasingly complete tree-of-life. We conclude that investment in reference tree construction provides greater long-term benefits for biodiversity assessment from barcodes than procedural standardisation of data filtering alone.

The choice of sequencing platform significantly affects this authentication framework's reliance on exact ASV calling. While our deep Illumina MiSeq approach yielded 58.0% technical artefacts, 92.1% were low-abundance reads (1-3 reads) effectively removed by our abundance threshold, preserving high-accuracy (>Q30) single-nucleotide resolution (Callahan, McMurdie, and Holmes 2017). Alternatively, Oxford Nanopore Technology (ONT) enables ultra-high throughput (Hebert et al. 2025), but its higher error rates necessitate consensus-based clustering rather than strict ASV-level analysis. PacBio HiFi achieves single-molecule, Sanger-equivalent accuracy (Congrains et al. 2025) but remains less practical than Illumina for processing tens of thousands of specimens due to higher costs and lower throughput. Consequently, Illumina offers cost-effective exact ASV resolution, PacBio HiFi provides premium reference-quality haplotypes, and ONT maximizes specimen throughput using OTUs. Although theoretically platform-agnostic, our framework's phylogenetic placement step relies on the nucleotide-level accuracy provided by Illumina or PacBio data and would require adaptation to handle ONT's persistent consensus errors.

### 4.5 Comparison with DADA2 Pipeline

To address whether alternative ASV inference tools would produce comparable results, we compared our VSEARCH/UNOISE-based pipeline with the DADA2 framework (Callahan et al. 2016). This comparison focused on three aspects: denoising approach, chimera detection, and the authentication gap -- i.e., problematic sequences that denoising alone cannot identify.

Both pipelines share core functionality: quality filtering, error correction (denoising), and chimera detection. DADA2 employs a parametric error model learned from the data, while our pipeline uses the UNOISE3 algorithm with an alpha parameter of 2 and a minimum read count threshold (MRCT) of 4. Analysis of our MRCT showed that at MRCT=4, all 15,921 Authenticated ASV records were retained (100% sensitivity) while all 83,709 Technical Artifacts (single-read and low-abundance errors) were removed. Lowering the threshold to MRCT=2, approximating DADA2's effective minimum, would retain an additional 27,690 Technical Artifact records without recovering any additional Authenticated sequences.

Critically, even after denoising with either approach, 75,481 records (82.6% of retained non-artifact sequences) represent biological contaminants -- including 36,214 Nuclear Pseudogene (NUMT) records and 39,267 Environmental Contamination records -- that neither DADA2 nor any standard denoising pipeline can identify. These sequences pass all quality filters because they are genuine DNA sequences; their problematic nature is only revealed through specimen-level authentication, phylogenetic placement, and taxonomic congruence assessment. The absence of chimeric sequences in the final ASV table confirms that our UCHIME3 de novo chimera detection step, applied during the bioinformatics pipeline (Section 2.3), effectively removed chimeric sequences prior to ASV table construction.

We emphasise that DADA2 and our authentication framework are complementary rather than competing approaches. DADA2 excels at statistical error correction, while our framework provides the downstream specimen-level authentication necessary for single-specimen barcoding studies. The authentication framework is pipeline-agnostic and could be applied after either DADA2 or VSEARCH/UNOISE denoising (Supplementary Figure S2; Supplementary Table S6).

---

## 5. Conclusions

By combining abundance- and phylogeny-based filtering, we provide a stringent authentication framework that addresses shortcomings of existing approaches: unrecognised contamination within the traps or the laboratory, inflated species diversity estimates, unused information about ecological interactions, and uncertainty of NUMT recognition. The systematic classification of ASVs by phylogenetic distance clarifies the biological and technical causes underlying complex amplicon variation. The framework also provides a robust approach to the identification of nuclear pseudogenes in (meta)barcoding, especially in cases of consistent co-amplification with a major variant. The authenticated dataset, placed in a phylogenetic tree and underpinned with high-resolution images and detailed collecting information, provides a valuable contribution to molecular taxonomic resources and demonstrates the feasibility of large-scale single-specimen authentication in taxonomically challenging tropical ecosystems. Not least, we are providing validated barcodes for >12,000 species of mostly tropical beetles, of which only a small proportion have been sequenced previously. The approach represents a significant step towards evolutionary context validation and improved standards for large-scale biodiversity assessment.

---

## Acknowledgements

We thank all collectors who contributed specimens to this study. SO was supported by the Institute for the Promotion of Teaching Science and Technology (IPST) through the Development and Promotion of Science and Technology Talents Project (DPST) PhD scholarship, and the SITE-100 charitable contribution to the NHM. Fieldwork was supported by the Biodiversity Initiative of the Natural History Museum, London (2013-2017). Computational analyses were performed using the Natural History Museum High Performance Computing facility. We thank the governments and institutions of Ecuador, Mexico, Honduras, Panama, French Guiana, Malaysia, Thailand, China, India, Equatorial Guinea, Mozambique, South Africa, the United Kingdom, and Palestine for providing collection and export permits. We declare no conflicts of interest.

---

## References

Alberdi, Antton, Ostaizka Aizpurua, M Thomas P Gilbert, and Kristine Bohmann. 2018. "Scrutinizing key steps for reliable metabarcoding of environmental samples." Methods in Ecology and Evolution 9 (1): 134-147.

Andujar, Carmelo, Thomas J Creedy, Paula Arribas, Heriberto Lopez, Antonia Salces-Castellano, Antonio Jose Perez-Delgado, Alfried P Vogler, and Brent C Emerson. 2021. "Validated removal of nuclear pseudogenes and sequencing artefacts from mitochondrial metabarcode data." Molecular Ecology Resources 21 (6): 1772-1787.

Arribas, Paula, Carmelo Andujar, Kevin Hopkins, Matthew Shepherd, and Alfried P Vogler. 2016. "Metabarcoding and mitochondrial metagenomics of endogean arthropods to unveil the mesofauna of the soil." Methods in Ecology and Evolution 7 (9): 1071-1081.

Bensasson, Douda, Daxi Zhang, Daniel L Hartl, and Godfrey M Hewitt. 2001. "Mitochondrial pseudogenes: evolution's misplaced witnesses." Trends in Ecology & Evolution 16 (6): 314-321.

Callahan, Benjamin J, Paul J McMurdie, Michael J Rosen, Andrew W Han, Amy Jo A Johnson, and Susan P Holmes. 2016. "DADA2: High-resolution sample inference from Illumina amplicon data." Nature Methods 13 (7): 581-583.

Callahan, Benjamin J, Paul J McMurdie, and Susan P Holmes. 2017. "Exact sequence variants should replace operational taxonomic units in marker-gene data analysis." The ISME journal 11 (12): 2639-2643.

Cameron, Stephen L. 2025. "Insect mitochondrial genomics: A decade of progress." Annual Review of Entomology 70 (1): 83-101.

Carpenter, Fiona L, and Alfried P Vogler. 2025. "If the tape were played again: convergent evolution of clade sizes and taxonomic composition in two tropical assemblages of Coleoptera." Ecography 2025 (6): e07786.

Chua, Physilia YS, Sarah J Bourlat, Cameron Ferguson, Petra Korlevic, Leia Zhao, Torbjorn Ekrem, Rudolf Meier, and Mara KN Lawniczak. 2023. "Future of DNA-based insect monitoring." Trends in Genetics 39 (7): 531-544.

Congrains, Carlos, Felix Bremer, Julian R Dupuis, Scott M Geib, Luc Leblanc, Daniel Rubinoff, and Michael San Jose. 2025. "CCS-Consensuser: A Haplotype-Aware Consensus Generator for PacBio Amplicon Sequences." Molecular Ecology Resources 25 (7): e14113.

Costello, Maura, Mark Fleharty, Justin Abreu, Yossi Farjoun, Sarah Ferriera, Lucinda Holmes, Brian Groves, Melanie Melber, Tracy Butler, and Niall Lennon. 2018. "Characterization and remediation of sample index swaps by non-redundant dual indexing on massively parallel sequencing platforms." BMC genomics 19 (1): 332.

Creedy, Thomas J, Hannah Norman, Cuong Q Tang, Kai Qing Chin, Carmelo Andujar, Paula Arribas, Rory S O'Connor, Claire Carvell, David G Notton, and Alfried P Vogler. 2020. "A validated workflow for rapid taxonomic assignment and monitoring of a national fauna of bees (Apiformes) using high throughput DNA barcoding." Molecular ecology resources 20 (1): 40-53.

Czech, Lucas, Pierre Barbera, and Alexandros Stamatakis. 2019. "Methods for automatic reference trees and multilevel phylogenetic placement." Bioinformatics 35 (7): 1151-1158.

Czech, Lucas, Alexandros Stamatakis, Micah Dunthorn, and Pierre Barbera. 2022. "Metagenomic analysis using phylogenetic placement--a review of the first decade." Frontiers in Bioinformatics 2: 871393.

Edgar, Robert C. 2016. "UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing." BioRxiv: 081257.

Elbrecht, Vasco, Bianca Peinert, and Florian Leese. 2017. "Sorting things out: Assessing effects of unequal specimen biomass on DNA metabarcoding." Ecology and evolution 7 (17): 6918-6926.

Elbrecht, Vasco, Ecaterina Edith Vamos, Kristian Meissner, Jukka Aroviita, and Florian Leese. 2017. "Assessing strengths and weaknesses of DNA metabarcoding-based macroinvertebrate identification for routine stream monitoring." Methods in Ecology and Evolution 8 (10): 1265-1275.

Ficetola, Gentile Francesco, Jerome Pansu, Aurelie Bonin, Eric Coissac, Clement Giguet-Covex, Marta De Barba, Lucie Gielly, Charline M Lopes, Frederic Boyer, Francois Pompanon, Gilles Raye, and Pierre Taberlet. 2015. "Replication levels, false presences and the estimation of the presence/absence from eDNA metabarcoding data." Molecular Ecology Resources 15 (3): 543-556.

Geller, Jonathan, Christopher Meyer, Michael Parker, and Robert Hawk. 2013. "Redesign of PCR primers for mitochondrial cytochrome c oxidase subunit I for marine invertebrates and application in all-taxa biotic surveys." Molecular Ecology Resources 13 (5): 851-861.

Hartop, Emily, Leshon Lee, Amrita Srivathsan, Mirkka Jones, Pablo Pena-Aguilera, Otso Ovaskainen, Tomas Roslin, and Rudolf Meier. 2024. "Resolving biology's dark matter: species richness, spatiotemporal distribution, and community composition of a dark taxon." BMC biology 22 (1): 215.

Hazkani-Covo, Einat, Raymond M Zeller, and William Martin. 2010. "Molecular poltergeists: mitochondrial DNA copies (numts) in sequenced nuclear genomes." PLoS Genetics 6 (2): e1000834.

Hajibabaei, Mehrdad, Michelle A Spall, Mehrnaz Shokralla, and Steven Van Rossum. 2012. "Environmental barcoding: a next-generation sequencing approach for biomonitoring applications using river benthos." PloS one 7 (4): e36298.

He, Yeyan, Siqin Ge, and Hongbin Liang. 2025. "A Genome-Wide Analysis of Nuclear Mitochondrial DNA Sequences (NUMTs) in Chrysomelidae Species (Coleoptera)." Insects 16 (2): 150.

Hebert, Paul DN, Dan G Bock, and Sean WJ Prosser. 2023. "Interrogating 1000 insect genomes for NUMTs: A risk assessment for estimates of species richness." PLoS One 18 (6): e0286620.

Hebert, Paul DN, Alina Cywinska, Shelley L Ball, and Jeremy R DeWaard. 2003. "Biological identifications through DNA barcodes." Proceedings of the Royal Society of London. Series B: Biological Sciences 270 (1512): 313-321.

Hebert, Paul DN, Robin Floyd, Saeideh Jafarpour, and Sean WJ Prosser. 2025. "Barcode 100k specimens: in a single nanopore run." Molecular Ecology Resources 25 (1): e14028.

Jabot, Franck, Gwenaelle Auger, Pauline Bonnal, Mathilde Pizaine, Marilyn Roncoroni, Sandrine Revaillot, and Julien Pottier. 2025. "Use of massive DNA barcoding to monitor biodiversity: A test on forest soil macrofauna." Forest Ecology and Management 595: 123004.

Katoh, Kazutaka, and Daron M Standley. 2013. "MAFFT multiple sequence alignment software version 7: improvements in performance and usability." Molecular biology and evolution 30 (4): 772-780.

Keck, Francois, Valentin Vasselon, Frederic Rimet, Agnes Bouchez, and Maria Kahlert. 2018. "Boosting DNA metabarcoding for biomonitoring with phylogenetic estimation of operational taxonomic units' ecological profiles." Molecular Ecology Resources 18 (6): 1299-1309.

Lamb, Philip D, Ewan Hunter, John K Pinnegar, Simon Creer, Richard G Davies, and Martin I Taylor. 2019. "How quantitative is metabarcoding: A meta-analytical approach." Molecular ecology 28 (2): 420-430.

Liu, Shanlin, Chentao Yang, Chengran Zhou, and Xin Zhou. 2017. "Filling reference gaps via assembling DNA barcodes using high-throughput sequencing--moving toward barcoding the world." GigaScience 6 (12): gix104.

MacConaill, Laura E, Robert T Burns, Anwesha Nag, Haley A Coleman, Michael K Sleber, Liam Grimmett, Matthew Testa, Christopher J Chatelle, Julia L Rizzolo, and Chandrani Mondal. 2018. "Unique, dual-indexed sequencing adapters with UMIs effectively eliminate index cross-talk and significantly improve sensitivity of massively parallel sequencing." BMC genomics 19 (1): 30.

Marquina, Daniel, Anders F Andersson, and Fredrik Ronquist. 2019. "New mitochondrial primers for metabarcoding of insects, designed and evaluated using in silico methods." Molecular Ecology Resources 19 (1): 90-104.

Martin, Marcel. 2011. "Cutadapt removes adapter sequences from high-throughput sequencing reads." EMBnet. journal 17 (1): 10-12.

Meier, Rudolf, Winghing Wong, Amrita Srivathsan, and Maosheng Foo. 2016. "$1 DNA barcodes for reconstructing complex phenomes and finding rare species in specimen-rich samples." Cladistics 32 (1): 100-110.

Noguerales, Victor, Emmanouil Meramveliotakis, Adrian Castro-Insua, Carmelo Andujar, Paula Arribas, Thomas J Creedy, Isaac Overcast, Helene Morlon, Brent C Emerson, and Alfried P Vogler. 2023. "Community metabarcoding reveals the relative role of environmental filtering and spatial processes in metacommunity dynamics of soil microarthropods across a mosaic of montane forests." Molecular Ecology 32 (23): 6110-6128.

Paradis, Emmanuel, and Klaus Schliep. 2019. "ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R." Bioinformatics 35 (3): 526-528.

Pons, Joan, and Alfried P Vogler. 2005. "Complex pattern of coalescence and fast evolution of a mitochondrial rRNA pseudogene in a recent radiation of tiger beetles." Molecular Biology and Evolution 22 (4): 991-1000.

Price, Morgan N, Paramvir S Dehal, and Adam P Arkin. 2010. "FastTree 2--approximately maximum-likelihood trees for large alignments." PloS one 5 (3): e9490.

QIAGEN. 2020. "DNeasy blood and tissue handbook." DNeasy Blood and Tissue Handbook.

Rognes, Torbjorn, Tomas Flouri, Ben Nichols, Christopher Quince, and Frederic Mahe. 2016. "VSEARCH: a versatile open source tool for metagenomics." PeerJ 4: e2584.

Schrader, Christin, Andreas Schielke, Linfa Ellerbroek, and Ralf Johne. 2012. "PCR inhibitors -- occurrence, properties and removal." Journal of Applied Microbiology 113 (5): 1014-1026.

Shokralla, Shadi, Teresita M Porter, Joel F Gibson, Rafal Dobosz, Daniel H Janzen, Winnie Hallwachs, G Brian Golding, and Mehrdad Hajibabaei. 2015. "Massively parallel multiplex DNA sequencing for specimen identification using an Illumina MiSeq platform." Scientific reports 5 (1): 9687.

Song, Hojun, Jennifer E Buhay, Michael F Whiting, and Keith A Crandall. 2008. "Many species in one: DNA barcoding overestimates the number of species when nuclear mitochondrial pseudogenes are coamplified." Proceedings of the national academy of sciences 105 (36): 13486-13491.

Souto-Vilaros, Daniel, Eduardo Navarro-Valencia, Ana Cecilia Zamora, Yahir Campusano, Greg PA Lamarre, Filonila Perez, Yacksecari Lopez, Ricardo Bobadilla, Jose Alejandro Ramirez Silva, and Milan Janda. 2025. "Navigating the seven seas of arthropod collection protocols: Metabarcoding arthropod diversity in a tropical forest." Methods in Ecology and Evolution 16 (10): 2395-2407.

Srivathsan, Amrita, Vivian Feng, Daniel Suarez, Brent Emerson, and Rudolf Meier. 2024. "ONTbarcoder 2.0: rapid species discovery and identification with real-time barcoding facilitated by Oxford Nanopore R10. 4." Cladistics 40 (2): 192-203.

Srivathsan, Amrita, Emily Hartop, Jayanthi Puniamoorthy, Wan Ting Lee, Sujatha Narayanan Kutty, Olavi Kurina, and Rudolf Meier. 2019. "Rapid, large-scale species discovery in hyperdiverse taxa using 1D MinION sequencing." BMC biology 17: 1-20.

Srivathsan, Amrita, Leshon Lee, Kazutaka Katoh, Emily Hartop, Sujatha Narayanan Kutty, Johnathan Wong, Darren Yeo, and Rudolf Meier. 2021. "ONTbarcoder and MinION barcodes aid biodiversity discovery and identification by everyone, for everyone." BMC biology 19: 1-21.

Taberlet, Pierre, Eric Coissac, Francois Pompanon, Christian Brochmann, and Eske Willerslev. 2012. "Towards next-generation biodiversity assessment using DNA metabarcoding." Molecular ecology 21 (8): 2045-2050.

Thermo Fisher Scientific. 2015. "User Guide: Qubit dsDNA HS Assay Kits." https://tools.thermofisher.com/content/sfs/manuals/Qubit_dsDNA_HS_Assay_UG.pdf.

Yeo, Darren, Amrita Srivathsan, Jayanthi Puniamoorthy, Foo Maosheng, Patrick Grootaert, Lena Chan, Benoit Guenard, Claas Damken, Rodzay A Wahab, and Ang Yuchen. 2021. "Mangroves are an overlooked hotspot of insect diversity despite low plant diversity." BMC biology 19 (1): 202.

Zhang, Jiajie, Kassian Kobert, Tomas Flouri, and Alexandros Stamatakis. 2014. "PEAR: a fast and accurate Illumina Paired-End reAd mergeR." Bioinformatics 30 (5): 614-620.

Zinger, Lucie, Anne-Sophie Benoiston, Yves Cuenot, Celine Leroy, Eliane Louisanna, Lucie Moreau, Frederic Petitclerc, Finn Piatscheck, Jerome Orivel, and Cecile Richard-Hansen. 2025. "Elusive tropical forest canopy diversity revealed through environmental DNA contained in rainwater." Science Advances 11 (33): eadx4909.

Zizka, Vera MA, Florian Leese, Baldwin Quadros, Alfred Gilbert, and Alexander Austin. 2019. "Assessing the influence of sample tagging and library preparation on DNA metabarcoding." Molecular Ecology Resources 19 (4): 893-899.

---

## Data Accessibility Statement

**Sequence Data:** Raw sequence reads, authenticated ASV sequences, specimen metadata, and phylogenetic data are deposited in Dryad: https://doi.org/10.5061/dryad.v41ns1s9q

The dataset includes quality-filtered ASV sequences (FASTA), ASV abundance tables, authentication classifications, specimen metadata (including collection localities, dates, trap types, preservation methods, and morphological identifications), phylogenetic trees (Newick format), COI alignments, and complete analysis outputs.

**Code Availability:** The authentication pipeline and all analysis scripts for sequence processing, ASV authentication, phylogenetic placement, and statistical analyses are available on GitHub: https://github.com/OunjaiS/ASV-Selection-Taxonomic-Assignment and archived with a permanent DOI on Zenodo: https://doi.org/10.5281/zenodo.XXXXXXX. [AUTHOR ACTION REQUIRED: Create a Zenodo release of the GitHub repository and replace XXXXXXX with the actual Zenodo DOI.] All third-party scripts used in this study, including those from the tjcreedy/biotools repository (https://github.com/tjcreedy/biotools, accessed March 2023) and catfasta2phyml (https://github.com/nylander/catfasta2phyml), are documented with version information in Supplementary Material 2 and 4.

**Specimen Images and Taxonomic Data:** High-resolution specimen images are archived on Flickr: https://www.flickr.com/photos/site-100/

Taxonomic determinations were made to family or lower taxonomic level where possible, aided by expert taxonomists via iNaturalist: https://www.inaturalist.org/people/site_100

---

## Benefit-Sharing Statement

The benefits of this research are shared through several mechanisms, in accordance with the principles of the Convention on Biological Diversity (CBD). First, through the establishment of international scientific collaboration, as reflected in the co-authorship of this study. Second, all genetic data, analysis scripts, and high-resolution specimen images are made publicly available (see Data Accessibility Statement), ensuring results are shared with the broader scientific community, including researchers and institutions in the countries of origin. Finally, the research addresses a priority concern by developing a scalable framework for the biodiversity assessment of under-characterised tropical fauna, which is essential for conservation efforts and sustainable resource management.

---

## Author Contributions

SO, APV designed research; SO, HL, ZZ, MPC, TJC, CA, PA conducted fieldwork and specimen processing; SO, HL performed laboratory work; SO, HL, MPC, TJC, ZZ performed bioinformatics analyses; SO analysed data; SO and APV wrote the paper with input from all authors. All authors reviewed and approved the final manuscript.

---

## Permits

All specimens were collected in accordance with international Access and Benefit Sharing regulations. Details are available at the NHM Data Portal.
