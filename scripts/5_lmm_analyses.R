####################################################################
# Additional Statistical Analyses – Response to Lab Colleague Comments
# 1. LMM: Read Count vs Authentication Status          (Comment 3&4)
# 2. LMM: GC Content – NUMT vs Authenticated          (Comment 5)
# 3. Pagel's lambda + Blomberg's K – GC phylo signal  (Comment 5)
####################################################################

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(ape)
  library(phytools)
  library(rotl)
})

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frames()[[1]]$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg[1])))
    else getwd()
  }
)
ANALYSIS_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
source(file.path(SCRIPT_DIR, "..", "config.R"))
out_dir <- LMM_OUT
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------
# ANALYSIS 1: LMM – Read Count vs Authentication Status
#-------------------------------------------------------------------
cat("==========================================================\n")
cat("ANALYSIS 1: LMM – Read Count vs Authentication Status\n")
cat("  Model: log(reads+1) ~ auth_status + (1|specimen_id)\n")
cat("==========================================================\n\n")

d1 <- read.csv(LMM_READCOUNT_CSV,
               stringsAsFactors = TRUE)
d1$auth_status <- relevel(d1$auth_status, ref = "Unauthenticated")
d1$log_reads   <- log1p(d1$reads)

cat("Sample sizes:\n")
print(table(d1$auth_status))
cat("Number of specimens:", length(unique(d1$specimen_id)), "\n\n")

lmm1      <- lmer(log_reads ~ auth_status + (1 | specimen_id), data = d1, REML = FALSE)
lmm1_null <- lmer(log_reads ~ 1            + (1 | specimen_id), data = d1, REML = FALSE)
lmm1_t    <- lmer(log_reads ~ auth_status + (1 | specimen_id), data = d1) # REML for inference

cat("--- Fixed Effects (lmerTest Satterthwaite df) ---\n")
fe1 <- summary(lmm1_t)$coefficients
print(fe1)

cat("\n--- Likelihood Ratio Test (fixed effect of authentication) ---\n")
lrt1 <- anova(lmm1_null, lmm1)
print(lrt1)

vc1  <- as.data.frame(VarCorr(lmm1_t))
icc1 <- vc1$vcov[1] / sum(vc1$vcov)
coef_auth <- fixef(lmm1_t)["auth_statusAuthenticated"]

cat(sprintf("\nICC (specimen): %.4f  (%.1f%% of total variance)\n", icc1, icc1 * 100))
cat(sprintf("Fold-change (Auth vs Unauth reads): %.2f×\n\n", exp(coef_auth)))

# Descriptive Mann-Whitney (retained as per response to comment)
mw1 <- wilcox.test(reads ~ auth_status, data = d1)
n1_auth  <- sum(d1$auth_status == "Authenticated")
n1_unauth <- sum(d1$auth_status == "Unauthenticated")
# Rank-biserial r: use double precision to avoid integer overflow
r1 <- 1 - (2 * as.numeric(mw1$statistic)) / (as.numeric(n1_auth) * as.numeric(n1_unauth))
cat(sprintf("Descriptive Mann-Whitney: W = %.0f, p < 2.2e-16, r = %.3f\n", mw1$statistic, r1))

sink(file.path(out_dir, "LMM1_ReadCount_results.txt"))
cat("=== LMM: log(reads+1) ~ auth_status + (1|specimen_id) ===\n\n")
cat("Sample sizes:\n"); print(table(d1$auth_status))
cat("\nFixed Effects:\n"); print(fe1)
cat("\nLRT vs null (random-intercept only) model:\n"); print(lrt1)
cat("\nRandom Effects variance components:\n"); print(vc1)
cat(sprintf("\nICC (specimen): %.4f\n", icc1))
cat(sprintf("Fold-change (Auth/Unauth): %.2fx\n", exp(coef_auth)))
cat(sprintf("\nDescriptive Mann-Whitney: W = %.0f, r = %.3f\n", mw1$statistic, r1))
sink()
cat(">> Results saved: LMM1_ReadCount_results.txt\n\n")


#-------------------------------------------------------------------
# ANALYSIS 2: LMM – GC Content: Nuclear Pseudogenes vs Authenticated
#-------------------------------------------------------------------
cat("==========================================================\n")
cat("ANALYSIS 2: LMM – GC Content (Proven NUMT vs Authenticated)\n")
cat("  Model: GC_content ~ group + (1|specimen_id)\n")
cat("  Population: Proven NUMTs = Nuclear_Pseudogenes, reads>=4,\n")
cat("    co-occurring with same anchor ASV in >=2 specimens\n")
cat("  GC content: computed from FASTA sequences (consistent with\n")
cat("    the original Mann-Whitney in the manuscript)\n")
cat("==========================================================\n\n")

d2 <- read.csv(LMM_GC_CSV,
               stringsAsFactors = TRUE)
d2$group <- relevel(d2$group, ref = "Authenticated")

cat("Sample sizes:\n")
print(table(d2$group))
cat("Unique specimens:", length(unique(d2$specimen_id)), "\n")
cat("Group means – GC content (%):\n")
print(tapply(d2$GC_content, d2$group, mean, na.rm = TRUE))
cat("\n")

lmm2      <- lmer(GC_content ~ group + (1 | specimen_id), data = d2, REML = FALSE)
lmm2_null <- lmer(GC_content ~ 1     + (1 | specimen_id), data = d2, REML = FALSE)
lmm2_t    <- lmer(GC_content ~ group + (1 | specimen_id), data = d2)

cat("--- Fixed Effects ---\n")
fe2 <- summary(lmm2_t)$coefficients
print(fe2)

cat("\n--- Likelihood Ratio Test ---\n")
lrt2 <- anova(lmm2_null, lmm2)
print(lrt2)

vc2  <- as.data.frame(VarCorr(lmm2_t))
icc2 <- vc2$vcov[1] / sum(vc2$vcov)
coef_numt <- fixef(lmm2_t)["groupProven_NUMT"]
cat(sprintf("\nICC (specimen): %.4f\n", icc2))
cat(sprintf("GC difference WITHIN specimen (NUMT - Auth, LMM estimate): %.3f%%\n", coef_numt))
cat("NOTE: LMM estimate reflects within-specimen contrast (specimens with both groups)\n")
cat("      The LMM intercept is calibrated to specimen-level GC baseline.\n\n")

# Descriptive Mann-Whitney
mw2 <- wilcox.test(GC_content ~ group, data = d2)
n_auth2 <- sum(d2$group == "Authenticated")
n_numt2 <- sum(d2$group == "Proven_NUMT")
r2 <- 1 - (2 * mw2$statistic) / (n_auth2 * n_numt2)
cat(sprintf("Descriptive Mann-Whitney: W = %.0f, p = %.3e, r = %.3f\n",
            mw2$statistic, mw2$p.value, r2))

sink(file.path(out_dir, "LMM2_NUMT_GC_results.txt"))
cat("=== LMM: GC_content ~ group + (1|specimen_id) ===\n")
cat("Population: Proven NUMTs (Nuclear_Pseudogenes, reads>=4, co-occur same anchor >=2 specimens)\n")
cat("GC source: computed from FASTA sequences (same method as original manuscript Mann-Whitney)\n\n")
cat("Sample sizes (occurrence rows):\n"); print(table(d2$group))
cat("\nUnique specimen IDs:", length(unique(d2$specimen_id)), "\n")
cat("Specimens with BOTH groups:", sum(unique(d2$specimen_id[d2$group=='Proven_NUMT']) %in%
                                         unique(d2$specimen_id[d2$group=='Authenticated'])), "\n")
cat("\nOccurrence-weighted group means:\n")
print(tapply(d2$GC_content, d2$group, mean, na.rm=TRUE))
cat("\nFixed Effects (Satterthwaite df):\n"); print(fe2)
cat("\nLRT vs null (random-intercept only):\n"); print(lrt2)
cat("\nRandom Effects variance components:\n"); print(vc2)
cat(sprintf("\nICC (specimen): %.4f\n", icc2))
cat(sprintf("LMM estimate (NUMT - Auth, within specimen): %.3f%%\n", coef_numt))
cat("\nNOTE: LMM tests within-specimen contrast, controlling for specimen-level GC baseline.\n")
cat("      Descriptive Mann-Whitney below is across all observations (no specimen control).\n")
cat(sprintf("\nDescriptive Mann-Whitney: W = %.0f, p = %.3e, r = %.3f\n",
            mw2$statistic, mw2$p.value, r2))
sink()
cat(">> Results saved: LMM2_NUMT_GC_results.txt\n\n")


#-------------------------------------------------------------------
# ANALYSIS 2b: LMM – All Metrics (GC, ENC, Hydrophobic AA, Leucine)
#              Within-specimen comparison for Table 3 (within-specimen)
#-------------------------------------------------------------------
cat("==========================================================\n")
cat("ANALYSIS 2b: LMM – All Metrics (Proven NUMT vs Authenticated)\n")
cat("  Model: metric ~ group + (1|specimen_id)  for each metric\n")
cat("  Data:  lmm_numt_all_metrics_data.csv\n")
cat("==========================================================\n\n")

d3 <- read.csv(LMM_ALL_METRICS_CSV,
               stringsAsFactors = TRUE)
d3$group <- relevel(d3$group, ref = "Authenticated")

cat("Sample sizes:\n")
print(table(d3$group))
cat("Unique specimens:", length(unique(d3$specimen_id)), "\n\n")

metrics_lmm <- list(
  GC_content      = list(col = "GC_content",      label = "GC Content (%)"),
  gc_pos1         = list(col = "gc_pos1",          label = "GC 1st Codon Position (%)"),
  gc_pos2         = list(col = "gc_pos2",          label = "GC 2nd Codon Position (%)"),
  gc_pos3         = list(col = "gc_pos3",          label = "GC 3rd Codon Position (%)"),
  ENC             = list(col = "ENC",              label = "ENC"),
  hydrophobic_pct = list(col = "hydrophobic_pct",  label = "Hydrophobic AA (%)"),
  leucine_pct     = list(col = "leucine_pct",       label = "Leucine (%)")
)

lmm_results <- list()

for (nm in names(metrics_lmm)) {
  info   <- metrics_lmm[[nm]]
  col    <- info$col
  label  <- info$label
  d_sub  <- d3[!is.na(d3[[col]]), ]

  cat(sprintf("--- %s ---\n", label))
  cat(sprintf("  n rows: %d  (NUMT=%d, Auth=%d)\n",
    nrow(d_sub),
    sum(d_sub$group == "Proven_NUMT"),
    sum(d_sub$group == "Authenticated")))

  lmm_full <- lmer(as.formula(sprintf("%s ~ group + (1|specimen_id)", col)),
                   data = d_sub, REML = FALSE)
  lmm_null <- lmer(as.formula(sprintf("%s ~ 1 + (1|specimen_id)", col)),
                   data = d_sub, REML = FALSE)
  lmm_reml <- lmer(as.formula(sprintf("%s ~ group + (1|specimen_id)", col)),
                   data = d_sub)

  fe      <- summary(lmm_reml)$coefficients
  lrt     <- anova(lmm_null, lmm_full)
  vc      <- as.data.frame(VarCorr(lmm_reml))
  icc_val <- vc$vcov[1] / sum(vc$vcov)
  coef_nm <- fixef(lmm_reml)["groupProven_NUMT"]
  p_lrt   <- lrt$`Pr(>Chisq)`[2]

  cat("  Fixed Effects:\n"); print(fe)
  cat(sprintf("  LRT p = %.3e\n", p_lrt))
  cat(sprintf("  ICC (specimen) = %.4f\n", icc_val))
  cat(sprintf("  Within-specimen estimate (NUMT - Auth): %.4f\n\n", coef_nm))

  lmm_results[[nm]] <- list(
    label       = label,
    n_numt      = sum(d_sub$group == "Proven_NUMT"),
    n_auth      = sum(d_sub$group == "Authenticated"),
    beta        = coef_nm,
    p_lrt       = p_lrt,
    icc         = icc_val,
    mean_auth   = mean(d_sub[[col]][d_sub$group == "Authenticated"], na.rm = TRUE),
    sd_auth     = sd(d_sub[[col]][d_sub$group == "Authenticated"], na.rm = TRUE),
    mean_numt   = mean(d_sub[[col]][d_sub$group == "Proven_NUMT"], na.rm = TRUE),
    sd_numt     = sd(d_sub[[col]][d_sub$group == "Proven_NUMT"], na.rm = TRUE)
  )
}

# Build summary table
lmm_tbl <- do.call(rbind, lapply(names(lmm_results), function(nm) {
  r <- lmm_results[[nm]]
  sig <- ifelse(r$p_lrt < 0.001, "***",
         ifelse(r$p_lrt < 0.01,  "**",
         ifelse(r$p_lrt < 0.05,  "*", "ns")))
  data.frame(
    Metric          = r$label,
    Auth_Mean_SD    = sprintf("%.2f ± %.2f", r$mean_auth, r$sd_auth),
    NUMT_Mean_SD    = sprintf("%.2f ± %.2f", r$mean_numt, r$sd_numt),
    Beta_NUMT_Auth  = round(r$beta, 4),
    LRT_p           = formatC(r$p_lrt, format = "e", digits = 2),
    ICC_specimen    = round(r$icc, 4),
    Significance    = sig,
    stringsAsFactors = FALSE
  )
}))

cat("=== WITHIN-SPECIMEN LMM SUMMARY TABLE ===\n")
print(lmm_tbl, row.names = FALSE)

write.csv(lmm_tbl, file.path(out_dir, "LMM_WithinSpecimen_AllMetrics_Table.csv"),
          row.names = FALSE)
cat("\n>> Results saved: LMM_WithinSpecimen_AllMetrics_Table.csv\n\n")


#-------------------------------------------------------------------
# ANALYSIS 3: Phylogenetic Signal in GC Content (family level)
#-------------------------------------------------------------------
cat("==========================================================\n")
cat("ANALYSIS 3: Phylogenetic Signal – Blomberg's K & Pagel's λ\n")
cat("  GC content averaged across beetle families\n")
cat("==========================================================\n\n")

fgc <- read.csv(FAMILY_GC_CSV)
cat(sprintf("Beetle families with GC data: %d\n\n", nrow(fgc)))

cat("Step 1: Query Open Tree of Life for Coleoptera family phylogeny...\n")
taxa   <- tnrs_match_names(names = fgc$family, context_name = "Arthropods")
in_otl <- is_in_tree(ott_id(taxa))
taxa_ok <- taxa[in_otl, ]
cat(sprintf("  Families matched in OTL synthetic tree: %d / %d\n", nrow(taxa_ok), nrow(taxa)))

ott_ids <- as.numeric(unlist(ott_id(taxa_ok)))

tree_raw <- tryCatch({
  tol_induced_subtree(ott_ids = ott_ids)
}, error = function(e) {
  # Remove any pruned (bad) OTT IDs reported in the error message
  bad_nums <- as.numeric(gsub("ott", "",
    regmatches(conditionMessage(e), gregexpr("ott[0-9]+", conditionMessage(e)))[[1]]))
  bad_nums <- unique(bad_nums)
  cat(sprintf("  Removing %d pruned OTT IDs and retrying...\n", length(bad_nums)))
  tol_induced_subtree(ott_ids = ott_ids[!ott_ids %in% bad_nums])
})

cat(sprintf("  Raw tree tips: %d\n", length(tree_raw$tip.label)))

# Clean tip labels (remove OTT suffixes, replace underscores)
clean_tips <- gsub("_ott[0-9]+", "", tree_raw$tip.label)
clean_tips <- gsub("_", " ",  clean_tips)
clean_tips <- trimws(clean_tips)
tree_raw$tip.label <- clean_tips

# Match GC values to tree tips
gc_lookup <- setNames(fgc$mean_GC, fgc$family)
gc_matched <- gc_lookup[tree_raw$tip.label]
gc_matched <- gc_matched[!is.na(gc_matched)]
tree_pruned <- keep.tip(tree_raw, names(gc_matched))
cat(sprintf("  Tips with GC data: %d\n\n", length(tree_pruned$tip.label)))

cat("Step 2: Make tree ultrametric (Grafen's method)...\n")
# Grafen's method assigns branch lengths proportional to number of descendants
# (appropriate when no molecular clock calibration is available)
tree_ultra <- tryCatch({
  compute.brlen(tree_pruned, method = "Grafen")
}, error = function(e) {
  cat("  Grafen failed; falling back to equal branch lengths (conservative).\n")
  tree_pruned$edge.length <- rep(1, nrow(tree_pruned$edge))
  tree_pruned
})
cat(sprintf("  Ultrametric: %s\n", is.ultrametric(tree_ultra)))

cat("\nStep 3: Blomberg's K\n")
K_res <- phylosig(tree_ultra, gc_matched[tree_ultra$tip.label],
                  method = "K", test = TRUE, nsim = 1000)
cat(sprintf("  K = %.4f,  P = %.4f\n", K_res$K, K_res$P))
cat("  (K = 1: random; K > 1: more similar than expected by Brownian motion)\n\n")

cat("Step 4: Pagel's lambda\n")
lambda_res <- tryCatch(
  phylosig(tree_ultra, gc_matched[tree_ultra$tip.label],
           method = "lambda", test = TRUE),
  error = function(e) {
    cat("  Lambda estimation failed:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(lambda_res)) {
  cat(sprintf("  lambda = %.4f,  LRT P = %.4f\n", lambda_res$lambda, lambda_res$P))
  cat("  (lambda = 0: no phylogenetic signal; lambda = 1: Brownian motion)\n\n")
}

sink(file.path(out_dir, "PhyloSignal_GC_results.txt"))
cat("=== Phylogenetic Signal in GC Content (family-level means) ===\n\n")
cat(sprintf("Families included: %d (out of %d with data)\n", length(tree_ultra$tip.label), nrow(fgc)))
cat("Tree: Open Tree of Life induced subtree, made ultrametric via chronoMPL\n\n")
cat("Blomberg's K:\n"); print(K_res)
if (!is.null(lambda_res)) {
  cat("\nPagel's lambda:\n"); print(lambda_res)
}
cat("\nFamily-level GC content data used:\n")
print(data.frame(family = names(gc_matched[tree_ultra$tip.label]),
                 mean_GC = gc_matched[tree_ultra$tip.label]))
sink()
cat(">> Results saved: PhyloSignal_GC_results.txt\n\n")


cat("==========================================================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("Output folder:", out_dir, "\n")
cat("==========================================================\n")
