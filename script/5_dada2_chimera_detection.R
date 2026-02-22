library(dada2)
library(Biostrings)

# Read ASV FASTA
seqs_file <- "/Users/sarawut/Desktop/Manuscript_ASV_selection/raw_data/ASV_sequences_64544.fasta"
seqs <- readDNAStringSet(seqs_file)
cat("Total ASVs loaded:", length(seqs), "\n")

# Get unique sequences with abundance as named integer vector
seq_chars <- as.character(seqs)
seq_names <- names(seqs)

# Remove duplicate sequences
uniq_seqs <- unique(seq_chars)
cat("Unique sequences:", length(uniq_seqs), "\n")

# Create named integer vector (abundance = 1 for each unique sequence)
abund <- rep(1L, length(uniq_seqs))
names(abund) <- uniq_seqs

# Run chimera detection
cat("Running DADA2 isBimeraDenovo...\n")
is_chimera <- isBimeraDenovo(abund, verbose = TRUE)

n_chimera <- sum(is_chimera)
n_nonchim <- sum(!is_chimera)
cat("\n=== DADA2 Chimera Detection Results ===\n")
cat("Input unique ASVs:", length(abund), "\n")
cat("Non-chimeric ASVs:", n_nonchim, "\n")
cat("Chimeric ASVs detected:", n_chimera, "\n")
cat("Chimera rate:", round(n_chimera / length(abund) * 100, 2), "%\n")

# Map back to original ASV IDs
chimera_seqs <- uniq_seqs[is_chimera]
chimera_mask <- seq_chars %in% chimera_seqs
chimera_ids <- seq_names[chimera_mask]
nonchim_ids <- seq_names[!chimera_mask]

cat("Chimera ASV IDs:", length(chimera_ids), "\n")
cat("Non-chimera ASV IDs:", length(nonchim_ids), "\n")

# Save
output_dir <- "/Users/sarawut/Desktop/Manuscript_ASV_selection/data_analysis/dada2_benchmarking"
write.csv(data.frame(asv_id = chimera_ids, chimera = TRUE),
          file.path(output_dir, "dada2_chimeras.csv"), row.names = FALSE)
write.csv(data.frame(asv_id = nonchim_ids, chimera = FALSE),
          file.path(output_dir, "dada2_nonchimeras.csv"), row.names = FALSE)
cat("Results saved.\n")
