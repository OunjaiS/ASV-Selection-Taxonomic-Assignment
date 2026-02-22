
library(dada2)
library(Biostrings)

# Read ASV FASTA
seqs_file <- "/Users/sarawut/Desktop/Manuscript_ASV_selection/raw_data/ASV_sequences_64544.fasta"
seqs <- readDNAStringSet(seqs_file)
cat("Total ASVs loaded:", length(seqs), "\n")

# Read abundance data
abund_df <- read.csv("/Users/sarawut/Desktop/Manuscript_ASV_selection/data_analysis/dada2_benchmarking/asv_abundances_for_chimera.csv")
cat("Abundance data loaded:", nrow(abund_df), "rows\n")

# Get sequences as character vector
seq_chars <- as.character(seqs)
seq_names <- names(seqs)

# Match abundances to sequences
abund_vec <- integer(length(seq_names))
names(abund_vec) <- seq_names
matched <- match(seq_names, abund_df$asv_id)
abund_vec[!is.na(matched)] <- abund_df$total_reads[matched[!is.na(matched)]]
# Set unmatched to 1 (minimum)
abund_vec[abund_vec == 0] <- 1L

cat("Abundance range:", min(abund_vec), "-", max(abund_vec), "\n")
cat("Median abundance:", median(abund_vec), "\n")

# Create a data.frame with $sequence and $abundance for isBimeraDenovo
# Remove duplicates (keep highest abundance per unique sequence)
seq_abund_df <- data.frame(sequence = seq_chars, abundance = abund_vec, stringsAsFactors = FALSE)
seq_abund_agg <- aggregate(abundance ~ sequence, data = seq_abund_df, FUN = max)
cat("Unique sequences:", nrow(seq_abund_agg), "\n")

# Create named integer vector (the format DADA2 expects)
uniq_abund <- as.integer(seq_abund_agg$abundance)
names(uniq_abund) <- seq_abund_agg$sequence

cat("Running DADA2 isBimeraDenovo with real abundances...\n")
cat("Class:", class(uniq_abund), "\n")
is_chimera <- isBimeraDenovo(uniq_abund, verbose = TRUE)

n_chimera <- sum(is_chimera)
n_nonchim <- sum(!is_chimera)
cat("\n=== DADA2 Chimera Detection Results ===\n")
cat("Input unique ASVs:", length(uniq_abund), "\n")
cat("Non-chimeric ASVs:", n_nonchim, "\n")
cat("Chimeric ASVs detected:", n_chimera, "\n")
cat("Chimera rate:", round(n_chimera / length(uniq_abund) * 100, 2), "%\n")

# Map chimera status back to all ASV IDs
chimera_seqs <- seq_abund_agg$sequence[is_chimera]
chimera_mask <- seq_chars %in% chimera_seqs
chimera_ids <- seq_names[chimera_mask]
nonchim_ids <- seq_names[!chimera_mask]

if (length(chimera_ids) > 0) {
    write.csv(data.frame(asv_id = chimera_ids, chimera = TRUE),
              "/Users/sarawut/Desktop/Manuscript_ASV_selection/data_analysis/dada2_benchmarking//dada2_chimeras.csv", row.names = FALSE)
} else {
    write.csv(data.frame(asv_id = character(0), chimera = logical(0)),
              "/Users/sarawut/Desktop/Manuscript_ASV_selection/data_analysis/dada2_benchmarking//dada2_chimeras.csv", row.names = FALSE)
}
write.csv(data.frame(asv_id = nonchim_ids, chimera = FALSE),
          "/Users/sarawut/Desktop/Manuscript_ASV_selection/data_analysis/dada2_benchmarking//dada2_nonchimeras.csv", row.names = FALSE)

cat("Results saved.\n")
