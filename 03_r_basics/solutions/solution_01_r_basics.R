#!/usr/bin/env Rscript
################################################################################
# Solution 01: R Basics for DNA Barcoding
#
# Complete solutions with explanations
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

cat("=== EXERCISE 1: WORKING WITH VECTORS ===\n\n")

# 1.1 Create vector of sequence lengths
set.seed(123)  # For reproducibility
sequence_lengths <- sample(500:700, 10, replace = TRUE)
cat("Sequence lengths:", sequence_lengths, "\n\n")

# 1.2 Calculate statistics
cat("Mean length:", mean(sequence_lengths), "\n")
cat("Median length:", median(sequence_lengths), "\n")
cat("Standard deviation:", sd(sequence_lengths), "\n")
cat("Min length:", min(sequence_lengths), "\n")
cat("Max length:", max(sequence_lengths), "\n\n")

# 1.3 Count sequences > 600 bp
n_long <- sum(sequence_lengths > 600)
cat("Sequences > 600 bp:", n_long, "\n\n")

# 1.4 Percentage > 600 bp
percent_long <- (sum(sequence_lengths > 600) / length(sequence_lengths)) * 100
cat("Percentage > 600 bp:", round(percent_long, 1), "%\n\n")

cat("=== EXERCISE 2: DATA FRAMES ===\n\n")

# 2.1 Create data frame
set.seed(123)
samples <- data.frame(
  sample_id = paste0("MOS", sprintf("%03d", 1:8)),
  genus = c("Aedes", "Aedes", "Culex", "Anopheles", "Aedes", "Culex", "Culex", "Anopheles"),
  location = c("Riverside", "San Diego", "Los Angeles", "Riverside",
               "Irvine", "San Francisco", "Sacramento", "Fresno"),
  seq_length = sample(600:700, 8, replace = TRUE),
  quality_score = sample(20:40, 8, replace = TRUE)
)

# 2.2 First 3 rows
cat("First 3 rows:\n")
print(head(samples, 3))
cat("\n")

# 2.3 Average length by genus
cat("Average length by genus:\n")
avg_by_genus <- aggregate(seq_length ~ genus, data = samples, FUN = mean)
print(avg_by_genus)
cat("\n")

# 2.4 Add QC status
samples$qc_status <- ifelse(samples$quality_score > 30 & samples$seq_length > 620,
                            "PASS", "FAIL")
cat("Data frame with QC status:\n")
print(samples)
cat("\n")

cat("=== EXERCISE 3: FILTERING DATA ===\n\n")

# 3.1 Only Aedes samples
aedes_samples <- samples[samples$genus == "Aedes", ]
cat("Aedes samples:\n")
print(aedes_samples)
cat("\n")

# 3.2 Samples that passed QC
passed_samples <- samples[samples$qc_status == "PASS", ]
cat("Passed QC:\n")
print(passed_samples)
cat("\n")

# 3.3 Culex from specific location
culex_riverside <- samples[samples$genus == "Culex" & samples$location == "Riverside", ]
cat("Culex from Riverside:\n")
print(culex_riverside)
cat("\n")

# 3.4 Pass rate by genus
pass_by_genus <- table(samples$genus, samples$qc_status)
cat("Pass/Fail by genus:\n")
print(pass_by_genus)
cat("\n")

cat("=== EXERCISE 4: FUNCTIONS ===\n\n")

# 4.1 GC content function
gc_content <- function(sequence) {
  sequence <- toupper(sequence)
  bases <- unlist(strsplit(sequence, ""))
  g_count <- sum(bases == "G")
  c_count <- sum(bases == "C")
  total_length <- length(bases)
  gc_percent <- ((g_count + c_count) / total_length) * 100
  return(gc_percent)
}

test_seq <- "ATCGATCGGCGCATATGCGC"
cat("GC content of", test_seq, ":", round(gc_content(test_seq), 1), "%\n\n")

# 4.2 Reverse complement function
reverse_complement <- function(sequence) {
  # Convert to uppercase
  sequence <- toupper(sequence)
  # Complement (swap A<->T, G<->C)
  complement <- chartr("ATGC", "TACG", sequence)
  # Reverse
  bases <- unlist(strsplit(complement, ""))
  reversed <- paste(rev(bases), collapse = "")
  return(reversed)
}

cat("Reverse complement of ATCG:", reverse_complement("ATCG"), "\n\n")

# 4.3 Assess quality function
assess_quality <- function(seq_length, quality_score) {
  if (seq_length > 650 & quality_score > 35) {
    return("Excellent")
  } else if (seq_length > 600 & quality_score > 30) {
    return("Good")
  } else {
    return("Poor")
  }
}

cat("Quality assessments:\n")
cat("(658, 38):", assess_quality(658, 38), "\n")
cat("(620, 32):", assess_quality(620, 32), "\n")
cat("(580, 25):", assess_quality(580, 25), "\n\n")

cat("=== EXERCISE 5: LOOPS ===\n\n")

# 5.1 Print processing messages
cat("Processing messages:\n")
for (id in samples$sample_id) {
  cat("Processing sample", id, "\n")
}
cat("\n")

# 5.2 Calculate GC for multiple sequences
sequences <- c("ATCGATCGATCG", "GCGCGCGCGCGC", "ATATATATATATAT")
gc_results <- c()
for (seq in sequences) {
  gc_results <- c(gc_results, gc_content(seq))
}
cat("GC contents:", round(gc_results, 1), "\n\n")

# 5.3 Simulate sequencing until high quality
cat("Simulating sequencing:\n")
attempt <- 1
quality <- 0
set.seed(123)
while (quality <= 35) {
  quality <- sample(20:45, 1)
  cat("Attempt", attempt, ": quality =", quality, "\n")
  attempt <- attempt + 1
}
cat("Success! High quality achieved.\n\n")

cat("=== EXERCISE 6: DATA ANALYSIS CHALLENGE ===\n\n")

# 6.1 Create larger dataset
set.seed(123)
large_samples <- data.frame(
  sample_id = paste0("SEQ", sprintf("%03d", 1:20)),
  genus = sample(c("Aedes", "Culex", "Anopheles"), 20, replace = TRUE),
  species = c("aegypti", "albopictus", "pipiens", "quinquefasciatus", "gambiae", "stephensi")[
    sample(1:6, 20, replace = TRUE)],
  location = sample(c("Riverside", "San Diego", "Los Angeles", "Irvine"), 20, replace = TRUE),
  seq_length = sample(580:680, 20, replace = TRUE),
  quality_score = sample(25:40, 20, replace = TRUE),
  year = sample(2020:2023, 20, replace = TRUE)
)

# 6.2 Answer questions

# a) Overall pass rate
pass_rate <- (sum(large_samples$quality_score > 30) / nrow(large_samples)) * 100
cat("a) Overall pass rate:", round(pass_rate, 1), "%\n\n")

# b) Genus with highest average length
avg_length_genus <- aggregate(seq_length ~ genus, data = large_samples, FUN = mean)
best_genus <- avg_length_genus[which.max(avg_length_genus$seq_length), ]
cat("b) Genus with highest avg length:", best_genus$genus,
    "(", round(best_genus$seq_length, 1), "bp )\n\n")

# c) Quality over years
avg_quality_year <- aggregate(quality_score ~ year, data = large_samples, FUN = mean)
cat("c) Average quality by year:\n")
print(avg_quality_year)
cat("\n")

# d) Location with most samples
location_counts <- table(large_samples$location)
cat("d) Samples by location:\n")
print(location_counts)
cat("\n")

# 6.3 Summary report
cat("=== SUMMARY REPORT ===\n")
cat("Total samples:", nrow(large_samples), "\n")
cat("\nSamples by genus:\n")
print(table(large_samples$genus))
cat("\nPass/Fail (quality > 30):\n")
print(table(ifelse(large_samples$quality_score > 30, "PASS", "FAIL")))
cat("\nAverage quality by genus:\n")
print(aggregate(quality_score ~ genus, data = large_samples, FUN = mean))
cat("\n")

cat("=== EXERCISE 7: STRING MANIPULATION ===\n\n")

# 7.1 Clean IDs
messy_ids <- c("mos001", "MOS-002", "Mos 003", "MOS004")
clean_ids <- gsub("[-\\s]", "_", toupper(messy_ids))
# Ensure format MOS_XXX
clean_ids <- gsub("MOS_?(\\d+)", "MOS_\\1", clean_ids)
cat("Cleaned IDs:", clean_ids, "\n\n")

# 7.2 Extract genus and species
labels <- c("Aedes_aegypti_CA_2023", "Culex_pipiens_NY_2023", "Anopheles_gambiae_KE_2022")
genera <- sapply(strsplit(labels, "_"), function(x) x[1])
species <- sapply(strsplit(labels, "_"), function(x) x[2])
cat("Genera:", genera, "\n")
cat("Species:", species, "\n\n")

# 7.3 Publication-ready labels
pub_labels <- paste(genera, species, "(Location)")
cat("Publication labels:", pub_labels, "\n\n")

cat("=== EXERCISE 8: APPLY FUNCTIONS ===\n\n")

# 8.1 sapply for length
sequences <- list("ATCG", "GCGCGC", "ATATAT", "GCGCGCGCGC")
seq_lengths <- sapply(sequences, nchar)
cat("Sequence lengths:", seq_lengths, "\n\n")

# 8.2 lapply for lowercase
lower_seqs <- lapply(sequences, tolower)
cat("Lowercase sequences:\n")
print(lower_seqs)
cat("\n")

# 8.3 apply for matrix
quality_matrix <- matrix(c(30, 35, 40, 25, 38, 42, 32, 36, 40, 28, 35, 38),
                        nrow = 3, ncol = 4)
row_means <- apply(quality_matrix, 1, mean)
col_means <- apply(quality_matrix, 2, mean)
cat("Row means (samples):", round(row_means, 1), "\n")
cat("Column means (positions):", round(col_means, 1), "\n\n")

cat("=== SOLUTIONS COMPLETE ===\n")
cat("\nKey Takeaways:\n")
cat("1. Vectors are the foundation - everything builds on them\n")
cat("2. Data frames are like spreadsheets - rows and columns\n")
cat("3. Functions encapsulate reusable logic\n")
cat("4. Loops repeat operations, but apply() is often more elegant\n")
cat("5. String manipulation is essential for bioinformatics\n")
cat("6. Always check your data - use summary(), str(), head()\n")
