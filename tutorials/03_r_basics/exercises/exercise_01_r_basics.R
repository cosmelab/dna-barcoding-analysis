#!/usr/bin/env Rscript
################################################################################
# Exercise 01: R Basics for DNA Barcoding
#
# Practice fundamental R skills using DNA sequence data
# Complete each section and check your answers with the solutions file
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

# EXERCISE 1: WORKING WITH VECTORS ===========================================

# 1.1 Create a vector of sequence lengths for 10 mosquito samples
# Use the sample() function to generate random lengths between 500 and 700
# Hint: sequence_lengths <- sample(500:700, 10, replace=TRUE)

# YOUR CODE HERE:


# 1.2 Calculate the following statistics for your sequence lengths:
# - Mean length
# - Median length
# - Standard deviation
# - Minimum and maximum lengths

# YOUR CODE HERE:


# 1.3 How many sequences are longer than 600 bp?
# Hint: Use sum() with a logical condition

# YOUR CODE HERE:


# 1.4 What percentage of sequences are longer than 600 bp?

# YOUR CODE HERE:


# EXERCISE 2: DATA FRAMES ===================================================

# 2.1 Create a data frame with information about 8 mosquito samples
# Include these columns:
# - sample_id: MOS001 through MOS008
# - genus: Mix of Aedes, Culex, Anopheles (your choice)
# - location: Different cities in California
# - seq_length: Random lengths between 600-700
# - quality_score: Random quality scores between 20-40

# YOUR CODE HERE:


# 2.2 Display the first 3 rows of your data frame

# YOUR CODE HERE:


# 2.3 What is the average sequence length by genus?
# Hint: Use aggregate() or tapply()

# YOUR CODE HERE:


# 2.4 Create a new column called "qc_status" that is "PASS" if
# quality_score > 30 AND seq_length > 620, otherwise "FAIL"

# YOUR CODE HERE:


# EXERCISE 3: FILTERING DATA ================================================

# 3.1 Filter your data frame to show only Aedes samples

# YOUR CODE HERE:


# 3.2 Filter to show samples that passed QC (qc_status == "PASS")

# YOUR CODE HERE:


# 3.3 Filter to show Culex samples from a specific location

# YOUR CODE HERE:


# 3.4 How many samples of each genus passed QC?
# Hint: Use table() on filtered data

# YOUR CODE HERE:


# EXERCISE 4: FUNCTIONS =====================================================

# 4.1 Write a function called "gc_content" that takes a DNA sequence
# and returns the GC percentage
# Test sequence: "ATCGATCGGCGCATATGCGC"

gc_content <- function(sequence) {
  # YOUR CODE HERE

}

# Test your function:
# gc_content("ATCGATCGGCGCATATGCGC")  # Should return 55%


# 4.2 Write a function called "reverse_complement" that takes a DNA sequence
# and returns its reverse complement (A<->T, G<->C, reversed)
# Hint: Use chartr() for substitution and strsplit/rev/paste for reversal

reverse_complement <- function(sequence) {
  # YOUR CODE HERE

}

# Test your function:
# reverse_complement("ATCG")  # Should return "CGAT"


# 4.3 Write a function called "assess_quality" that takes seq_length and
# quality_score and returns:
# - "Excellent" if length > 650 and quality > 35
# - "Good" if length > 600 and quality > 30
# - "Poor" otherwise

assess_quality <- function(seq_length, quality_score) {
  # YOUR CODE HERE

}

# Test your function with different values


# EXERCISE 5: LOOPS ========================================================

# 5.1 Use a for loop to print "Processing sample X" for each sample_id
# in your data frame

# YOUR CODE HERE:


# 5.2 Use a for loop to calculate the GC content for these sequences:
sequences <- c("ATCGATCGATCG", "GCGCGCGCGCGC", "ATATATATATATAT")
# Store results in a vector

# YOUR CODE HERE:


# 5.3 Use a while loop to simulate sequencing until you get a
# high-quality sequence (quality > 35)
# Hint: Use sample() to generate random quality scores

# YOUR CODE HERE:


# EXERCISE 6: DATA ANALYSIS CHALLENGE =======================================

# 6.1 Load (or create) a data frame with 20 mosquito samples including:
# - Sample ID
# - Genus (Aedes, Culex, Anopheles)
# - Species
# - Location
# - Sequence length
# - Quality score
# - Year collected (2020-2023)

# YOUR CODE HERE:


# 6.2 Answer these questions about your data:
# a) What is the overall pass rate (samples with quality > 30)?


# b) Which genus has the highest average sequence length?


# c) Has sequence quality improved over the years?


# d) Which location has the most samples?


# 6.3 Create a summary report showing:
# - Total samples
# - Samples by genus
# - Pass/fail counts
# - Average quality by genus

# YOUR CODE HERE:


# EXERCISE 7: STRING MANIPULATION ===========================================

# 7.1 You have these sample IDs with different formats:
# Clean them to a standard format (e.g., "MOS_001")
messy_ids <- c("mos001", "MOS-002", "Mos 003", "MOS004")

# YOUR CODE HERE:


# 7.2 Extract genus and species from these labels:
labels <- c("Aedes_aegypti_CA_2023", "Culex_pipiens_NY_2023", "Anopheles_gambiae_KE_2022")
# Create two vectors: genera and species

# YOUR CODE HERE:


# 7.3 Create publication-ready labels in the format "Genus species (Location)"
# Example: "Aedes aegypti (California)"

# YOUR CODE HERE:


# EXERCISE 8: APPLY FUNCTIONS ==============================================

# 8.1 Use sapply() to calculate the length of each sequence:
sequences <- list("ATCG", "GCGCGC", "ATATAT", "GCGCGCGCGC")

# YOUR CODE HERE:


# 8.2 Use lapply() to convert all sequences to lowercase

# YOUR CODE HERE:


# 8.3 Use apply() to calculate mean quality per sample (rows)
# and mean quality per position (columns)
quality_matrix <- matrix(c(30, 35, 40, 25, 38, 42, 32, 36, 40, 28, 35, 38),
                        nrow = 3, ncol = 4)

# YOUR CODE HERE:


# BONUS CHALLENGE ==========================================================

# Write a complete analysis function that takes a data frame of samples
# and produces a comprehensive QC report including:
# - Number of samples
# - Pass/fail summary
# - Statistics by genus
# - Quality distribution plot
# - Recommendations for re-sequencing

comprehensive_qc_report <- function(samples_df) {
  # YOUR CODE HERE

}

# Test your function on your sample data

################################################################################
# REFLECTION QUESTIONS
#
# After completing these exercises, answer these questions:
#
# 1. What's the difference between a vector and a list in R?
#
# 2. When would you use a for loop vs sapply()?
#
# 3. Why are functions useful in data analysis?
#
# 4. How would you handle missing data (NA values) in your analyses?
#
# 5. What strategies would you use to ensure your data is clean and consistent?
#
################################################################################

# Check your answers with:
# source("../solutions/solution_01_r_basics.R")
