#!/usr/bin/env Rscript
################################################################################
# R Basics for DNA Barcoding Analysis
# Script 01: R Syntax and Data Structures
#
# Learning Objectives:
# - Understand R syntax and basic operations
# - Work with vectors, data frames, and lists
# - Write functions and loops
# - Apply these skills to DNA sequence analysis
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

# SECTION 1: R BASICS AND ARITHMETIC =========================================

# R can be used as a calculator
2 + 2                    # Addition
10 - 3                   # Subtraction
4 * 5                    # Multiplication
20 / 4                   # Division
2^10                     # Exponentiation (2 to the power of 10)
10 %% 3                  # Modulo (remainder: 10 divided by 3)

# Assignment operators (create variables)
sequence_length <- 658   # This is the standard COI barcode length
sequence_length          # Print the value

# You can also use = for assignment, but <- is preferred in R
quality_threshold = 20   # Phred quality score threshold

# SECTION 2: VECTORS ========================================================

# A vector is a sequence of elements of the same type
# Creating numeric vectors
sequence_lengths <- c(658, 640, 655, 650, 658)  # c() combines values
quality_scores <- c(35, 38, 42, 30, 25, 20)

# Access elements by position (R uses 1-based indexing!)
sequence_lengths[1]      # First element
sequence_lengths[3]      # Third element
sequence_lengths[1:3]    # Elements 1 through 3

# Vector operations (element-wise)
sequence_lengths * 3     # Multiply each length by 3
sequence_lengths + 100   # Add 100 to each length
mean(sequence_lengths)   # Average length
median(sequence_lengths) # Median length
sd(sequence_lengths)     # Standard deviation
summary(sequence_lengths)# Summary statistics

# Character vectors (for DNA sequences, species names, etc.)
mosquito_genera <- c("Aedes", "Culex", "Anopheles", "Psorophora")
mosquito_genera[2]       # Access second genus

# Logical vectors (TRUE/FALSE)
high_quality <- quality_scores > 30
high_quality             # Which scores are above 30?

# Using logical vectors for filtering
quality_scores[high_quality]  # Only get scores above 30

# SECTION 3: DATA FRAMES ====================================================

# Data frames are like spreadsheets - rows and columns of different types
# This is how we'll store sequence metadata!

# Create a data frame with mosquito sample information
samples <- data.frame(
  sample_id = c("MOS001", "MOS002", "MOS003", "MOS004", "MOS005"),
  genus = c("Aedes", "Aedes", "Culex", "Anopheles", "Aedes"),
  species = c("aegypti", "albopictus", "pipiens", "gambiae", "aegypti"),
  location = c("Riverside", "San Diego", "Los Angeles", "Riverside", "Irvine"),
  seq_length = c(658, 640, 655, 650, 658),
  quality = c(35, 38, 42, 30, 35)
)

# View the data frame
print(samples)
View(samples)            # Opens in a viewer (in RStudio)

# Access columns by name (two methods)
samples$genus            # Using $
samples[["genus"]]       # Using [[]]

# Access specific rows and columns
samples[1, ]             # First row, all columns
samples[, "genus"]       # All rows, genus column
samples[1:3, c("sample_id", "genus", "species")]  # Specific rows & columns

# Get data frame dimensions
nrow(samples)            # Number of rows (samples)
ncol(samples)            # Number of columns (variables)
dim(samples)             # Both dimensions
colnames(samples)        # Column names

# SECTION 4: FILTERING AND SUBSETTING ========================================

# Filter rows based on conditions
# Get only Aedes samples
aedes_samples <- samples[samples$genus == "Aedes", ]
print(aedes_samples)

# Get samples with high quality scores
high_qual_samples <- samples[samples$quality > 35, ]
print(high_qual_samples)

# Multiple conditions (AND: &, OR: |)
aedes_riverside <- samples[samples$genus == "Aedes" &
                           samples$location == "Riverside", ]
print(aedes_riverside)

# Get samples from Riverside OR San Diego
sd_riv_samples <- samples[samples$location == "Riverside" |
                          samples$location == "San Diego", ]
print(sd_riv_samples)

# SECTION 5: FUNCTIONS ======================================================

# Functions take inputs and return outputs
# Built-in functions we've already used: mean(), c(), print()

# Create our own function: calculate GC content
# GC content is important for primer design and sequence quality
calculate_gc_content <- function(sequence) {
  # Count G and C bases
  g_count <- sum(unlist(strsplit(sequence, "")) == "G")
  c_count <- sum(unlist(strsplit(sequence, "")) == "C")
  total_length <- nchar(sequence)

  # Calculate percentage
  gc_percent <- ((g_count + c_count) / total_length) * 100

  return(gc_percent)
}

# Test the function
test_sequence <- "ATCGATCGATCGATCGGGCCCATATA"
calculate_gc_content(test_sequence)

# Function with multiple parameters
assess_sequence_quality <- function(length, quality, min_length = 500, min_quality = 30) {
  # Default values: min_length=500, min_quality=30

  if (length >= min_length & quality >= min_quality) {
    return("PASS")
  } else if (length < min_length) {
    return("FAIL: Too short")
  } else {
    return("FAIL: Low quality")
  }
}

# Test the function
assess_sequence_quality(658, 35)         # Should pass
assess_sequence_quality(450, 35)         # Too short
assess_sequence_quality(658, 25)         # Low quality
assess_sequence_quality(658, 25, min_quality = 20)  # Override default

# SECTION 6: LOOPS ==========================================================

# For loops: repeat code for each element
# Print each genus
for (genus in mosquito_genera) {
  print(paste("Genus:", genus))
}

# Loop through a sequence of numbers
for (i in 1:5) {
  print(paste("Sample", i, "has length", sequence_lengths[i]))
}

# Apply quality assessment to all samples
for (i in 1:nrow(samples)) {
  result <- assess_sequence_quality(samples$seq_length[i],
                                    samples$quality[i])
  print(paste(samples$sample_id[i], ":", result))
}

# While loops: repeat while a condition is true
count <- 1
while (count <= 5) {
  print(paste("Count is", count))
  count <- count + 1
}

# SECTION 7: LISTS ==========================================================

# Lists can contain different types of objects (unlike vectors)
# This is useful for storing complex analysis results

analysis_results <- list(
  samples = samples,
  average_length = mean(samples$seq_length),
  average_quality = mean(samples$quality),
  genera_count = table(samples$genus),
  high_quality_samples = samples[samples$quality > 35, ]
)

# Access list elements
analysis_results$average_length
analysis_results[[2]]            # Same as above (by position)
analysis_results$genera_count

# SECTION 8: CONDITIONAL STATEMENTS ==========================================

# If-else statements for decision making
sample_quality <- 42

if (sample_quality >= 40) {
  print("Excellent quality!")
} else if (sample_quality >= 30) {
  print("Good quality")
} else if (sample_quality >= 20) {
  print("Acceptable quality")
} else {
  print("Poor quality - consider re-sequencing")
}

# Vectorized ifelse (works on vectors)
samples$quality_category <- ifelse(samples$quality >= 35,
                                   "High",
                                   "Medium")
print(samples)

# SECTION 9: WORKING WITH STRINGS (DNA SEQUENCES) ===========================

# Strings are essential for working with DNA sequences
dna_sequence <- "ATCGATCGATCG"

# String length
nchar(dna_sequence)      # 12 bases

# Convert to lowercase/uppercase
tolower(dna_sequence)
toupper(dna_sequence)

# Substring (extract part of string)
substr(dna_sequence, 1, 3)    # First 3 bases (ATG)
substr(dna_sequence, 4, 6)    # Bases 4-6 (GAT)

# String concatenation
genus <- "Aedes"
species <- "aegypti"
full_name <- paste(genus, species)  # "Aedes aegypti"
full_name_no_space <- paste0(genus, "_", species)  # "Aedes_aegypti"

# Split string into characters
strsplit(dna_sequence, "")[[1]]

# SECTION 10: APPLY FAMILY ==================================================

# More elegant than loops for applying functions to data

# lapply: apply function to list, return list
sequence_list <- list("ATCG", "GCTA", "TACG")
lapply(sequence_list, nchar)  # Get length of each sequence

# sapply: apply function to list, return vector (simpler output)
sapply(sequence_list, nchar)

# apply: apply function to margins of array/matrix
# Create matrix of quality scores (samples x positions)
quality_matrix <- matrix(c(30, 35, 40, 25, 38, 42), nrow = 2, ncol = 3)
apply(quality_matrix, 1, mean)  # Mean of each row (sample)
apply(quality_matrix, 2, mean)  # Mean of each column (position)

# SECTION 11: PRACTICAL EXAMPLE - BATCH SEQUENCE ANALYSIS ===================

# Simulate reading multiple sequence files
# In real analysis, you'd read from FASTA files

sequences_data <- data.frame(
  id = paste0("SEQ", 1:10),
  length = c(658, 645, 660, 580, 670, 655, 640, 590, 665, 658),
  quality = c(38, 35, 40, 25, 42, 37, 33, 28, 39, 36),
  gc_content = c(48.2, 47.8, 49.1, 46.5, 48.9, 47.2, 48.5, 45.8, 49.3, 48.0)
)

# Analyze which sequences pass quality control
sequences_data$qc_status <- apply(sequences_data, 1, function(row) {
  assess_sequence_quality(as.numeric(row["length"]),
                         as.numeric(row["quality"]),
                         min_length = 600,
                         min_quality = 30)
})

print(sequences_data)

# Calculate summary statistics
cat("\n=== SEQUENCE ANALYSIS SUMMARY ===\n")
cat("Total sequences:", nrow(sequences_data), "\n")
cat("Passed QC:", sum(sequences_data$qc_status == "PASS"), "\n")
cat("Failed QC:", sum(sequences_data$qc_status != "PASS"), "\n")
cat("Average length:", round(mean(sequences_data$length), 1), "bp\n")
cat("Average quality:", round(mean(sequences_data$quality), 1), "\n")
cat("Average GC%:", round(mean(sequences_data$gc_content), 1), "%\n")

# SECTION 12: READING AND WRITING FILES =====================================

# In real analysis, you'll read data from files
# Common formats: CSV, TSV, FASTA

# Write data frame to CSV
write.csv(samples, "mosquito_samples.csv", row.names = FALSE)

# Read CSV file
# samples_from_file <- read.csv("mosquito_samples.csv")

# Read table with custom separator
# data <- read.table("file.tsv", sep="\t", header=TRUE)

# SECTION 13: GETTING HELP ==================================================

# R has excellent built-in documentation
?mean                    # Help for mean function
help(data.frame)         # Same as ?data.frame
??phylogenetic           # Search for topics containing "phylogenetic"
example(plot)            # Show examples for plot function

# EXERCISES TO TRY ==========================================================
#
# 1. Create a vector of 10 random sequence lengths between 500-700
#    Hint: sample(500:700, 10, replace=TRUE)
#
# 2. Calculate what percentage of your samples have quality > 30
#
# 3. Write a function that converts a DNA sequence to its complement
#    (A<->T, G<->C)
#
# 4. Create a data frame with 20 hypothetical samples including:
#    - Sample ID
#    - Genus (choose from: Aedes, Culex, Anopheles)
#    - Collection date
#    - Latitude/Longitude
#
# 5. Filter your data frame to get only samples from a specific genus
#    and calculate the mean quality score for that genus
#
################################################################################
# END OF SCRIPT
#
# Summary:
# - R uses <- for assignment
# - Vectors hold elements of the same type
# - Data frames are like spreadsheets
# - Functions encapsulate reusable code
# - Loops and apply() for repetitive tasks
# - Everything in R is 1-indexed (not 0-indexed like Python)
#
# Next: Learn to work with phylogenetic trees using the ape package!
################################################################################
