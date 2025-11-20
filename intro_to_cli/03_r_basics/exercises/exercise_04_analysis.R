#!/usr/bin/env Rscript
################################################################################
# Exercise 04: Phylogenetic Statistical Analysis
#
# Practice calculating tree statistics and performing analyses
# Apply statistical methods to phylogenetic data
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)

# EXERCISE 1: DISTANCE CALCULATIONS =========================================

# 1.1 Read your tree
# tree <- read.tree("../data/example_tree.newick")

# YOUR CODE HERE:


# 1.2 Calculate the cophenetic distance matrix

# YOUR CODE HERE:


# 1.3 What is the mean pairwise distance across all samples?

# YOUR CODE HERE:


# 1.4 What is the maximum pairwise distance?

# YOUR CODE HERE:


# 1.5 Create a histogram of all pairwise distances

# YOUR CODE HERE:


# EXERCISE 2: WITHIN VS BETWEEN GROUP DISTANCES =============================

# 2.1 Separate samples by genus
# Extract genus from tip labels

# YOUR CODE HERE:


# 2.2 Calculate mean within-genus distances for each genus

# YOUR CODE HERE:


# 2.3 Calculate mean between-genus distances

# YOUR CODE HERE:


# 2.4 Create a boxplot comparing within vs between genus distances

# YOUR CODE HERE:


# 2.5 Perform a statistical test (t-test or Wilcoxon test)
# Are within-genus distances significantly different from between-genus?

# YOUR CODE HERE:


# 2.6 What does this tell you about species boundaries?

# YOUR ANSWER:


# EXERCISE 3: BRANCH LENGTH ANALYSIS ========================================

# 3.1 Extract all branch lengths from the tree

# YOUR CODE HERE:


# 3.2 Calculate summary statistics:
# - Mean, median, SD, min, max

# YOUR CODE HERE:


# 3.3 Test if branch lengths follow a normal distribution
# Use shapiro.test()

# YOUR CODE HERE:


# 3.4 Create a density plot of branch lengths

# YOUR CODE HERE:


# 3.5 Identify outlier branches (>2 SD from mean)

# YOUR CODE HERE:


# EXERCISE 4: BOOTSTRAP SUPPORT ANALYSIS ===================================

# 4.1 Extract bootstrap values from node labels

# YOUR CODE HERE:


# 4.2 Calculate percentage of nodes with:
# - >95% support (very strong)
# - 90-95% support (strong)
# - 70-90% support (moderate)
# - <70% support (weak)

# YOUR CODE HERE:


# 4.3 Create a bar chart showing support categories

# YOUR CODE HERE:


# 4.4 Is the overall tree well-supported?

# YOUR ANSWER:


# EXERCISE 5: ROOT-TO-TIP ANALYSIS ==========================================

# 5.1 Calculate root-to-tip distances for all tips

# YOUR CODE HERE:


# 5.2 Is the tree ultrametric (all tips equidistant from root)?
# Check variance in root-to-tip distances

# YOUR CODE HERE:


# 5.3 Create a barplot of root-to-tip distances

# YOUR CODE HERE:


# 5.4 Which samples are most divergent from the root?

# YOUR CODE HERE:


# EXERCISE 6: GENETIC DIVERSITY METRICS =====================================

# 6.1 Calculate Faith's Phylogenetic Diversity (total tree length)

# YOUR CODE HERE:


# 6.2 Calculate mean pairwise distance (MPD)

# YOUR CODE HERE:


# 6.3 Compare phylogenetic diversity between genera
# (sum of branch lengths for each genus)

# YOUR CODE HERE:


# 6.4 Which genus is most genetically diverse?

# YOUR ANSWER:


# EXERCISE 7: IDENTIFYING OUTLIERS ==========================================

# 7.1 For each sample, calculate its mean distance to other samples
# in the same genus

# YOUR CODE HERE:


# 7.2 Identify samples that are >2 SD away from their genus mean

# YOUR CODE HERE:


# 7.3 These could be:
# - Mislabeled samples
# - Contamination
# - New species
# What would you do to investigate?

# YOUR ANSWER:


# EXERCISE 8: COMPARING GROUPS ==============================================

# 8.1 Create a function that compares two groups (e.g., two genera)
# and returns:
# - Mean within-group distance for each
# - Mean between-group distance
# - Statistical test result

compare_groups <- function(tree, group1_pattern, group2_pattern) {
  # YOUR CODE HERE

}

# Test your function:
# compare_groups(tree, "Aedes", "Culex")


# 8.2 Which pair of genera is most divergent?

# YOUR CODE HERE:


# 8.3 Which pair is least divergent?

# YOUR CODE HERE:


# EXERCISE 9: TEMPORAL ANALYSIS =============================================

# If your tree has sampling dates in tip labels:

# 9.1 Extract years from tip labels

# YOUR CODE HERE:


# 9.2 Test for temporal signal:
# Is there a correlation between sampling date and root-to-tip distance?

# YOUR CODE HERE:


# 9.3 Create a scatter plot of date vs root-to-tip distance

# YOUR CODE HERE:


# 9.4 What does this tell you about evolutionary rate?

# YOUR ANSWER:


# EXERCISE 10: COMPREHENSIVE TREE STATISTICS ================================

# 10.1 Create a function that generates a complete statistical report:
# - Tree properties (N tips, nodes, etc.)
# - Branch length statistics
# - Distance statistics (within/between groups)
# - Bootstrap support summary
# - Diversity metrics
# - Identification of outliers

comprehensive_tree_stats <- function(tree) {
  # YOUR CODE HERE

}

# Test your function


# EXERCISE 11: QUALITY ASSESSMENT ===========================================

# 11.1 Write a function that assesses tree quality:
# Returns TRUE if:
# - >80% of nodes have >70% bootstrap support
# - No branches >5x median branch length (potential errors)
# - Mean within-genus distance < mean between-genus distance

assess_tree_quality <- function(tree) {
  # YOUR CODE HERE

}


# 11.2 What would you do if a tree fails quality assessment?

# YOUR ANSWER:


# EXERCISE 12: DISTANCE-BASED SPECIES DELIMITATION =========================

# 12.1 Calculate the "barcoding gap": the difference between maximum
# within-species distance and minimum between-species distance

# YOUR CODE HERE:


# 12.2 Plot the distribution of within vs between species distances
# Is there a clear gap?

# YOUR CODE HERE:


# 12.3 At what distance threshold would you delimit species?

# YOUR ANSWER:


# EXERCISE 13: VISUALIZATION OF STATISTICS ==================================

# 13.1 Create a multi-panel figure showing:
# - Tree
# - Branch length distribution
# - Pairwise distance distribution
# - Bootstrap support distribution

# YOUR CODE HERE:


# 13.2 Add statistical annotations to your tree:
# - Color tips by their mean distance to others
# - Size nodes by bootstrap support

# YOUR CODE HERE:


# EXERCISE 14: COMPARATIVE ANALYSIS =========================================

# 14.1 If you have two trees (e.g., COI vs 16S), calculate:
# - Robinson-Foulds distance (topological difference)
# - Correlation of branch lengths

# YOUR CODE HERE:


# 14.2 Are the trees concordant or discordant?

# YOUR ANSWER:


# BONUS CHALLENGE ==========================================================

# Create a complete quality control pipeline that:
# 1. Reads a tree file
# 2. Calculates all relevant statistics
# 3. Identifies potential problems:
#    - Low bootstrap support
#    - Outlier sequences
#    - Unexpected relationships
#    - Very long branches
# 4. Generates visualizations
# 5. Produces a PDF report with recommendations

phylogenetic_qc_pipeline <- function(tree_file, output_dir = "qc_results") {
  # YOUR CODE HERE

}


# ADVANCED CHALLENGE: PHYLOGENETIC SIGNAL ==================================

# If you have trait data:

# 15.1 Create trait data for your samples
# (e.g., body size, biting rate, disease transmission)

# YOUR CODE HERE:


# 15.2 Test for phylogenetic signal using Blomberg's K

# YOUR CODE HERE:


# 15.3 Visualize traits on the tree
# Are closely related species similar?

# YOUR CODE HERE:


# 15.4 Perform ancestral state reconstruction

# YOUR CODE HERE:


################################################################################
# REFLECTION QUESTIONS
#
# After completing these exercises, answer these questions:
#
# 1. What is the "barcoding gap" and why is it important?
#
# 2. How do you interpret low bootstrap support values?
#
# 3. What causes some samples to be outliers from their expected group?
#
# 4. Why is it important to compare within-group vs between-group distances?
#
# 5. What statistical tests are appropriate for phylogenetic data?
#
# 6. How would you use these analyses to identify mislabeled samples?
#
################################################################################

# Check your answers with:
# source("../solutions/solution_04_analysis.R")
