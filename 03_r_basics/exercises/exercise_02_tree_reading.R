#!/usr/bin/env Rscript
################################################################################
# Exercise 02: Reading and Exploring Phylogenetic Trees
#
# Practice reading trees, examining their structure, and basic manipulation
# Uses the ape package for phylogenetic analysis
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

library(ape)

# EXERCISE 1: READING TREES =================================================

# 1.1 Read the example tree from the data folder
# tree <- read.tree("../data/example_tree.newick")

# YOUR CODE HERE:


# 1.2 Print the tree object to see its structure

# YOUR CODE HERE:


# 1.3 Use str() to examine the detailed structure of the tree

# YOUR CODE HERE:


# EXERCISE 2: BASIC TREE PROPERTIES =========================================

# 2.1 How many tips (sequences/species) are in the tree?

# YOUR CODE HERE:


# 2.2 How many internal nodes does the tree have?

# YOUR CODE HERE:


# 2.3 Is the tree rooted?

# YOUR CODE HERE:


# 2.4 Is the tree binary (fully resolved)?

# YOUR CODE HERE:


# 2.5 Get the tip labels (names of all sequences)

# YOUR CODE HERE:


# EXERCISE 3: EXPLORING TREE STRUCTURE ======================================

# 3.1 What is the total tree length (sum of all branch lengths)?

# YOUR CODE HERE:


# 3.2 What is the longest branch in the tree?

# YOUR CODE HERE:


# 3.3 What is the shortest branch in the tree?

# YOUR CODE HERE:


# 3.4 What is the mean branch length?

# YOUR CODE HERE:


# 3.5 Create a histogram of branch lengths

# YOUR CODE HERE:


# EXERCISE 4: WORKING WITH TIP LABELS =======================================

# 4.1 Extract all unique genera from the tip labels
# Hint: Tip labels are formatted as "Genus_species_location"
# Use strsplit() and sapply()

# YOUR CODE HERE:


# 4.2 How many samples are there for each genus?
# Use table()

# YOUR CODE HERE:


# 4.3 Find all tip labels containing "Aedes"
# Hint: Use grep() or grepl()

# YOUR CODE HERE:


# 4.4 Find all tip labels from a specific location (e.g., "CA" or "USA")

# YOUR CODE HERE:


# EXERCISE 5: PHYLOGENETIC DISTANCES ========================================

# 5.1 Calculate the cophenetic distance matrix
# (tip-to-tip distances through the tree)

# YOUR CODE HERE:


# 5.2 What is the genetic distance between:
# - The first and second tips?


# - Two Aedes aegypti samples?


# - An Aedes sample and a Culex sample?


# 5.3 Which two tips are the most divergent (furthest apart)?
# Hint: Use which.max() on the distance matrix

# YOUR CODE HERE:


# 5.4 Which two tips are the most similar (closest)?
# Exclude diagonal and lower triangle

# YOUR CODE HERE:


# 5.5 Calculate the mean pairwise distance between all tips

# YOUR CODE HERE:


# EXERCISE 6: ROOTING TREES =================================================

# 6.1 Check if your tree is already rooted

# YOUR CODE HERE:


# 6.2 Root the tree using Anopheles as the outgroup
# Hint: Find all Anopheles tips first

# YOUR CODE HERE:


# 6.3 Root the tree at the midpoint

# YOUR CODE HERE:


# 6.4 Compare the two rooting methods by plotting them side by side
# Hint: Use par(mfrow=c(1,2))

# YOUR CODE HERE:


# EXERCISE 7: TREE MANIPULATION =============================================

# 7.1 Create a subtree containing only Aedes samples
# Hint: Use keep.tip()

# YOUR CODE HERE:


# 7.2 Create a tree with one genus removed
# Hint: Use drop.tip()

# YOUR CODE HERE:


# 7.3 Ladderize the tree (organize clades for easier viewing)

# YOUR CODE HERE:


# 7.4 Plot the original and ladderized trees side by side

# YOUR CODE HERE:


# EXERCISE 8: BOOTSTRAP VALUES ==============================================

# 8.1 If your tree has bootstrap values, extract them
# They are stored in tree$node.label

# YOUR CODE HERE:


# 8.2 Convert bootstrap values to numeric

# YOUR CODE HERE:


# 8.3 Calculate summary statistics for bootstrap support:
# - Mean support


# - How many nodes have >90% support?


# - How many nodes have <70% support?


# 8.4 Create a histogram of bootstrap values

# YOUR CODE HERE:


# EXERCISE 9: NODE DEPTHS AND HEIGHTS =======================================

# 9.1 Calculate node depths (distance from root to each node)

# YOUR CODE HERE:


# 9.2 What is the maximum tree depth?

# YOUR CODE HERE:


# 9.3 Get the root-to-tip distance for each tip

# YOUR CODE HERE:


# 9.4 Is the tree ultrametric (all tips at same distance from root)?
# Hint: Check if all root-to-tip distances are equal (within tolerance)

# YOUR CODE HERE:


# EXERCISE 10: WRITING TREES ================================================

# 10.1 Save your tree in Newick format

# YOUR CODE HERE:


# 10.2 Save your tree in Nexus format

# YOUR CODE HERE:


# 10.3 Save only the Aedes subtree you created earlier

# YOUR CODE HERE:


# EXERCISE 11: ANALYSIS CHALLENGE ===========================================

# 11.1 Write a function that takes a tree and genus name,
# and returns summary statistics for that genus:
# - Number of samples
# - Mean within-genus distance
# - Mean distance to other genera

analyze_genus <- function(tree, genus_name) {
  # YOUR CODE HERE

}

# Test your function:
# analyze_genus(tree, "Aedes")


# 11.2 Write a function that identifies potential contamination
# (tips that are very far from their expected genus)

find_outliers <- function(tree, distance_threshold = 0.1) {
  # YOUR CODE HERE

}


# 11.3 Create a comprehensive tree summary report that includes:
# - Basic tree properties
# - Branch length statistics
# - Distance statistics
# - Genus composition
# - Bootstrap support summary

tree_summary_report <- function(tree) {
  # YOUR CODE HERE

}

# Test your function


# EXERCISE 12: COMPARING TREES ==============================================

# 12.1 Read a second tree file (or create a slightly different tree)

# YOUR CODE HERE:


# 12.2 Calculate the Robinson-Foulds distance between the two trees
# This measures topological difference

# YOUR CODE HERE:


# 12.3 Are the trees significantly different?

# YOUR CODE HERE:


# BONUS CHALLENGE ==========================================================

# Write a function that:
# 1. Reads a tree file
# 2. Checks for common problems:
#    - Negative branch lengths
#    - Missing tip labels
#    - Low bootstrap support
#    - Unusual outliers
# 3. Produces a QC report with recommendations

tree_quality_control <- function(tree_file) {
  # YOUR CODE HERE

}

################################################################################
# REFLECTION QUESTIONS
#
# After completing these exercises, answer these questions:
#
# 1. What is the difference between a rooted and unrooted tree?
#
# 2. Why are cophenetic distances useful for DNA barcoding?
#
# 3. What do low bootstrap values indicate about tree nodes?
#
# 4. When should you use midpoint rooting vs outgroup rooting?
#
# 5. How would you identify potential mislabeled or contaminated samples?
#
################################################################################

# Check your answers with:
# source("../solutions/solution_02_tree_reading.R")
