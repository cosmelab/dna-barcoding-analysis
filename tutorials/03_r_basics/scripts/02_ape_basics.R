#!/usr/bin/env Rscript
################################################################################
# APE Package Basics for Phylogenetic Analysis
# Script 02: Reading, Manipulating, and Plotting Trees
#
# Learning Objectives:
# - Install and load the ape package
# - Read phylogenetic trees from Newick/Nexus files
# - Understand tree structure in R
# - Manipulate trees (root, prune, extract clades)
# - Create basic tree plots
# - Calculate phylogenetic distances
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

# SECTION 1: INSTALLING AND LOADING APE ======================================

# Install ape (only need to do this once)
# Uncomment the line below to install
# install.packages("ape")

# Load the library (do this every R session)
library(ape)

# Check version
packageVersion("ape")

# SECTION 2: PHYLOGENETIC TREE BASICS ========================================

# What is a phylogenetic tree?
# - Hypothesis of evolutionary relationships
# - Nodes: represent common ancestors
# - Tips (leaves): represent taxa (species, sequences)
# - Branches (edges): represent evolutionary time or change

# In R, trees are stored as "phylo" objects
# A phylo object contains:
# - edge: matrix defining parent-child relationships
# - edge.length: vector of branch lengths
# - tip.label: vector of tip names
# - Nnode: number of internal nodes

# SECTION 3: CREATING SIMPLE TREES ===========================================

# Create a simple tree manually (for learning)
simple_tree <- read.tree(text = "(((Aedes_aegypti:0.1,Aedes_albopictus:0.1):0.2,Culex_pipiens:0.3):0.1,Anopheles_gambiae:0.4);")

# Plot the tree
plot(simple_tree, main = "Simple Mosquito Tree")

# Examine the tree structure
print(simple_tree)
str(simple_tree)  # Detailed structure

# Access tree components
simple_tree$tip.label        # Names of tips
simple_tree$edge            # Edge matrix (parent-child pairs)
simple_tree$edge.length     # Branch lengths
simple_tree$Nnode           # Number of internal nodes

# Get tree information
Ntip(simple_tree)           # Number of tips (species)
is.rooted(simple_tree)      # Is the tree rooted?
is.ultrametric(simple_tree) # Are all tips equidistant from root?

# SECTION 4: READING TREES FROM FILES ========================================

# Most commonly, you'll read trees from files
# Two main formats: Newick (.tree, .treefile, .newick) and Nexus (.nex)

# Read Newick format tree
# tree <- read.tree("../data/example_tree.newick")

# Read multiple trees from one file
# trees <- read.tree("multiple_trees.newick")  # Returns a list of trees

# Read Nexus format
# tree <- read.nexus("tree.nex")

# For this demonstration, let's create a more realistic mosquito tree
mosquito_tree <- read.tree(text = "((((Aedes_aegypti_1:0.01,Aedes_aegypti_2:0.01):0.02,
                                        (Aedes_albopictus_1:0.015,Aedes_albopictus_2:0.015):0.015):0.05,
                                      ((Culex_pipiens_1:0.02,Culex_pipiens_2:0.02):0.03,
                                       (Culex_quinquefasciatus_1:0.025,Culex_quinquefasciatus_2:0.025):0.025):0.02):0.1,
                                     ((Anopheles_gambiae_1:0.03,Anopheles_gambiae_2:0.03):0.05,
                                      (Anopheles_stephensi_1:0.035,Anopheles_stephensi_2:0.035):0.045):0.05);")

# Plot this more complex tree
plot(mosquito_tree,
     main = "COI Gene Tree - Mosquito Species",
     cex = 0.8)

# SECTION 5: BASIC TREE PLOTTING =============================================

# ape provides flexible plotting options

# Different tree layouts
par(mfrow = c(2, 2))  # Create 2x2 panel

# Phylogram (default - shows branch lengths)
plot(mosquito_tree, main = "Phylogram", cex = 0.6)

# Cladogram (ignore branch lengths)
plot(mosquito_tree, type = "cladogram", main = "Cladogram", cex = 0.6)

# Unrooted tree
plot(mosquito_tree, type = "unrooted", main = "Unrooted", cex = 0.6)

# Fan/radial tree
plot(mosquito_tree, type = "fan", main = "Fan Tree", cex = 0.6)

par(mfrow = c(1, 1))  # Reset to single panel

# SECTION 6: CUSTOMIZING TREE PLOTS ==========================================

# Control tip labels
plot(mosquito_tree,
     show.tip.label = TRUE,    # Show tip labels (default)
     cex = 0.7,                # Label size
     font = 2,                 # Bold font
     tip.color = "darkblue",   # Label color
     main = "Customized Tree")

# Add scale bar
add.scale.bar()

# Add node labels (node numbers)
nodelabels(cex = 0.5, bg = "lightblue")

# Add tip labels with numbers
tiplabels(cex = 0.5, bg = "lightgreen")

# Color branches by clade
plot(mosquito_tree, edge.color = "gray50", edge.width = 2)

# Highlight specific tips
plot(mosquito_tree, cex = 0.7)
tips_to_highlight <- which(grepl("Aedes", mosquito_tree$tip.label))
tiplabels(pch = 19, col = "red", cex = 1.5,
          tip = tips_to_highlight)

# SECTION 7: ROOTING TREES ===================================================

# Many phylogenetic analyses require rooted trees
# Root by specifying an outgroup

# Check if tree is rooted
is.rooted(mosquito_tree)

# Root at the midpoint (equal distance to all tips)
midpoint_rooted <- midpoint.root(mosquito_tree)
plot(midpoint_rooted, main = "Midpoint Rooted Tree", cex = 0.7)

# Root by outgroup (Anopheles species are outgroup to Aedes/Culex)
outgroup_tips <- mosquito_tree$tip.label[grepl("Anopheles", mosquito_tree$tip.label)]
outgroup_rooted <- root(mosquito_tree, outgroup = outgroup_tips)
plot(outgroup_rooted, main = "Outgroup Rooted Tree", cex = 0.7)

# Unroot a tree
unrooted_tree <- unroot(outgroup_rooted)
is.rooted(unrooted_tree)  # FALSE

# SECTION 8: MANIPULATING TREES ==============================================

# Drop specific tips
tree_no_aedes <- drop.tip(mosquito_tree,
                          mosquito_tree$tip.label[grepl("Aedes", mosquito_tree$tip.label)])
plot(tree_no_aedes, main = "Tree with Aedes Removed", cex = 0.7)

# Keep only specific tips
culex_only <- keep.tip(mosquito_tree,
                       mosquito_tree$tip.label[grepl("Culex", mosquito_tree$tip.label)])
plot(culex_only, main = "Culex Species Only", cex = 0.8)

# Extract a clade (subtree)
# First, identify the node number
plot(mosquito_tree, cex = 0.6)
nodelabels(cex = 0.6, bg = "yellow")

# Extract clade from specific node (example: node 15)
# Note: actual node numbers depend on tree structure
# subtree <- extract.clade(mosquito_tree, node = 15)
# plot(subtree, main = "Extracted Clade")

# SECTION 9: TREE STATISTICS =================================================

# Calculate phylogenetic distances between all tips
# This creates a distance matrix
dist_matrix <- cophenetic.phylo(mosquito_tree)
print(round(dist_matrix[1:5, 1:5], 4))  # Show first 5x5

# Get distance between two specific tips
dist_matrix["Aedes_aegypti_1", "Aedes_albopictus_1"]
dist_matrix["Aedes_aegypti_1", "Culex_pipiens_1"]
dist_matrix["Aedes_aegypti_1", "Anopheles_gambiae_1"]

# Tree-wide statistics
# Total tree length (sum of all branches)
sum(mosquito_tree$edge.length)

# Mean branch length
mean(mosquito_tree$edge.length)

# Number of tips and nodes
cat("Number of tips:", Ntip(mosquito_tree), "\n")
cat("Number of nodes:", mosquito_tree$Nnode, "\n")

# SECTION 10: COMPARING TREES ================================================

# Create a slightly different tree for comparison
alternative_tree <- read.tree(text = "((((Aedes_aegypti_1:0.01,Aedes_aegypti_2:0.01):0.02,
                                          (Aedes_albopictus_1:0.015,Aedes_albopictus_2:0.015):0.015):0.04,
                                        (Culex_pipiens_1:0.02,Culex_pipiens_2:0.02):0.05):0.1,
                                       ((Anopheles_gambiae_1:0.03,Anopheles_gambiae_2:0.03):0.05,
                                        (Anopheles_stephensi_1:0.035,Anopheles_stephensi_2:0.035):0.045):0.05);")

# Compare tree topologies (Robinson-Foulds distance)
# 0 means identical topology, higher means more different
rf_distance <- treedist(mosquito_tree, alternative_tree, method = "RF")
cat("Robinson-Foulds distance:", rf_distance, "\n")

# SECTION 11: WORKING WITH BOOTSTRAP VALUES ==================================

# Bootstrap values indicate support for each node
# Often stored as node labels

# Create tree with bootstrap values
bootstrap_tree <- read.tree(text = "((Aedes_aegypti:0.1,Aedes_albopictus:0.1)95:0.2,(Culex_pipiens:0.15,Anopheles_gambiae:0.15)78:0.15);")

# Plot with bootstrap values
plot(bootstrap_tree, main = "Tree with Bootstrap Support")
nodelabels(bootstrap_tree$node.label, cex = 0.8, bg = "lightblue")

# Filter nodes by bootstrap support
# Get bootstrap values as numeric
boot_values <- as.numeric(bootstrap_tree$node.label)
cat("Nodes with >90% bootstrap support:", sum(boot_values > 90, na.rm = TRUE), "\n")

# SECTION 12: WRITING TREES TO FILES =========================================

# Save tree in Newick format
write.tree(mosquito_tree, "mosquito_phylogeny.newick")

# Save tree in Nexus format
write.nexus(mosquito_tree, file = "mosquito_phylogeny.nex")

# Save multiple trees
# write.tree(list(tree1, tree2, tree3), "multiple_trees.newick")

# SECTION 13: PRACTICAL EXAMPLE - ANALYZING YOUR COI TREE ===================

# Function to analyze a phylogenetic tree
analyze_phylogeny <- function(tree, tree_name = "Phylogenetic Tree") {
  cat("\n========================================\n")
  cat("PHYLOGENETIC TREE ANALYSIS:", tree_name, "\n")
  cat("========================================\n\n")

  # Basic information
  cat("Number of sequences (tips):", Ntip(tree), "\n")
  cat("Number of internal nodes:", tree$Nnode, "\n")
  cat("Is rooted:", is.rooted(tree), "\n")
  cat("Is ultrametric:", is.ultrametric(tree), "\n\n")

  # Branch length statistics
  cat("--- Branch Length Statistics ---\n")
  cat("Total tree length:", round(sum(tree$edge.length), 4), "\n")
  cat("Mean branch length:", round(mean(tree$edge.length), 4), "\n")
  cat("Median branch length:", round(median(tree$edge.length), 4), "\n")
  cat("Min branch length:", round(min(tree$edge.length), 4), "\n")
  cat("Max branch length:", round(max(tree$edge.length), 4), "\n\n")

  # Distance statistics
  cat("--- Pairwise Distance Statistics ---\n")
  dist_mat <- cophenetic.phylo(tree)
  cat("Mean pairwise distance:", round(mean(dist_mat[lower.tri(dist_mat)]), 4), "\n")
  cat("Max pairwise distance:", round(max(dist_mat), 4), "\n\n")

  # Genus composition
  cat("--- Sample Composition ---\n")
  genera <- sapply(strsplit(tree$tip.label, "_"), function(x) x[1])
  genus_table <- table(genera)
  print(genus_table)

  return(list(
    n_tips = Ntip(tree),
    n_nodes = tree$Nnode,
    total_length = sum(tree$edge.length),
    mean_distance = mean(dist_mat[lower.tri(dist_mat)]),
    genus_counts = genus_table
  ))
}

# Analyze our mosquito tree
results <- analyze_phylogeny(mosquito_tree, "Mosquito COI Gene Tree")

# SECTION 14: LADDERIZING TREES =============================================

# Ladderizing makes trees easier to read by rotating nodes
# so that clades with more tips are on one side

par(mfrow = c(1, 2))

# Original tree
plot(mosquito_tree, main = "Original Tree", cex = 0.6)

# Ladderized tree (right)
ladderized_right <- ladderize(mosquito_tree, right = TRUE)
plot(ladderized_right, main = "Ladderized (Right)", cex = 0.6)

par(mfrow = c(1, 1))

# SECTION 15: COMPUTING CONSENSUS TREES ======================================

# When you have multiple trees (e.g., from bootstrap or Bayesian analysis)
# you can compute a consensus tree

# Create some example trees (variations on mosquito tree)
tree_list <- list(mosquito_tree, alternative_tree)

# Consensus tree (majority-rule: keep splits in >50% of trees)
consensus <- consensus(tree_list, p = 0.5)
plot(consensus, main = "Consensus Tree", cex = 0.7)

# Strict consensus (keep only splits in 100% of trees)
strict_consensus <- consensus(tree_list, p = 1)
plot(strict_consensus, main = "Strict Consensus Tree", cex = 0.7)

# EXERCISES TO TRY ==========================================================
#
# 1. Read the example tree from ../data/example_tree.newick
#    (you'll create this file in the next section)
#
# 2. Root the tree using Anopheles as the outgroup
#
# 3. Create a plot showing:
#    - Tree with tip labels colored by genus
#    - Scale bar
#    - Node labels showing bootstrap support
#
# 4. Calculate the genetic distance between:
#    - Two Aedes species
#    - Aedes and Culex
#    - Aedes and Anopheles
#
# 5. Extract only the Culex clade and plot it separately
#
# 6. Write a function that takes a tree and returns:
#    - Number of each genus
#    - Average branch length
#    - Total tree depth
#
################################################################################
# END OF SCRIPT
#
# Summary:
# - ape is the foundational package for phylogenetics in R
# - Trees are "phylo" objects with edges, tip labels, and branch lengths
# - read.tree() and read.nexus() read tree files
# - plot() provides flexible tree visualization
# - Trees can be rooted, pruned, and manipulated
# - cophenetic.phylo() calculates pairwise distances
# - write.tree() saves trees to files
#
# Next: Learn to create beautiful visualizations with ggtree!
################################################################################
