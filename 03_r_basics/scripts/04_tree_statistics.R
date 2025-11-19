#!/usr/bin/env Rscript
################################################################################
# Phylogenetic Tree Statistics and Analysis
# Script 04: Quantitative Analysis of Phylogenies
#
# Learning Objectives:
# - Calculate phylogenetic distances and diversity metrics
# - Analyze branch lengths and tree topology
# - Assess node support and confidence
# - Compute tree-based diversity indices
# - Compare multiple trees statistically
# - Test evolutionary hypotheses
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

# Load required libraries
library(ape)           # Basic phylogenetic functions
library(phangorn)      # Additional phylogenetic methods
library(ggtree)        # For visualization
library(ggplot2)       # For plotting statistics
library(dplyr)         # For data manipulation

# SECTION 1: CREATING EXAMPLE TREES ==========================================

# Create a mosquito tree with bootstrap values
mosquito_tree <- read.tree(text = "((((Aedes_aegypti_1:0.005,Aedes_aegypti_2:0.005)100:0.010,
                                        (Aedes_albopictus_1:0.008,Aedes_albopictus_2:0.008)98:0.012)95:0.020,
                                      ((Culex_pipiens_1:0.015,Culex_pipiens_2:0.015)100:0.010,
                                       (Culex_quinquefasciatus_1:0.013,Culex_quinquefasciatus_2:0.013)97:0.012)92:0.015)88:0.040,
                                     ((Anopheles_gambiae_1:0.018,Anopheles_gambiae_2:0.018)100:0.020,
                                      (Anopheles_stephensi_1:0.016,Anopheles_stephensi_2:0.016)99:0.022)85:0.030);")

# Visualize the tree
ggtree(mosquito_tree) +
  geom_tiplab(size = 3, fontface = "italic") +
  geom_nodelab(aes(label = label), size = 2.5, color = "red", hjust = 1.3) +
  theme_tree2() +
  ggtitle("Mosquito COI Gene Tree with Bootstrap Support")

# SECTION 2: BASIC TREE STATISTICS ==========================================

# Number of tips (taxa/sequences)
n_tips <- Ntip(mosquito_tree)
cat("Number of tips:", n_tips, "\n")

# Number of internal nodes
n_nodes <- mosquito_tree$Nnode
cat("Number of internal nodes:", n_nodes, "\n")

# Number of edges (branches)
n_edges <- nrow(mosquito_tree$edge)
cat("Number of edges:", n_edges, "\n")

# Check if tree is binary (fully resolved)
is_binary <- is.binary(mosquito_tree)
cat("Is binary tree:", is_binary, "\n")

# Check if tree is rooted
is_rooted <- is.rooted(mosquito_tree)
cat("Is rooted:", is_rooted, "\n")

# Check if tree is ultrametric (all tips equidistant from root)
is_ultra <- is.ultrametric(mosquito_tree)
cat("Is ultrametric:", is_ultra, "\n")

# SECTION 3: BRANCH LENGTH STATISTICS ========================================

# Extract all branch lengths
branch_lengths <- mosquito_tree$edge.length

# Summary statistics
cat("\n=== BRANCH LENGTH STATISTICS ===\n")
cat("Total tree length:", round(sum(branch_lengths), 4), "\n")
cat("Mean branch length:", round(mean(branch_lengths), 4), "\n")
cat("Median branch length:", round(median(branch_lengths), 4), "\n")
cat("Standard deviation:", round(sd(branch_lengths), 4), "\n")
cat("Min branch length:", round(min(branch_lengths), 4), "\n")
cat("Max branch length:", round(max(branch_lengths), 4), "\n")

# Visualize branch length distribution
hist(branch_lengths,
     breaks = 20,
     main = "Distribution of Branch Lengths",
     xlab = "Branch Length (substitutions/site)",
     col = "steelblue",
     border = "white")

# Boxplot
boxplot(branch_lengths,
        main = "Branch Length Distribution",
        ylab = "Branch Length",
        col = "lightblue")

# SECTION 4: PHYLOGENETIC DISTANCES ==========================================

# Cophenetic distance matrix (tip-to-tip distances through the tree)
dist_matrix <- cophenetic.phylo(mosquito_tree)

# View part of the distance matrix
cat("\n=== PAIRWISE DISTANCES (first 4 samples) ===\n")
print(round(dist_matrix[1:4, 1:4], 4))

# Summary of all pairwise distances
all_distances <- dist_matrix[lower.tri(dist_matrix)]
cat("\n=== PAIRWISE DISTANCE STATISTICS ===\n")
cat("Mean pairwise distance:", round(mean(all_distances), 4), "\n")
cat("Median pairwise distance:", round(median(all_distances), 4), "\n")
cat("Min pairwise distance:", round(min(all_distances), 4), "\n")
cat("Max pairwise distance:", round(max(all_distances), 4), "\n")

# Visualize distance distribution
hist(all_distances,
     breaks = 20,
     main = "Distribution of Pairwise Distances",
     xlab = "Genetic Distance",
     col = "darkseagreen",
     border = "white")

# SECTION 5: WITHIN-GROUP VS BETWEEN-GROUP DISTANCES ========================

# Function to calculate within-group and between-group distances
calculate_group_distances <- function(tree, group_pattern) {
  # Get distance matrix
  dist_mat <- cophenetic.phylo(tree)

  # Identify groups
  tips <- tree$tip.label
  groups <- sapply(strsplit(tips, "_"), function(x) x[1])

  # Calculate within-group distances
  within_distances <- c()
  for (group in unique(groups)) {
    group_tips <- tips[groups == group]
    if (length(group_tips) > 1) {
      group_dist <- dist_mat[group_tips, group_tips]
      within_distances <- c(within_distances, group_dist[lower.tri(group_dist)])
    }
  }

  # Calculate between-group distances
  between_distances <- c()
  group_pairs <- combn(unique(groups), 2)
  for (i in 1:ncol(group_pairs)) {
    g1 <- group_pairs[1, i]
    g2 <- group_pairs[2, i]
    g1_tips <- tips[groups == g1]
    g2_tips <- tips[groups == g2]
    between_dist <- as.vector(dist_mat[g1_tips, g2_tips])
    between_distances <- c(between_distances, between_dist)
  }

  return(list(
    within = within_distances,
    between = between_distances
  ))
}

# Calculate distances
group_dist <- calculate_group_distances(mosquito_tree)

# Compare within vs between group distances
cat("\n=== WITHIN-GROUP VS BETWEEN-GROUP DISTANCES ===\n")
cat("Mean within-genus distance:", round(mean(group_dist$within), 4), "\n")
cat("Mean between-genus distance:", round(mean(group_dist$between), 4), "\n")

# Visualize comparison
boxplot(list(Within = group_dist$within, Between = group_dist$between),
        main = "Within vs Between Genus Distances",
        ylab = "Genetic Distance",
        col = c("lightblue", "lightcoral"))

# Statistical test
t_test_result <- t.test(group_dist$within, group_dist$between)
cat("T-test p-value:", t_test_result$p.value, "\n")
cat("Interpretation:", ifelse(t_test_result$p.value < 0.05,
                              "Significantly different",
                              "Not significantly different"), "\n")

# SECTION 6: NODE SUPPORT VALUES ============================================

# Extract bootstrap values
bootstrap_values <- as.numeric(mosquito_tree$node.label)
bootstrap_values <- bootstrap_values[!is.na(bootstrap_values)]

cat("\n=== BOOTSTRAP SUPPORT STATISTICS ===\n")
cat("Number of nodes:", length(bootstrap_values), "\n")
cat("Mean bootstrap support:", round(mean(bootstrap_values), 1), "%\n")
cat("Median bootstrap support:", round(median(bootstrap_values), 1), "%\n")
cat("Nodes with >90% support:", sum(bootstrap_values > 90), "\n")
cat("Nodes with >70% support:", sum(bootstrap_values > 70), "\n")
cat("Nodes with <50% support:", sum(bootstrap_values < 50), "\n")

# Visualize bootstrap distribution
hist(bootstrap_values,
     breaks = 10,
     main = "Distribution of Bootstrap Support Values",
     xlab = "Bootstrap Support (%)",
     col = "darkorange",
     border = "white",
     xlim = c(0, 100))
abline(v = 70, col = "red", lty = 2, lwd = 2)
abline(v = 90, col = "darkgreen", lty = 2, lwd = 2)
legend("topleft",
       legend = c("70% threshold", "90% threshold"),
       col = c("red", "darkgreen"),
       lty = 2,
       lwd = 2)

# SECTION 7: TREE DEPTH AND HEIGHT ==========================================

# Calculate node depths (distance from root)
node_depths <- node.depth.edgelength(mosquito_tree)

cat("\n=== TREE DEPTH STATISTICS ===\n")
cat("Maximum tree depth:", round(max(node_depths), 4), "\n")
cat("Mean tip depth:", round(mean(node_depths[1:Ntip(mosquito_tree)]), 4), "\n")

# Get root-to-tip distances for each tip
root_to_tip <- node_depths[1:Ntip(mosquito_tree)]
names(root_to_tip) <- mosquito_tree$tip.label

# Display root-to-tip distances
cat("\n=== ROOT-TO-TIP DISTANCES ===\n")
print(round(sort(root_to_tip, decreasing = TRUE), 4))

# Visualize root-to-tip variation
barplot(sort(root_to_tip),
        las = 2,
        main = "Root-to-Tip Distances",
        ylab = "Distance from Root",
        col = "steelblue",
        cex.names = 0.6)

# SECTION 8: PHYLOGENETIC DIVERSITY METRICS =================================

# Faith's Phylogenetic Diversity (PD): sum of branch lengths
phylo_diversity <- sum(mosquito_tree$edge.length)
cat("\n=== PHYLOGENETIC DIVERSITY ===\n")
cat("Total phylogenetic diversity (Faith's PD):", round(phylo_diversity, 4), "\n")

# Mean pairwise distance (MPD)
mpd <- mean(all_distances)
cat("Mean pairwise distance (MPD):", round(mpd, 4), "\n")

# SECTION 9: TREE BALANCE ===================================================

# Tree balance measures how symmetrical the tree is
# Colless' Index: measures tree imbalance (0 = perfectly balanced)

# Function to calculate Colless index
colless_index <- function(tree) {
  if (!is.binary(tree)) {
    warning("Tree is not binary, Colless index may not be meaningful")
  }

  n <- Ntip(tree)
  nodes <- (n + 1):(n + tree$Nnode)

  imbalance <- 0
  for (node in nodes) {
    descendants <- prop.part(tree)[[node - n]]
    if (length(descendants) > 0) {
      children <- tree$edge[tree$edge[, 1] == node, 2]
      if (length(children) == 2) {
        n_left <- length(prop.part(tree)[[which(sapply(prop.part(tree), function(x) children[1] %in% x))]])
        n_right <- length(prop.part(tree)[[which(sapply(prop.part(tree), function(x) children[2] %in% x))]])
        imbalance <- imbalance + abs(n_left - n_right)
      }
    }
  }

  return(imbalance)
}

# Note: For more sophisticated balance metrics, use the 'apTreeshape' package

# SECTION 10: COMPARING TREES ===============================================

# Create an alternative tree for comparison
alternative_tree <- read.tree(text = "((((Aedes_aegypti_1:0.005,Aedes_aegypti_2:0.005):0.010,
                                          (Aedes_albopictus_1:0.008,Aedes_albopictus_2:0.008):0.012):0.015,
                                        ((Culex_pipiens_1:0.015,Culex_pipiens_2:0.015):0.010,
                                         (Culex_quinquefasciatus_1:0.013,Culex_quinquefasciatus_2:0.013):0.012):0.020):0.040,
                                       ((Anopheles_gambiae_1:0.018,Anopheles_gambiae_2:0.018):0.020,
                                        (Anopheles_stephensi_1:0.016,Anopheles_stephensi_2:0.016):0.022):0.030);")

# Robinson-Foulds distance (topological difference)
# 0 = identical topology, higher = more different
rf_dist <- treedist(mosquito_tree, alternative_tree, method = "RF")
cat("\n=== TREE COMPARISON ===\n")
cat("Robinson-Foulds distance:", rf_dist, "\n")

# Normalized RF distance (0-1 scale)
max_rf <- 2 * (Ntip(mosquito_tree) - 3)
rf_normalized <- rf_dist / max_rf
cat("Normalized RF distance:", round(rf_normalized, 3), "\n")

# Branch score difference (considers branch lengths)
branch_score <- treedist(mosquito_tree, alternative_tree, method = "score")
cat("Branch score difference:", round(branch_score, 4), "\n")

# SECTION 11: PHYLOGENETIC SIGNAL ===========================================

# Test for phylogenetic signal in a trait
# Phylogenetic signal: tendency for related species to resemble each other

# Create example trait data (e.g., body size)
trait_data <- data.frame(
  species = mosquito_tree$tip.label,
  body_size = c(2.5, 2.6, 2.2, 2.3, 3.1, 3.2, 3.0, 3.15, 2.8, 2.9, 2.7, 2.75)
)

# Convert to named vector
body_size <- trait_data$body_size
names(body_size) <- trait_data$species

# Calculate Blomberg's K (measure of phylogenetic signal)
# K = 1: trait evolves as expected under Brownian motion
# K < 1: less phylogenetic signal than expected
# K > 1: more phylogenetic signal than expected
library(phytools)
k_stat <- phylosig(mosquito_tree, body_size, method = "K", test = TRUE)

cat("\n=== PHYLOGENETIC SIGNAL ===\n")
cat("Blomberg's K:", round(k_stat$K, 3), "\n")
cat("P-value:", round(k_stat$P, 4), "\n")
cat("Interpretation:",
    ifelse(k_stat$P < 0.05,
           "Significant phylogenetic signal detected",
           "No significant phylogenetic signal"), "\n")

# SECTION 12: COMPREHENSIVE TREE ANALYSIS FUNCTION ==========================

analyze_tree_statistics <- function(tree, tree_name = "Phylogenetic Tree") {
  cat("\n")
  cat("================================================================================\n")
  cat("COMPREHENSIVE PHYLOGENETIC TREE ANALYSIS:", tree_name, "\n")
  cat("================================================================================\n\n")

  # Basic properties
  cat("--- BASIC TREE PROPERTIES ---\n")
  cat("Number of tips:", Ntip(tree), "\n")
  cat("Number of internal nodes:", tree$Nnode, "\n")
  cat("Is rooted:", is.rooted(tree), "\n")
  cat("Is binary:", is.binary(tree), "\n")
  cat("Is ultrametric:", is.ultrametric(tree), "\n\n")

  # Branch lengths
  cat("--- BRANCH LENGTH STATISTICS ---\n")
  cat("Total tree length:", round(sum(tree$edge.length), 4), "\n")
  cat("Mean branch length:", round(mean(tree$edge.length), 4), "\n")
  cat("Median branch length:", round(median(tree$edge.length), 4), "\n")
  cat("SD branch length:", round(sd(tree$edge.length), 4), "\n")
  cat("Range:", round(min(tree$edge.length), 4), "-", round(max(tree$edge.length), 4), "\n\n")

  # Distances
  cat("--- PAIRWISE DISTANCE STATISTICS ---\n")
  dist_mat <- cophenetic.phylo(tree)
  all_dist <- dist_mat[lower.tri(dist_mat)]
  cat("Mean pairwise distance:", round(mean(all_dist), 4), "\n")
  cat("Median pairwise distance:", round(median(all_dist), 4), "\n")
  cat("Max pairwise distance:", round(max(all_dist), 4), "\n\n")

  # Bootstrap support (if available)
  if (!is.null(tree$node.label)) {
    boot_vals <- as.numeric(tree$node.label)
    boot_vals <- boot_vals[!is.na(boot_vals)]
    if (length(boot_vals) > 0) {
      cat("--- BOOTSTRAP SUPPORT ---\n")
      cat("Mean support:", round(mean(boot_vals), 1), "%\n")
      cat("Nodes >90%:", sum(boot_vals > 90), "\n")
      cat("Nodes >70%:", sum(boot_vals > 70), "\n")
      cat("Nodes <50%:", sum(boot_vals < 50), "\n\n")
    }
  }

  # Tree depth
  cat("--- TREE DEPTH ---\n")
  depths <- node.depth.edgelength(tree)
  cat("Maximum depth:", round(max(depths), 4), "\n")
  cat("Mean tip depth:", round(mean(depths[1:Ntip(tree)]), 4), "\n\n")

  # Composition
  cat("--- TAXONOMIC COMPOSITION ---\n")
  genera <- sapply(strsplit(tree$tip.label, "_"), function(x) x[1])
  genus_table <- table(genera)
  print(genus_table)

  cat("\n")
  cat("================================================================================\n\n")

  # Return statistics as list
  return(invisible(list(
    n_tips = Ntip(tree),
    n_nodes = tree$Nnode,
    total_length = sum(tree$edge.length),
    mean_branch_length = mean(tree$edge.length),
    mean_pairwise_distance = mean(all_dist),
    genus_composition = genus_table
  )))
}

# Use the function
tree_stats <- analyze_tree_statistics(mosquito_tree, "Mosquito COI Gene Tree")

# SECTION 13: VISUAL SUMMARY OF STATISTICS ==================================

# Create a multi-panel figure summarizing tree statistics
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# 1. Branch length distribution
hist(mosquito_tree$edge.length,
     breaks = 15,
     main = "Branch Lengths",
     xlab = "Length",
     col = "steelblue")

# 2. Pairwise distance distribution
hist(all_distances,
     breaks = 15,
     main = "Pairwise Distances",
     xlab = "Distance",
     col = "darkseagreen")

# 3. Bootstrap support
if (!is.null(mosquito_tree$node.label)) {
  boot_vals <- as.numeric(mosquito_tree$node.label)
  boot_vals <- boot_vals[!is.na(boot_vals)]
  hist(boot_vals,
       breaks = 10,
       main = "Bootstrap Support",
       xlab = "Support (%)",
       col = "darkorange",
       xlim = c(0, 100))
}

# 4. Root-to-tip distances
barplot(sort(root_to_tip, decreasing = TRUE),
        main = "Root-to-Tip Distances",
        ylab = "Distance",
        col = "coral",
        las = 2,
        cex.names = 0.5)

par(mfrow = c(1, 1))

# EXERCISES TO TRY ==========================================================
#
# 1. Load your own phylogenetic tree and calculate:
#    - Total tree length
#    - Mean pairwise distance
#    - Bootstrap support summary
#
# 2. Compare within-species vs between-species distances
#    Are they significantly different?
#
# 3. Create a plot showing:
#    - Branch length distribution
#    - Pairwise distance distribution
#    - Bootstrap support histogram
#    - Root-to-tip distance barplot
#
# 4. Write a function that takes a tree and reports:
#    - Whether it passes quality criteria (e.g., >70% nodes with >70% bootstrap)
#    - Most divergent pair of sequences
#    - Most closely related pair of sequences
#
# 5. If you have trait data, test for phylogenetic signal using Blomberg's K
#
################################################################################
# END OF SCRIPT
#
# Summary:
# - Trees have quantifiable statistical properties
# - Branch lengths and distances are key metrics
# - Bootstrap values indicate node confidence
# - Within-group vs between-group distances reveal population structure
# - Robinson-Foulds distance compares tree topologies
# - Phylogenetic signal tests whether traits are conserved
# - These statistics help interpret biological meaning of trees
#
# Next: Create publication-quality figures for papers!
################################################################################
