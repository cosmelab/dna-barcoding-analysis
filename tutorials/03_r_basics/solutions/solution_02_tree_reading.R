#!/usr/bin/env Rscript
################################################################################
# Solution 02: Reading and Exploring Phylogenetic Trees
#
# Complete solutions with explanations
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

library(ape)

cat("=== SOLUTION 02: PHYLOGENETIC TREE ANALYSIS ===\n\n")

# EXERCISE 1: READING TREES =================================================

cat("--- Exercise 1: Reading Trees ---\n\n")

# 1.1 Read tree
tree <- read.tree("../data/example_tree.newick")
cat("Tree successfully loaded\n\n")

# 1.2 Print tree
print(tree)
cat("\n")

# 1.3 Detailed structure
cat("Detailed structure:\n")
str(tree)
cat("\n")

# EXERCISE 2: BASIC TREE PROPERTIES =========================================

cat("--- Exercise 2: Basic Tree Properties ---\n\n")

# 2.1 Number of tips
n_tips <- Ntip(tree)
cat("Number of tips:", n_tips, "\n")

# 2.2 Number of nodes
n_nodes <- tree$Nnode
cat("Number of internal nodes:", n_nodes, "\n")

# 2.3 Is rooted?
cat("Is rooted:", is.rooted(tree), "\n")

# 2.4 Is binary?
cat("Is binary:", is.binary(tree), "\n")

# 2.5 Tip labels
cat("\nTip labels:\n")
print(tree$tip.label)
cat("\n")

# EXERCISE 3: EXPLORING TREE STRUCTURE ======================================

cat("--- Exercise 3: Tree Structure ---\n\n")

# 3.1 Total tree length
total_length <- sum(tree$edge.length)
cat("Total tree length:", round(total_length, 4), "\n")

# 3.2 Longest branch
max_branch <- max(tree$edge.length)
cat("Longest branch:", round(max_branch, 4), "\n")

# 3.3 Shortest branch
min_branch <- min(tree$edge.length)
cat("Shortest branch:", round(min_branch, 4), "\n")

# 3.4 Mean branch length
mean_branch <- mean(tree$edge.length)
cat("Mean branch length:", round(mean_branch, 4), "\n\n")

# 3.5 Histogram of branch lengths
hist(tree$edge.length,
     breaks = 15,
     main = "Branch Length Distribution",
     xlab = "Branch Length (substitutions/site)",
     col = "steelblue",
     border = "white")

# EXERCISE 4: WORKING WITH TIP LABELS =======================================

cat("\n--- Exercise 4: Tip Labels ---\n\n")

# 4.1 Extract genera
genera <- sapply(strsplit(tree$tip.label, "_"), function(x) x[1])
unique_genera <- unique(genera)
cat("Unique genera:", paste(unique_genera, collapse = ", "), "\n")

# 4.2 Count samples per genus
genus_counts <- table(genera)
cat("\nSamples per genus:\n")
print(genus_counts)

# 4.3 Find Aedes tips
aedes_tips <- tree$tip.label[grepl("Aedes", tree$tip.label)]
cat("\nAedes samples:\n")
print(aedes_tips)

# 4.4 Find tips from specific location (CA)
ca_tips <- tree$tip.label[grepl("CA", tree$tip.label)]
cat("\nSamples from California:\n")
print(ca_tips)
cat("\n")

# EXERCISE 5: PHYLOGENETIC DISTANCES ========================================

cat("--- Exercise 5: Phylogenetic Distances ---\n\n")

# 5.1 Distance matrix
dist_matrix <- cophenetic.phylo(tree)
cat("Distance matrix calculated (", nrow(dist_matrix), "x", ncol(dist_matrix), ")\n\n")

# 5.2 Specific distances
cat("Example distances:\n")
cat("First to second tip:", round(dist_matrix[1, 2], 4), "\n")

# Distance between two Aedes aegypti
aedes_aegypti_tips <- tree$tip.label[grepl("Aedes_aegypti", tree$tip.label)]
if (length(aedes_aegypti_tips) >= 2) {
  cat("Between Aedes aegypti samples:",
      round(dist_matrix[aedes_aegypti_tips[1], aedes_aegypti_tips[2]], 4), "\n")
}

# Distance between Aedes and Culex
aedes_sample <- tree$tip.label[grepl("Aedes", tree$tip.label)][1]
culex_sample <- tree$tip.label[grepl("Culex", tree$tip.label)][1]
cat("Aedes to Culex:", round(dist_matrix[aedes_sample, culex_sample], 4), "\n\n")

# 5.3 Most divergent tips
max_dist <- max(dist_matrix)
max_pair <- which(dist_matrix == max_dist, arr.ind = TRUE)[1, ]
cat("Most divergent pair:\n")
cat(" ", tree$tip.label[max_pair[1]], "\n")
cat(" ", tree$tip.label[max_pair[2]], "\n")
cat("  Distance:", round(max_dist, 4), "\n\n")

# 5.4 Most similar tips (excluding self-comparisons)
dist_matrix_copy <- dist_matrix
diag(dist_matrix_copy) <- NA
min_dist <- min(dist_matrix_copy, na.rm = TRUE)
min_pair <- which(dist_matrix_copy == min_dist, arr.ind = TRUE)[1, ]
cat("Most similar pair:\n")
cat(" ", tree$tip.label[min_pair[1]], "\n")
cat(" ", tree$tip.label[min_pair[2]], "\n")
cat("  Distance:", round(min_dist, 4), "\n\n")

# 5.5 Mean pairwise distance
mean_dist <- mean(dist_matrix[lower.tri(dist_matrix)])
cat("Mean pairwise distance:", round(mean_dist, 4), "\n\n")

# EXERCISE 6: ROOTING TREES =================================================

cat("--- Exercise 6: Rooting Trees ---\n\n")

# 6.1 Check if rooted
cat("Original tree is rooted:", is.rooted(tree), "\n\n")

# 6.2 Root by outgroup (Anopheles)
anopheles_tips <- tree$tip.label[grepl("Anopheles", tree$tip.label)]
tree_outgroup <- root(tree, outgroup = anopheles_tips)
cat("Tree rooted using Anopheles as outgroup\n")

# 6.3 Midpoint rooting
tree_midpoint <- midpoint.root(tree)
cat("Tree rooted at midpoint\n\n")

# 6.4 Compare rooting methods
par(mfrow = c(1, 2))
plot(tree_outgroup, main = "Outgroup Rooted", cex = 0.6)
plot(tree_midpoint, main = "Midpoint Rooted", cex = 0.6)
par(mfrow = c(1, 1))

cat("\n")

# EXERCISE 7: TREE MANIPULATION =============================================

cat("--- Exercise 7: Tree Manipulation ---\n\n")

# 7.1 Aedes subtree
aedes_tips <- tree$tip.label[grepl("Aedes", tree$tip.label)]
aedes_tree <- keep.tip(tree, aedes_tips)
cat("Aedes subtree created with", Ntip(aedes_tree), "tips\n")

# 7.2 Remove Culex
culex_tips <- tree$tip.label[grepl("Culex", tree$tip.label)]
tree_no_culex <- drop.tip(tree, culex_tips)
cat("Tree without Culex:", Ntip(tree_no_culex), "tips remaining\n\n")

# 7.3 Ladderize
tree_ladderized <- ladderize(tree)
cat("Tree ladderized\n\n")

# 7.4 Compare original and ladderized
par(mfrow = c(1, 2))
plot(tree, main = "Original", cex = 0.6)
plot(tree_ladderized, main = "Ladderized", cex = 0.6)
par(mfrow = c(1, 1))

cat("\n")

# EXERCISE 8: BOOTSTRAP VALUES ==============================================

cat("--- Exercise 8: Bootstrap Values ---\n\n")

# 8.1 Extract bootstrap values
if (!is.null(tree$node.label)) {
  boot_vals <- tree$node.label
  cat("Bootstrap values extracted\n")

  # 8.2 Convert to numeric
  boot_numeric <- as.numeric(boot_vals)
  boot_numeric <- boot_numeric[!is.na(boot_numeric)]

  # 8.3 Statistics
  cat("Mean bootstrap support:", round(mean(boot_numeric), 1), "%\n")
  cat("Nodes >90% support:", sum(boot_numeric > 90), "\n")
  cat("Nodes <70% support:", sum(boot_numeric < 70), "\n\n")

  # 8.4 Histogram
  hist(boot_numeric,
       breaks = 10,
       main = "Bootstrap Support Distribution",
       xlab = "Bootstrap Value (%)",
       col = "darkorange",
       border = "white",
       xlim = c(0, 100))
  abline(v = 70, col = "red", lty = 2, lwd = 2)
  abline(v = 90, col = "darkgreen", lty = 2, lwd = 2)
} else {
  cat("No bootstrap values in this tree\n")
}

cat("\n")

# EXERCISE 9: NODE DEPTHS AND HEIGHTS =======================================

cat("--- Exercise 9: Node Depths ---\n\n")

# 9.1 Calculate node depths
node_depths <- node.depth.edgelength(tree)
cat("Node depths calculated\n")

# 9.2 Maximum depth
max_depth <- max(node_depths)
cat("Maximum tree depth:", round(max_depth, 4), "\n")

# 9.3 Root-to-tip distances
root_to_tip <- node_depths[1:Ntip(tree)]
names(root_to_tip) <- tree$tip.label
cat("\nRoot-to-tip distances:\n")
print(round(sort(root_to_tip, decreasing = TRUE), 4))

# 9.4 Is ultrametric?
is_ultra <- is.ultrametric(tree)
cat("\nIs ultrametric:", is_ultra, "\n")
if (!is_ultra) {
  cat("Root-to-tip distance variance:", round(var(root_to_tip), 6), "\n")
}

cat("\n")

# EXERCISE 10: WRITING TREES ================================================

cat("--- Exercise 10: Writing Trees ---\n\n")

# 10.1 Save as Newick
write.tree(tree, "output_tree.newick")
cat("Tree saved as output_tree.newick\n")

# 10.2 Save as Nexus
write.nexus(tree, file = "output_tree.nex")
cat("Tree saved as output_tree.nex\n")

# 10.3 Save Aedes subtree
write.tree(aedes_tree, "aedes_subtree.newick")
cat("Aedes subtree saved\n\n")

# EXERCISE 11: ANALYSIS FUNCTIONS ===========================================

cat("--- Exercise 11: Analysis Functions ---\n\n")

# 11.1 Genus analysis function
analyze_genus <- function(tree, genus_name) {
  # Get distance matrix
  dist_mat <- cophenetic.phylo(tree)

  # Identify genus samples
  genera <- sapply(strsplit(tree$tip.label, "_"), function(x) x[1])
  genus_tips <- tree$tip.label[genera == genus_name]

  if (length(genus_tips) == 0) {
    cat("No samples found for genus:", genus_name, "\n")
    return(NULL)
  }

  # Calculate within-genus distances
  within_dist <- c()
  if (length(genus_tips) > 1) {
    within_mat <- dist_mat[genus_tips, genus_tips]
    within_dist <- within_mat[lower.tri(within_mat)]
  }

  # Calculate to other genera
  other_tips <- tree$tip.label[genera != genus_name]
  between_dist <- as.vector(dist_mat[genus_tips, other_tips])

  cat("\n=== Analysis for", genus_name, "===\n")
  cat("Number of samples:", length(genus_tips), "\n")
  if (length(within_dist) > 0) {
    cat("Mean within-genus distance:", round(mean(within_dist), 4), "\n")
  }
  cat("Mean distance to other genera:", round(mean(between_dist), 4), "\n\n")

  return(list(
    n_samples = length(genus_tips),
    within_dist = mean(within_dist),
    between_dist = mean(between_dist)
  ))
}

# Test function
for (genus in unique_genera) {
  analyze_genus(tree, genus)
}

# 11.2 Find outliers
find_outliers <- function(tree, distance_threshold = 0.05) {
  dist_mat <- cophenetic.phylo(tree)
  genera <- sapply(strsplit(tree$tip.label, "_"), function(x) x[1])

  outliers <- c()
  for (i in 1:Ntip(tree)) {
    same_genus <- genera == genera[i]
    same_genus[i] <- FALSE  # Exclude self

    if (sum(same_genus) > 0) {
      mean_dist_to_genus <- mean(dist_mat[i, same_genus])
      if (mean_dist_to_genus > distance_threshold) {
        outliers <- c(outliers, tree$tip.label[i])
      }
    }
  }

  cat("Potential outliers (mean distance to genus >", distance_threshold, "):\n")
  if (length(outliers) > 0) {
    print(outliers)
  } else {
    cat("None detected\n")
  }
  cat("\n")

  return(outliers)
}

find_outliers(tree, distance_threshold = 0.04)

# 11.3 Comprehensive report
tree_summary_report <- function(tree) {
  cat("\n")
  cat("==================================================\n")
  cat("PHYLOGENETIC TREE SUMMARY REPORT\n")
  cat("==================================================\n\n")

  # Basic properties
  cat("--- BASIC PROPERTIES ---\n")
  cat("Tips:", Ntip(tree), "\n")
  cat("Nodes:", tree$Nnode, "\n")
  cat("Rooted:", is.rooted(tree), "\n")
  cat("Binary:", is.binary(tree), "\n\n")

  # Branch statistics
  cat("--- BRANCH LENGTHS ---\n")
  cat("Total:", round(sum(tree$edge.length), 4), "\n")
  cat("Mean:", round(mean(tree$edge.length), 4), "\n")
  cat("Range:", round(min(tree$edge.length), 4), "-",
      round(max(tree$edge.length), 4), "\n\n")

  # Distances
  dist_mat <- cophenetic.phylo(tree)
  all_dist <- dist_mat[lower.tri(dist_mat)]
  cat("--- PAIRWISE DISTANCES ---\n")
  cat("Mean:", round(mean(all_dist), 4), "\n")
  cat("Range:", round(min(all_dist), 4), "-",
      round(max(all_dist), 4), "\n\n")

  # Composition
  genera <- sapply(strsplit(tree$tip.label, "_"), function(x) x[1])
  cat("--- GENUS COMPOSITION ---\n")
  print(table(genera))

  # Bootstrap
  if (!is.null(tree$node.label)) {
    boot_vals <- as.numeric(tree$node.label)
    boot_vals <- boot_vals[!is.na(boot_vals)]
    if (length(boot_vals) > 0) {
      cat("\n--- BOOTSTRAP SUPPORT ---\n")
      cat("Mean:", round(mean(boot_vals), 1), "%\n")
      cat(">90%:", sum(boot_vals > 90), "nodes\n")
    }
  }

  cat("\n==================================================\n\n")
}

tree_summary_report(tree)

cat("=== SOLUTIONS COMPLETE ===\n\n")
cat("Key Takeaways:\n")
cat("1. Trees are 'phylo' objects with specific structure\n")
cat("2. cophenetic.phylo() gives tip-to-tip distances\n")
cat("3. Rooting affects interpretation but not topology\n")
cat("4. Bootstrap values indicate node reliability\n")
cat("5. Always check for outliers and unusual patterns\n")
