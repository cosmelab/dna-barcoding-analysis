#!/usr/bin/env Rscript
################################################################################
# Solution 04: Phylogenetic Statistical Analysis
#
# Complete solutions with explanations
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)

cat("=== SOLUTION 04: PHYLOGENETIC ANALYSIS ===\n\n")

# Read tree
tree <- read.tree("../data/example_tree.newick")
cat("Tree loaded with", Ntip(tree), "tips\n\n")

# EXERCISE 1: DISTANCE CALCULATIONS =========================================

cat("--- Exercise 1: Distance Calculations ---\n\n")

# 1.2 Distance matrix
dist_matrix <- cophenetic.phylo(tree)
cat("Distance matrix calculated (", nrow(dist_matrix), "x", ncol(dist_matrix), ")\n")

# 1.3 Mean pairwise distance
all_distances <- dist_matrix[lower.tri(dist_matrix)]
mean_dist <- mean(all_distances)
cat("Mean pairwise distance:", round(mean_dist, 4), "\n")

# 1.4 Maximum distance
max_dist <- max(all_distances)
cat("Maximum pairwise distance:", round(max_dist, 4), "\n")

# 1.5 Histogram
hist(all_distances,
     breaks = 20,
     main = "Distribution of Pairwise Distances",
     xlab = "Genetic Distance",
     col = "steelblue",
     border = "white")

cat("\n")

# EXERCISE 2: WITHIN VS BETWEEN GROUP DISTANCES =============================

cat("--- Exercise 2: Within vs Between Group Distances ---\n\n")

# 2.1 Extract genera
genera <- sapply(strsplit(tree$tip.label, "_"), function(x) x[1])

# 2.2 Within-genus distances
within_distances <- list()
for (genus in unique(genera)) {
  genus_tips <- tree$tip.label[genera == genus]
  if (length(genus_tips) > 1) {
    genus_dist <- dist_matrix[genus_tips, genus_tips]
    within_distances[[genus]] <- genus_dist[lower.tri(genus_dist)]
  }
}

cat("Mean within-genus distances:\n")
for (genus in names(within_distances)) {
  cat(" ", genus, ":", round(mean(within_distances[[genus]]), 4), "\n")
}
cat("\n")

# 2.3 Between-genus distances
between_distances <- c()
genus_pairs <- combn(unique(genera), 2)
for (i in 1:ncol(genus_pairs)) {
  g1 <- genus_pairs[1, i]
  g2 <- genus_pairs[2, i]
  g1_tips <- tree$tip.label[genera == g1]
  g2_tips <- tree$tip.label[genera == g2]
  pair_dist <- as.vector(dist_matrix[g1_tips, g2_tips])
  between_distances <- c(between_distances, pair_dist)
}

cat("Mean between-genus distance:", round(mean(between_distances), 4), "\n\n")

# 2.4 Boxplot comparison
all_within <- unlist(within_distances)
comparison_data <- data.frame(
  distance = c(all_within, between_distances),
  category = c(rep("Within genus", length(all_within)),
               rep("Between genus", length(between_distances)))
)

boxplot(distance ~ category, data = comparison_data,
        main = "Within vs Between Genus Distances",
        ylab = "Genetic Distance",
        col = c("lightblue", "lightcoral"))

# 2.5 Statistical test
t_result <- t.test(all_within, between_distances)
cat("T-test results:\n")
cat("  t-statistic:", round(t_result$statistic, 3), "\n")
cat("  p-value:", format(t_result$p.value, scientific = TRUE), "\n")
cat("  Interpretation:", ifelse(t_result$p.value < 0.05,
                                "Significantly different",
                                "Not significantly different"), "\n\n")

cat("Within-genus distances are significantly lower than between-genus,\n")
cat("indicating clear species boundaries (barcoding gap)\n\n")

# EXERCISE 3: BRANCH LENGTH ANALYSIS ========================================

cat("--- Exercise 3: Branch Length Analysis ---\n\n")

# 3.1 Extract branch lengths
branch_lengths <- tree$edge.length

# 3.2 Summary statistics
cat("Branch length statistics:\n")
cat("  Mean:", round(mean(branch_lengths), 4), "\n")
cat("  Median:", round(median(branch_lengths), 4), "\n")
cat("  SD:", round(sd(branch_lengths), 4), "\n")
cat("  Min:", round(min(branch_lengths), 4), "\n")
cat("  Max:", round(max(branch_lengths), 4), "\n\n")

# 3.3 Normality test
shapiro_result <- shapiro.test(branch_lengths)
cat("Shapiro-Wilk test for normality:\n")
cat("  W:", round(shapiro_result$statistic, 4), "\n")
cat("  p-value:", format(shapiro_result$p.value, scientific = TRUE), "\n")
cat("  Interpretation:", ifelse(shapiro_result$p.value < 0.05,
                                "Not normally distributed",
                                "Normally distributed"), "\n\n")

# 3.4 Density plot
plot(density(branch_lengths),
     main = "Branch Length Distribution",
     xlab = "Branch Length",
     col = "darkblue",
     lwd = 2)
polygon(density(branch_lengths), col = rgb(0, 0, 1, 0.2))

# 3.5 Identify outliers
mean_bl <- mean(branch_lengths)
sd_bl <- sd(branch_lengths)
outliers <- branch_lengths[branch_lengths > mean_bl + 2*sd_bl]
cat("\nOutlier branches (>2 SD from mean):\n")
cat("  Count:", length(outliers), "\n")
if (length(outliers) > 0) {
  cat("  Values:", round(outliers, 4), "\n")
}
cat("\n")

# EXERCISE 4: BOOTSTRAP SUPPORT =============================================

cat("--- Exercise 4: Bootstrap Support ---\n\n")

# 4.1 Extract bootstrap values
if (!is.null(tree$node.label)) {
  boot_vals <- as.numeric(tree$node.label)
  boot_vals <- boot_vals[!is.na(boot_vals)]

  # 4.2 Support categories
  very_strong <- sum(boot_vals > 95)
  strong <- sum(boot_vals > 90 & boot_vals <= 95)
  moderate <- sum(boot_vals >= 70 & boot_vals <= 90)
  weak <- sum(boot_vals < 70)

  cat("Bootstrap support categories:\n")
  cat("  Very strong (>95%):", very_strong,
      paste0("(", round(very_strong/length(boot_vals)*100, 1), "%)"), "\n")
  cat("  Strong (90-95%):", strong,
      paste0("(", round(strong/length(boot_vals)*100, 1), "%)"), "\n")
  cat("  Moderate (70-90%):", moderate,
      paste0("(", round(moderate/length(boot_vals)*100, 1), "%)"), "\n")
  cat("  Weak (<70%):", weak,
      paste0("(", round(weak/length(boot_vals)*100, 1), "%)"), "\n\n")

  # 4.3 Bar chart
  support_data <- data.frame(
    category = c("Very Strong\n(>95%)", "Strong\n(90-95%)",
                 "Moderate\n(70-90%)", "Weak\n(<70%)"),
    count = c(very_strong, strong, moderate, weak)
  )

  barplot(support_data$count,
          names.arg = support_data$category,
          main = "Bootstrap Support Distribution",
          ylab = "Number of Nodes",
          col = c("darkgreen", "green", "orange", "red"),
          las = 2)

  # 4.4 Overall assessment
  well_supported <- (very_strong + strong) / length(boot_vals) * 100
  cat("Overall tree support:", round(well_supported, 1), "% of nodes >90%\n")
  cat("Assessment:", ifelse(well_supported > 70,
                           "Well-supported tree",
                           "Poorly supported tree - interpret with caution"), "\n\n")
} else {
  cat("No bootstrap values in this tree\n\n")
}

# EXERCISE 5: ROOT-TO-TIP ANALYSIS ==========================================

cat("--- Exercise 5: Root-to-Tip Analysis ---\n\n")

# 5.1 Calculate distances
node_depths <- node.depth.edgelength(tree)
root_to_tip <- node_depths[1:Ntip(tree)]
names(root_to_tip) <- tree$tip.label

# 5.2 Check ultrametricity
rtt_var <- var(root_to_tip)
is_ultra <- is.ultrametric(tree)
cat("Is tree ultrametric:", is_ultra, "\n")
cat("Root-to-tip variance:", format(rtt_var, scientific = TRUE), "\n")
cat("Coefficient of variation:", round(sd(root_to_tip)/mean(root_to_tip)*100, 2), "%\n\n")

# 5.3 Barplot
barplot(sort(root_to_tip, decreasing = TRUE),
        las = 2,
        main = "Root-to-Tip Distances",
        ylab = "Distance from Root",
        col = "steelblue",
        cex.names = 0.6)

# 5.4 Most divergent
most_divergent <- names(root_to_tip)[which.max(root_to_tip)]
cat("\nMost divergent sample:", most_divergent, "\n")
cat("Distance:", round(max(root_to_tip), 4), "\n\n")

# EXERCISE 6: GENETIC DIVERSITY =============================================

cat("--- Exercise 6: Genetic Diversity Metrics ---\n\n")

# 6.1 Faith's PD
faiths_pd <- sum(tree$edge.length)
cat("Faith's Phylogenetic Diversity:", round(faiths_pd, 4), "\n")

# 6.2 Mean Pairwise Distance
mpd <- mean(all_distances)
cat("Mean Pairwise Distance (MPD):", round(mpd, 4), "\n\n")

# 6.3 Diversity by genus
cat("Phylogenetic diversity by genus:\n")
for (genus in unique(genera)) {
  genus_tips <- tree$tip.label[genera == genus]
  if (length(genus_tips) > 1) {
    genus_tree <- keep.tip(tree, genus_tips)
    genus_pd <- sum(genus_tree$edge.length)
    cat(" ", genus, ":", round(genus_pd, 4), "\n")
  }
}
cat("\n")

# EXERCISE 7: IDENTIFYING OUTLIERS ==========================================

cat("--- Exercise 7: Identifying Outliers ---\n\n")

# Calculate mean distance to same genus for each sample
outlier_analysis <- data.frame(
  sample = tree$tip.label,
  genus = genera,
  mean_to_genus = NA,
  is_outlier = FALSE
)

for (i in 1:Ntip(tree)) {
  same_genus <- genera == genera[i]
  same_genus[i] <- FALSE

  if (sum(same_genus) > 0) {
    mean_dist_genus <- mean(dist_matrix[i, same_genus])
    outlier_analysis$mean_to_genus[i] <- mean_dist_genus

    # Check if outlier (>2 SD from mean of genus)
    genus_tips <- which(genera == genera[i])
    if (length(genus_tips) > 2) {
      genus_dists <- c()
      for (j in genus_tips) {
        if (j != i) {
          genus_dists <- c(genus_dists, mean(dist_matrix[j, genus_tips[genus_tips != j]]))
        }
      }
      threshold <- mean(genus_dists) + 2*sd(genus_dists)
      outlier_analysis$is_outlier[i] <- mean_dist_genus > threshold
    }
  }
}

outliers <- outlier_analysis[outlier_analysis$is_outlier, ]
cat("Potential outliers detected:\n")
if (nrow(outliers) > 0) {
  print(outliers)
  cat("\nThese samples could be:\n")
  cat("  - Mislabeled\n")
  cat("  - Contamination\n")
  cat("  - New/cryptic species\n")
  cat("  - Sequencing errors\n")
  cat("\nRecommended actions:\n")
  cat("  1. Verify sample identification\n")
  cat("  2. Check chromatogram quality\n")
  cat("  3. Re-sequence if possible\n")
  cat("  4. BLAST sequence for confirmation\n")
} else {
  cat("  None detected\n")
}
cat("\n")

# EXERCISE 10: COMPREHENSIVE STATISTICS =====================================

cat("--- Exercise 10: Comprehensive Statistics Function ---\n\n")

comprehensive_tree_stats <- function(tree) {
  cat("\n")
  cat("========================================================\n")
  cat("COMPREHENSIVE PHYLOGENETIC TREE STATISTICS\n")
  cat("========================================================\n\n")

  # Basic properties
  cat("=== TREE PROPERTIES ===\n")
  cat("Number of tips:", Ntip(tree), "\n")
  cat("Number of internal nodes:", tree$Nnode, "\n")
  cat("Is rooted:", is.rooted(tree), "\n")
  cat("Is binary:", is.binary(tree), "\n")
  cat("Is ultrametric:", is.ultrametric(tree), "\n\n")

  # Branch lengths
  cat("=== BRANCH LENGTHS ===\n")
  bl <- tree$edge.length
  cat("Total tree length:", round(sum(bl), 4), "\n")
  cat("Mean:", round(mean(bl), 4), "\n")
  cat("Median:", round(median(bl), 4), "\n")
  cat("SD:", round(sd(bl), 4), "\n")
  cat("Range:", round(min(bl), 4), "-", round(max(bl), 4), "\n\n")

  # Distances
  cat("=== PAIRWISE DISTANCES ===\n")
  dm <- cophenetic.phylo(tree)
  all_d <- dm[lower.tri(dm)]
  cat("Mean:", round(mean(all_d), 4), "\n")
  cat("Median:", round(median(all_d), 4), "\n")
  cat("Range:", round(min(all_d), 4), "-", round(max(all_d), 4), "\n\n")

  # Within vs between group
  genera <- sapply(strsplit(tree$tip.label, "_"), function(x) x[1])
  within_d <- c()
  between_d <- c()

  for (i in 1:(Ntip(tree)-1)) {
    for (j in (i+1):Ntip(tree)) {
      if (genera[i] == genera[j]) {
        within_d <- c(within_d, dm[i,j])
      } else {
        between_d <- c(between_d, dm[i,j])
      }
    }
  }

  cat("=== GROUP DISTANCES ===\n")
  cat("Mean within-genus:", round(mean(within_d), 4), "\n")
  cat("Mean between-genus:", round(mean(between_d), 4), "\n")
  cat("Barcoding gap:", round(mean(between_d) - mean(within_d), 4), "\n\n")

  # Bootstrap
  if (!is.null(tree$node.label)) {
    boot <- as.numeric(tree$node.label)
    boot <- boot[!is.na(boot)]
    cat("=== BOOTSTRAP SUPPORT ===\n")
    cat("Mean:", round(mean(boot), 1), "%\n")
    cat(">90%:", sum(boot > 90), "nodes (", round(sum(boot > 90)/length(boot)*100, 1), "%)\n")
    cat("70-90%:", sum(boot >= 70 & boot <= 90), "nodes\n")
    cat("<70%:", sum(boot < 70), "nodes\n\n")
  }

  # Diversity
  cat("=== DIVERSITY METRICS ===\n")
  cat("Faith's PD:", round(sum(tree$edge.length), 4), "\n")
  cat("Mean pairwise distance:", round(mean(all_d), 4), "\n\n")

  # Composition
  cat("=== TAXONOMIC COMPOSITION ===\n")
  print(table(genera))

  cat("\n========================================================\n\n")

  # Return invisibly
  return(invisible(list(
    n_tips = Ntip(tree),
    total_length = sum(bl),
    mean_pairwise_dist = mean(all_d),
    within_genus_dist = mean(within_d),
    between_genus_dist = mean(between_d),
    barcoding_gap = mean(between_d) - mean(within_d)
  )))
}

# Run comprehensive analysis
stats <- comprehensive_tree_stats(tree)

# EXERCISE 12: SPECIES DELIMITATION =========================================

cat("--- Exercise 12: Barcoding Gap Analysis ---\n\n")

# Calculate barcoding gap
within_max <- max(unlist(within_distances))
between_min <- min(between_distances)
barcoding_gap <- between_min - within_max

cat("Barcoding gap analysis:\n")
cat("  Max within-species distance:", round(within_max, 4), "\n")
cat("  Min between-species distance:", round(between_min, 4), "\n")
cat("  Barcoding gap:", round(barcoding_gap, 4), "\n")

if (barcoding_gap > 0) {
  cat("\nClear barcoding gap exists!\n")
  cat("Suggested threshold:", round((within_max + between_min)/2, 4), "\n")
} else {
  cat("\nNo clear barcoding gap - species delimitation difficult\n")
}

# Plot distribution
hist(c(all_within, between_distances),
     breaks = 30,
     main = "Within vs Between Species Distances",
     xlab = "Genetic Distance",
     col = "gray")
hist(all_within, breaks = 30, col = "lightblue", add = TRUE)
hist(between_distances, breaks = 30, col = "lightcoral", add = TRUE)
legend("topright",
       legend = c("Within species", "Between species"),
       fill = c("lightblue", "lightcoral"))

abline(v = (within_max + between_min)/2, col = "red", lwd = 2, lty = 2)

cat("\n=== SOLUTIONS COMPLETE ===\n\n")
cat("Key Takeaways:\n")
cat("1. Barcoding gap = difference between within/between species distances\n")
cat("2. Bootstrap <70% = low confidence, >90% = high confidence\n")
cat("3. Outliers may indicate mislabeling or cryptic species\n")
cat("4. Faith's PD measures total evolutionary diversity\n")
cat("5. Statistical tests validate observed patterns\n")
cat("6. Comprehensive QC is essential for reliable conclusions\n")
