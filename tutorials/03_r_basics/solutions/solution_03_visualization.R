#!/usr/bin/env Rscript
################################################################################
# Solution 03: Phylogenetic Tree Visualization
#
# Complete solutions with explanations
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(patchwork)

cat("=== SOLUTION 03: TREE VISUALIZATION ===\n\n")

# Read tree and prepare metadata
tree <- read.tree("../data/example_tree.newick")

metadata <- data.frame(
  label = tree$tip.label,
  genus = sapply(strsplit(tree$tip.label, "_"), function(x) x[1]),
  species = sapply(strsplit(tree$tip.label, "_"), function(x) x[2]),
  location = sapply(strsplit(tree$tip.label, "_"), function(x) x[3]),
  year = as.numeric(sapply(strsplit(tree$tip.label, "_"), function(x) x[4]))
)

cat("Tree and metadata loaded\n\n")

# EXERCISE 1: BASIC PLOTTING ================================================

cat("--- Exercise 1: Basic Plotting ---\n\n")

# 1.1-1.5 Basic plot with all elements
basic_plot <- ggtree(tree) +
  geom_tiplab(size = 3) +
  theme_tree2() +
  labs(title = "Mosquito COI Phylogeny",
       x = "Substitutions per site")

print(basic_plot)
cat("Basic plot created\n\n")

# EXERCISE 2: TREE LAYOUTS ==================================================

cat("--- Exercise 2: Tree Layouts ---\n\n")

# Create different layouts
p_rect <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(size = 2.5) +
  ggtitle("Rectangular")

p_circ <- ggtree(tree, layout = "circular") +
  geom_tiplab(size = 2.5) +
  ggtitle("Circular")

p_fan <- ggtree(tree, layout = "fan", open.angle = 120) +
  geom_tiplab(size = 2.5) +
  ggtitle("Fan")

p_slant <- ggtree(tree, layout = "slanted") +
  geom_tiplab(size = 2.5) +
  ggtitle("Slanted")

# Display all layouts
layout_comparison <- (p_rect | p_circ) / (p_fan | p_slant)
print(layout_comparison)

cat("\nLayout comparison created\n")
cat("Best layout depends on: number of tips, purpose, and audience\n\n")

# EXERCISE 3: STYLING TIPS ==================================================

cat("--- Exercise 3: Styling Tips ---\n\n")

# 3.3 Color by genus
styled_plot <- ggtree(tree) %<+% metadata +
  geom_tiplab(aes(color = genus),
              size = 4,
              fontface = "italic",
              offset = 0.002) +
  scale_color_manual(
    values = c("Aedes" = "#E69F00",
               "Culex" = "#0072B2",
               "Anopheles" = "#009E73")
  ) +
  theme_tree2() +
  labs(title = "Phylogeny with Genus Colors")

print(styled_plot)
cat("Styled plot with colored labels created\n\n")

# 3.4 Add tip points by location
points_plot <- ggtree(tree) %<+% metadata +
  geom_tiplab(size = 3, fontface = "italic", offset = 0.003) +
  geom_tippoint(aes(color = location), size = 3) +
  theme_tree2() +
  theme(legend.position = "right")

print(points_plot)
cat("Plot with tip points created\n\n")

# EXERCISE 4: BOOTSTRAP VALUES ==============================================

cat("--- Exercise 4: Bootstrap Values ---\n\n")

# Display bootstrap with color coding
bootstrap_plot <- ggtree(tree) +
  geom_tiplab(size = 3, fontface = "italic") +
  # Green points for >90
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),
                 color = "darkgreen", size = 3) +
  # Orange for 70-90
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) &
                      as.numeric(label) <= 90 & as.numeric(label) > 70),
                 color = "orange", size = 3) +
  # Red for <70
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) <= 70),
                 color = "red", size = 3) +
  theme_tree2() +
  labs(title = "Bootstrap Support",
       subtitle = "Green >90%, Orange 70-90%, Red <70%")

print(bootstrap_plot)
cat("Bootstrap visualization created\n\n")

# EXERCISE 5: METADATA ======================================================

cat("--- Exercise 5: Working with Metadata ---\n\n")

# Comprehensive metadata plot
metadata_plot <- ggtree(tree) %<+% metadata +
  geom_tree(aes(color = genus), size = 1) +
  geom_tiplab(aes(label = species),
              size = 3,
              fontface = "italic",
              offset = 0.003) +
  geom_tippoint(aes(shape = location), size = 3) +
  scale_color_manual(
    values = c("Aedes" = "#E69F00",
               "Culex" = "#0072B2",
               "Anopheles" = "#009E73"),
    name = "Genus"
  ) +
  theme_tree2() +
  theme(legend.position = "right") +
  labs(title = "Phylogeny with Comprehensive Metadata")

print(metadata_plot)
cat("Metadata plot created\n\n")

# EXERCISE 6: HIGHLIGHTING CLADES ===========================================

cat("--- Exercise 6: Highlighting Clades ---\n\n")

# First, identify nodes (you may need to adjust these)
node_plot <- ggtree(tree) +
  geom_tiplab(size = 2) +
  geom_text(aes(label = node), hjust = -0.3, size = 2.5, color = "red") +
  theme_tree2()

print(node_plot)
cat("\nNode numbers displayed - use these to identify clades\n\n")

# Example with hypothetical node numbers (adjust based on your tree)
# Assuming Aedes clade is node 18, Culex is 20, Anopheles is 24
clade_plot <- ggtree(tree) %<+% metadata +
  geom_hilight(node = 18, fill = "#E69F00", alpha = 0.2, extend = 0.03) +
  geom_hilight(node = 20, fill = "#0072B2", alpha = 0.2, extend = 0.03) +
  geom_hilight(node = 24, fill = "#009E73", alpha = 0.2, extend = 0.03) +
  geom_tiplab(size = 3, fontface = "italic", offset = 0.002) +
  geom_cladelab(node = 18, label = "Aedes", offset = 0.025, fontsize = 4) +
  geom_cladelab(node = 20, label = "Culex", offset = 0.025, fontsize = 4) +
  geom_cladelab(node = 24, label = "Anopheles", offset = 0.025, fontsize = 4) +
  theme_tree2()

print(clade_plot)
cat("Clade highlighting complete\n\n")

# EXERCISE 9: PUBLICATION FIGURE ============================================

cat("--- Exercise 9: Publication-Quality Figure ---\n\n")

pub_figure <- ggtree(tree) %<+% metadata +
  geom_tree(aes(color = genus), size = 1) +
  geom_tiplab(aes(label = paste(genus, species), color = genus),
              size = 3.5,
              fontface = "italic",
              offset = 0.003) +
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),
                 size = 2.5, color = "black") +
  scale_color_manual(
    values = c("Aedes" = "#E69F00",
               "Culex" = "#0072B2",
               "Anopheles" = "#009E73"),
    name = "Genus"
  ) +
  theme_tree2() +
  labs(
    title = "Phylogenetic relationships of mosquito species",
    subtitle = "Based on COI gene sequences",
    x = "Substitutions per site",
    caption = "Black circles indicate >90% bootstrap support"
  ) +
  xlim(0, 0.35) +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(face = "italic"),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

print(pub_figure)
cat("Publication figure created\n\n")

# EXERCISE 10: MULTI-PANEL ==================================================

cat("--- Exercise 10: Multi-Panel Figure ---\n\n")

# Panel 1: Rectangular tree
panel1 <- ggtree(tree) %<+% metadata +
  geom_tree(aes(color = genus)) +
  geom_tiplab(size = 2.5, fontface = "italic") +
  scale_color_manual(values = c("Aedes" = "#E69F00",
                                 "Culex" = "#0072B2",
                                 "Anopheles" = "#009E73")) +
  theme_tree2() +
  ggtitle("A) Rectangular Layout") +
  theme(legend.position = "none")

# Panel 2: Circular tree
panel2 <- ggtree(tree, layout = "circular") %<+% metadata +
  geom_tree(aes(color = genus)) +
  geom_tiplab(size = 2.5, fontface = "italic") +
  scale_color_manual(values = c("Aedes" = "#E69F00",
                                 "Culex" = "#0072B2",
                                 "Anopheles" = "#009E73")) +
  ggtitle("B) Circular Layout") +
  theme(legend.position = "none")

# Combine
multi_panel <- panel1 | panel2
print(multi_panel)
cat("Multi-panel figure created\n\n")

# EXERCISE 11: HEATMAP ======================================================

cat("--- Exercise 11: Heatmap ---\n\n")

# Create trait data
trait_data <- data.frame(
  label = tree$tip.label,
  wing_length = rnorm(Ntip(tree), mean = 3, sd = 0.3),
  biting_rate = sample(5:10, Ntip(tree), replace = TRUE),
  resistance = runif(Ntip(tree), min = 0.2, max = 0.9)
)

# Create tree
p_tree <- ggtree(tree) +
  geom_tiplab(size = 3, fontface = "italic", offset = 0.002) +
  theme_tree2()

# Add heatmap
heatmap_plot <- gheatmap(p_tree,
                         trait_data[, c("wing_length", "biting_rate", "resistance")],
                         offset = 0.04,
                         width = 0.25,
                         colnames_angle = 90,
                         colnames_offset_y = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0.5) +
  labs(title = "Phylogeny with Trait Data")

print(heatmap_plot)
cat("Heatmap created\n\n")

# EXERCISE 12: SAVING FIGURES ===============================================

cat("--- Exercise 12: Saving Figures ---\n\n")

# Create output directory
dir.create("figures", showWarnings = FALSE)

# Save publication figure in multiple formats
ggsave("figures/publication_tree.pdf",
       pub_figure,
       width = 7, height = 8, units = "in", dpi = 300)

ggsave("figures/publication_tree.png",
       pub_figure,
       width = 7, height = 8, units = "in", dpi = 300)

ggsave("figures/publication_tree.svg",
       pub_figure,
       width = 7, height = 8, units = "in")

# Single column version
ggsave("figures/publication_tree_single_col.pdf",
       pub_figure,
       width = 3.5, height = 5, units = "in", dpi = 300)

cat("Figures saved in multiple formats\n\n")

# EXERCISE 14: COMPREHENSIVE FUNCTION =======================================

cat("--- Exercise 14: Complete Figure Function ---\n\n")

create_complete_figure <- function(tree_file,
                                   metadata_file = NULL,
                                   trait_file = NULL,
                                   output_prefix = "phylogeny") {

  # Read tree
  tree <- read.tree(tree_file)

  # Read metadata if provided
  if (!is.null(metadata_file)) {
    metadata <- read.csv(metadata_file)
  } else {
    # Create basic metadata from tip labels
    metadata <- data.frame(
      label = tree$tip.label,
      genus = sapply(strsplit(tree$tip.label, "_"), function(x) x[1])
    )
  }

  # Main phylogeny
  p1 <- ggtree(tree) %<+% metadata +
    geom_tree(aes(color = genus), size = 1) +
    geom_tiplab(size = 3, fontface = "italic", offset = 0.002) +
    theme_tree2() +
    ggtitle("A) Phylogenetic Tree") +
    theme(legend.position = "right")

  # Bootstrap histogram
  if (!is.null(tree$node.label)) {
    boot_vals <- as.numeric(tree$node.label)
    boot_vals <- boot_vals[!is.na(boot_vals)]

    p2 <- ggplot(data.frame(bootstrap = boot_vals), aes(x = bootstrap)) +
      geom_histogram(bins = 10, fill = "steelblue", color = "white") +
      geom_vline(xintercept = 70, linetype = "dashed", color = "red") +
      labs(title = "B) Bootstrap Support",
           x = "Bootstrap (%)", y = "Count") +
      theme_classic()
  } else {
    p2 <- ggplot() + ggtitle("B) No Bootstrap Data")
  }

  # Distance distribution
  dist_mat <- cophenetic.phylo(tree)
  all_dist <- dist_mat[lower.tri(dist_mat)]

  p3 <- ggplot(data.frame(distance = all_dist), aes(x = distance)) +
    geom_histogram(bins = 20, fill = "darkseagreen", color = "white") +
    labs(title = "C) Distance Distribution",
         x = "Pairwise Distance", y = "Count") +
    theme_classic()

  # Combine
  combined <- p1 | (p2 / p3)

  # Save
  ggsave(paste0(output_prefix, ".pdf"),
         combined, width = 12, height = 8, dpi = 300)

  cat("Complete figure saved as", paste0(output_prefix, ".pdf\n"))

  return(combined)
}

# Test function
fig <- create_complete_figure("../data/example_tree.newick",
                             output_prefix = "figures/complete_analysis")

print(fig)

cat("\n=== SOLUTIONS COMPLETE ===\n\n")
cat("Key Takeaways:\n")
cat("1. ggtree uses ggplot2 grammar - build plots layer by layer\n")
cat("2. Metadata joins with %<+% operator\n")
cat("3. Different layouts serve different purposes\n")
cat("4. Color schemes should be accessible (colorblind-friendly)\n")
cat("5. Publication figures need: scale, support values, clear labels\n")
cat("6. Save as PDF for publications, PNG for presentations\n")
