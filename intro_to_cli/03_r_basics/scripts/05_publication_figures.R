#!/usr/bin/env Rscript
################################################################################
# Publication-Quality Phylogenetic Figures
# Script 05: Creating Figures for Scientific Papers
#
# Learning Objectives:
# - Design publication-ready phylogenetic figures
# - Combine trees with metadata and annotations
# - Use professional color schemes and layouts
# - Export figures in multiple formats (PDF, PNG, SVG)
# - Follow journal formatting guidelines
# - Create multi-panel figures
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

# Load required libraries
library(ape)           # Phylogenetic analysis
library(ggtree)        # Tree visualization
library(ggplot2)       # Grammar of graphics
library(tidyverse)     # Data manipulation
library(patchwork)     # Combine plots
library(RColorBrewer)  # Color palettes
library(scales)        # Scale functions

# SECTION 1: PUBLICATION STANDARDS ===========================================

# Key requirements for publication figures:
# 1. High resolution (300 DPI minimum)
# 2. Readable labels (8-12 pt font)
# 3. Clear legends and scale bars
# 4. Appropriate color schemes (colorblind-friendly)
# 5. Vector format when possible (PDF, SVG)
# 6. Proper figure dimensions (usually 3.5" or 7" wide)
# 7. All text editable (avoid rasterization)

# Common journal figure widths:
# - Single column: 3.5 inches (89 mm)
# - Double column: 7 inches (178 mm)
# - Full page width: 7.5 inches (190 mm)

# SECTION 2: PREPARE DATA ===================================================

# Load example mosquito tree
mosquito_tree <- read.tree(text = "((((Aedes_aegypti_CA_2023:0.005,Aedes_aegypti_FL_2023:0.005)100:0.010,
                                        (Aedes_albopictus_HI_2023:0.008,Aedes_albopictus_TX_2023:0.008)98:0.012)95:0.020,
                                      ((Culex_pipiens_NY_2023:0.015,Culex_pipiens_IL_2023:0.015)100:0.010,
                                       (Culex_quinquefasciatus_CA_2023:0.013,Culex_quinquefasciatus_TX_2023:0.013)97:0.012)92:0.015)88:0.040,
                                     ((Anopheles_gambiae_Kenya_2022:0.018,Anopheles_gambiae_Lab_2022:0.018)100:0.020,
                                      (Anopheles_stephensi_India_2022:0.016,Anopheles_stephensi_Lab_2022:0.016)99:0.022)85:0.030);")

# Create comprehensive metadata
metadata <- data.frame(
  label = mosquito_tree$tip.label,
  genus = c("Aedes", "Aedes", "Aedes", "Aedes",
            "Culex", "Culex", "Culex", "Culex",
            "Anopheles", "Anopheles", "Anopheles", "Anopheles"),
  species = c("aegypti", "aegypti", "albopictus", "albopictus",
              "pipiens", "pipiens", "quinquefasciatus", "quinquefasciatus",
              "gambiae", "gambiae", "stephensi", "stephensi"),
  location = c("California", "Florida", "Hawaii", "Texas",
               "New York", "Illinois", "California", "Texas",
               "Kenya", "Lab", "India", "Lab"),
  vector_status = c("Major", "Major", "Major", "Major",
                    "Major", "Major", "Major", "Major",
                    "Major", "Major", "Major", "Major"),
  disease = c("Dengue", "Dengue", "Dengue", "Dengue",
              "West Nile", "West Nile", "West Nile", "West Nile",
              "Malaria", "Malaria", "Malaria", "Malaria"),
  year = c(2023, 2023, 2023, 2023, 2023, 2023, 2023, 2023, 2022, 2022, 2022, 2022)
)

# SECTION 3: BASIC PUBLICATION FIGURE ========================================

# Create a clean, professional tree
basic_pub_tree <- ggtree(mosquito_tree) %<+% metadata +
  # Main tree structure
  geom_tree(size = 0.8, color = "gray30") +
  # Tip labels with proper formatting
  geom_tiplab(aes(label = paste(genus, species)),
              size = 3.5,
              fontface = "italic",
              offset = 0.003) +
  # Bootstrap values at nodes
  geom_nodelab(aes(label = label,
                   subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
               size = 2.8,
               color = "black",
               hjust = 1.3,
               vjust = -0.5) +
  # Professional theme with axis
  theme_tree2() +
  # Labels and title
  labs(
    title = "Phylogenetic relationships of mosquito species",
    subtitle = "Based on COI gene sequences",
    x = "Nucleotide substitutions per site",
    caption = "Bootstrap support values shown at nodes (>70%)"
  ) +
  # Adjust axis limits for labels
  xlim(0, 0.35) +
  # Theme customizations
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, hjust = 0),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 9),
    plot.caption = element_text(size = 8, hjust = 0)
  )

print(basic_pub_tree)

# SECTION 4: COLOR-CODED BY TAXONOMY =========================================

# Colorblind-friendly palette
genus_colors <- c(
  "Aedes" = "#E69F00",      # Orange
  "Culex" = "#0072B2",      # Blue
  "Anopheles" = "#009E73"   # Green
)

colored_tree <- ggtree(mosquito_tree) %<+% metadata +
  # Colored branches by genus
  geom_tree(aes(color = genus), size = 1.2) +
  # Tip labels
  geom_tiplab(aes(label = paste(genus, species), color = genus),
              size = 3.5,
              fontface = "italic",
              offset = 0.003) +
  # Bootstrap support
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),
                 size = 2.5,
                 color = "black",
                 fill = "white",
                 shape = 21) +
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) <= 90 & as.numeric(label) > 70),
                 size = 2.5,
                 color = "black",
                 fill = "gray70",
                 shape = 21) +
  # Color scale
  scale_color_manual(
    values = genus_colors,
    name = "Genus",
    labels = c("Aedes", "Anopheles", "Culex")
  ) +
  # Theme
  theme_tree2() +
  labs(
    title = "Mosquito phylogeny colored by genus",
    x = "Substitutions per site"
  ) +
  xlim(0, 0.35) +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(face = "italic"),
    plot.title = element_text(size = 12, face = "bold")
  )

print(colored_tree)

# SECTION 5: TREE WITH HIGHLIGHTED CLADES ====================================

# Identify clade nodes (you may need to adjust these)
# Use: ggtree(tree) + geom_text(aes(label=node), hjust=-.3) to find nodes

highlighted_tree <- ggtree(mosquito_tree) %<+% metadata +
  # Highlight clades
  geom_hilight(node = 15, fill = "#E69F00", alpha = 0.2, extend = 0.032) +
  geom_hilight(node = 18, fill = "#0072B2", alpha = 0.2, extend = 0.032) +
  geom_hilight(node = 21, fill = "#009E73", alpha = 0.2, extend = 0.032) +
  # Tree
  geom_tree(size = 0.8) +
  # Labels
  geom_tiplab(aes(label = paste(genus, species)),
              size = 3.5,
              fontface = "italic",
              offset = 0.003) +
  # Clade labels
  geom_cladelab(node = 15, label = "Aedes", offset = 0.028, fontsize = 4, color = "#E69F00") +
  geom_cladelab(node = 18, label = "Culex", offset = 0.028, fontsize = 4, color = "#0072B2") +
  geom_cladelab(node = 21, label = "Anopheles", offset = 0.028, fontsize = 4, color = "#009E73") +
  # Bootstrap
  geom_nodelab(aes(label = label,
                   subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
               size = 2.5, hjust = 1.3, vjust = -0.5) +
  # Theme
  theme_tree2() +
  labs(
    title = "Mosquito genera form distinct clades",
    x = "Substitutions per site"
  ) +
  xlim(0, 0.40)

print(highlighted_tree)

# SECTION 6: CIRCULAR TREE FOR PRESENTATIONS =================================

circular_tree <- ggtree(mosquito_tree, layout = "circular") %<+% metadata +
  # Colored tree
  geom_tree(aes(color = genus), size = 1.2) +
  # Tip labels
  geom_tiplab(aes(label = species, color = genus),
              size = 3.5,
              fontface = "italic") +
  # Support values
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),
                 size = 3,
                 color = "black") +
  # Colors
  scale_color_manual(values = genus_colors, name = "Genus") +
  # Theme
  theme_tree() +
  labs(title = "Mosquito COI Gene Phylogeny") +
  theme(
    legend.position = "right",
    legend.text = element_text(face = "italic"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

print(circular_tree)

# SECTION 7: TREE WITH METADATA HEATMAP ======================================

# Add quantitative trait data
trait_data <- data.frame(
  label = mosquito_tree$tip.label,
  wing_length_mm = c(2.5, 2.6, 2.2, 2.3, 3.1, 3.2, 3.0, 3.15, 2.8, 2.9, 2.7, 2.75),
  biting_rate = c(8, 7, 9, 8, 6, 5, 6, 7, 9, 8, 9, 9),
  resistance = c(0.8, 0.7, 0.9, 0.85, 0.6, 0.4, 0.5, 0.65, 0.3, 0.2, 0.35, 0.25)
)

# Create tree
p_tree <- ggtree(mosquito_tree) %<+% metadata +
  geom_tree(size = 0.8) +
  geom_tiplab(aes(label = paste(genus, species)),
              size = 3,
              fontface = "italic",
              offset = 0.002) +
  theme_tree2() +
  theme(legend.position = "right")

# Add heatmap
tree_with_heatmap <- gheatmap(p_tree,
                              trait_data[, c("wing_length_mm", "biting_rate", "resistance")],
                              offset = 0.04,
                              width = 0.25,
                              colnames_angle = 90,
                              colnames_offset_y = 0.5,
                              colnames_position = "top",
                              font.size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0.5,
                       name = "Value") +
  labs(title = "Mosquito phylogeny with trait data")

print(tree_with_heatmap)

# SECTION 8: MULTI-PANEL FIGURE ==============================================

# Create multiple related figures for a comprehensive panel

# Panel A: Main phylogeny
panel_a <- ggtree(mosquito_tree) %<+% metadata +
  geom_tree(aes(color = genus), size = 1) +
  geom_tiplab(aes(label = species, color = genus),
              size = 3,
              fontface = "italic",
              offset = 0.002) +
  geom_nodelab(aes(label = label,
                   subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
               size = 2.5, hjust = 1.3, vjust = -0.5) +
  scale_color_manual(values = genus_colors) +
  theme_tree2() +
  labs(title = "A) Phylogenetic tree",
       x = "Substitutions per site") +
  xlim(0, 0.35) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"))

# Panel B: Distance distribution
dist_matrix <- cophenetic.phylo(mosquito_tree)
genera <- sapply(strsplit(mosquito_tree$tip.label, "_"), function(x) x[1])

# Calculate within vs between genus distances
within_dist <- c()
between_dist <- c()
for (i in 1:(length(genera)-1)) {
  for (j in (i+1):length(genera)) {
    if (genera[i] == genera[j]) {
      within_dist <- c(within_dist, dist_matrix[i,j])
    } else {
      between_dist <- c(between_dist, dist_matrix[i,j])
    }
  }
}

dist_df <- data.frame(
  distance = c(within_dist, between_dist),
  category = c(rep("Within genus", length(within_dist)),
               rep("Between genus", length(between_dist)))
)

panel_b <- ggplot(dist_df, aes(x = category, y = distance, fill = category)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Within genus" = "lightblue",
                                "Between genus" = "lightcoral")) +
  labs(title = "B) Genetic distances",
       x = "",
       y = "Pairwise distance") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"))

# Panel C: Bootstrap support distribution
boot_vals <- as.numeric(mosquito_tree$node.label)
boot_vals <- boot_vals[!is.na(boot_vals)]

panel_c <- ggplot(data.frame(bootstrap = boot_vals), aes(x = bootstrap)) +
  geom_histogram(bins = 10, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 70, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 90, linetype = "dashed", color = "darkgreen") +
  labs(title = "C) Bootstrap support",
       x = "Bootstrap value (%)",
       y = "Number of nodes") +
  theme_classic() +
  theme(plot.title = element_text(size = 11, face = "bold"))

# Combine panels using patchwork
multi_panel <- (panel_a | (panel_b / panel_c)) +
  plot_annotation(
    title = "Phylogenetic analysis of mosquito species",
    subtitle = "Based on COI gene sequences",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11))
  )

print(multi_panel)

# SECTION 9: SAVING FIGURES IN PUBLICATION FORMATS ==========================

# Create output directory
dir.create("publication_figures", showWarnings = FALSE)

# Save as PDF (vector format - best for publications)
ggsave("publication_figures/figure1_basic.pdf",
       basic_pub_tree,
       width = 7,
       height = 8,
       units = "in",
       dpi = 300)

# Save as high-resolution PNG (for web/presentations)
ggsave("publication_figures/figure1_basic.png",
       basic_pub_tree,
       width = 7,
       height = 8,
       units = "in",
       dpi = 300)

# Save as SVG (editable in Illustrator/Inkscape)
ggsave("publication_figures/figure1_basic.svg",
       basic_pub_tree,
       width = 7,
       height = 8,
       units = "in")

# Save colored tree
ggsave("publication_figures/figure2_colored.pdf",
       colored_tree,
       width = 7,
       height = 8,
       units = "in",
       dpi = 300)

# Save multi-panel figure
ggsave("publication_figures/figure3_multipanel.pdf",
       multi_panel,
       width = 10,
       height = 8,
       units = "in",
       dpi = 300)

# Save circular tree (good for presentations)
ggsave("publication_figures/figure4_circular.pdf",
       circular_tree,
       width = 8,
       height = 8,
       units = "in",
       dpi = 300)

cat("\nAll figures saved to publication_figures/ directory\n")

# SECTION 10: FIGURE FOR SPECIFIC JOURNAL FORMATS ===========================

# Example: Create figure for Nature-style journal (single column width)
nature_tree <- ggtree(mosquito_tree) %<+% metadata +
  geom_tree(aes(color = genus), size = 0.8) +
  geom_tiplab(aes(label = species, color = genus),
              size = 2.5,
              fontface = "italic",
              offset = 0.002) +
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),
                 size = 1.5, color = "black") +
  scale_color_manual(values = genus_colors) +
  theme_tree2() +
  labs(x = "Substitutions/site") +
  xlim(0, 0.35) +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_blank(),
    legend.text = element_text(size = 7, face = "italic"),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8)
  )

# Save at Nature single-column width (89mm = 3.5 inches)
ggsave("publication_figures/figure_nature_format.pdf",
       nature_tree,
       width = 3.5,
       height = 5,
       units = "in",
       dpi = 300)

# SECTION 11: SUPPLEMENTARY FIGURE (HIGH DETAIL) =============================

# Create detailed supplementary figure with all information
supp_tree <- ggtree(mosquito_tree) %<+% metadata +
  geom_tree(size = 0.8) +
  # Full tip labels with location
  geom_tiplab(aes(label = paste(genus, species, "-", location)),
              size = 2.5,
              fontface = "italic",
              offset = 0.002) +
  # All bootstrap values
  geom_nodelab(aes(label = label),
               size = 2,
               hjust = 1.3,
               vjust = -0.5,
               color = "red") +
  # Branch lengths
  geom_text(aes(x = branch, label = round(branch.length, 4)),
            size = 1.8,
            vjust = -0.5,
            color = "blue") +
  theme_tree2() +
  labs(
    title = "Supplementary Figure S1: Detailed mosquito phylogeny",
    subtitle = "All bootstrap values and branch lengths shown",
    x = "Substitutions per site",
    caption = "Red: bootstrap support; Blue: branch lengths"
  ) +
  xlim(0, 0.40)

# Save supplementary figure
ggsave("publication_figures/supplementary_figure_s1.pdf",
       supp_tree,
       width = 10,
       height = 10,
       units = "in",
       dpi = 300)

# SECTION 12: FIGURE LEGEND TEMPLATE =========================================

# Generate a detailed figure legend
figure_legend <- "
Figure 1. Phylogenetic relationships of mosquito species based on COI gene sequences.

Maximum likelihood phylogenetic tree constructed from COI gene sequences (658 bp) of
12 mosquito samples representing three genera (Aedes, Culex, and Anopheles). The tree
was built using IQ-TREE 2 with automatic model selection (GTR+I+G) and 1000 bootstrap
replicates. Branch lengths represent the number of nucleotide substitutions per site.
Bootstrap support values >70% are shown at nodes. Different genera are shown in
different colors: Aedes (orange), Culex (blue), and Anopheles (green). The tree is
rooted using Anopheles as the outgroup. Scale bar represents 0.05 substitutions per site.
"

cat(figure_legend)
writeLines(figure_legend, "publication_figures/figure1_legend.txt")

# SECTION 13: COMPREHENSIVE PUBLICATION FUNCTION =============================

create_publication_figure <- function(tree_file,
                                     metadata_file = NULL,
                                     output_prefix = "phylogeny",
                                     width = 7,
                                     height = 8,
                                     format = c("pdf", "png", "svg"),
                                     color_by = NULL,
                                     show_bootstrap = TRUE,
                                     bootstrap_cutoff = 70,
                                     title = "Phylogenetic tree") {

  # Read tree
  tree <- read.tree(tree_file)

  # Read metadata if provided
  if (!is.null(metadata_file)) {
    metadata <- read.csv(metadata_file)
    p <- ggtree(tree) %<+% metadata
  } else {
    p <- ggtree(tree)
  }

  # Build plot
  p <- p + geom_tree(size = 0.8)

  # Add tip labels
  if (!is.null(color_by) && !is.null(metadata_file)) {
    p <- p + geom_tiplab(aes_string(color = color_by),
                        size = 3.5,
                        fontface = "italic",
                        offset = 0.003)
  } else {
    p <- p + geom_tiplab(size = 3.5,
                        fontface = "italic",
                        offset = 0.003)
  }

  # Add bootstrap values
  if (show_bootstrap) {
    p <- p + geom_nodelab(
      aes(label = label,
          subset = !is.na(as.numeric(label)) & as.numeric(label) > bootstrap_cutoff),
      size = 2.5,
      hjust = 1.3,
      vjust = -0.5
    )
  }

  # Theme
  p <- p + theme_tree2() +
    labs(title = title,
         x = "Substitutions per site")

  # Save in requested formats
  for (fmt in format) {
    filename <- paste0(output_prefix, ".", fmt)
    ggsave(filename, p,
           width = width,
           height = height,
           units = "in",
           dpi = 300)
    cat("Saved:", filename, "\n")
  }

  return(p)
}

# Example usage:
# fig <- create_publication_figure(
#   tree_file = "../data/example_tree.newick",
#   metadata_file = "../data/tree_metadata.csv",
#   output_prefix = "publication_figures/mosquito_phylogeny",
#   format = c("pdf", "png"),
#   color_by = "genus",
#   title = "Mosquito COI Gene Phylogeny"
# )

# EXERCISES TO TRY ==========================================================
#
# 1. Create a publication figure for your own tree with:
#    - Appropriate title and labels
#    - Bootstrap support values
#    - Color-coded tips by taxonomy or location
#    - Save as PDF at journal-appropriate size
#
# 2. Make a multi-panel figure combining:
#    - Main phylogenetic tree
#    - Distance distribution plot
#    - Bootstrap support histogram
#
# 3. Create both rectangular and circular versions of your tree
#    Which is more appropriate for your data?
#
# 4. Add a heatmap showing trait data next to your tree
#
# 5. Write a complete figure legend following journal guidelines
#    Include: methods, sample size, software, interpretation
#
# 6. Prepare figures in multiple formats for:
#    - Print publication (PDF, 300 DPI)
#    - Web presentation (PNG, 72 DPI)
#    - Poster (PNG, 150 DPI, larger size)
#
################################################################################
# END OF SCRIPT
#
# Summary:
# - Publication figures require high resolution (300 DPI)
# - Use vector formats (PDF, SVG) when possible
# - Follow journal guidelines for width and formatting
# - Use colorblind-friendly palettes
# - Include clear labels, legends, and scale bars
# - Multi-panel figures tell a complete story
# - Always include detailed figure legends
# - Test figures at final print size
#
# You now have complete skills for R-based phylogenetic visualization!
################################################################################
