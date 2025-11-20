#!/usr/bin/env Rscript
################################################################################
# ggtree: Beautiful Phylogenetic Tree Visualizations
# Script 03: Creating Publication-Quality Tree Figures
#
# Learning Objectives:
# - Install and use ggtree package
# - Understand the grammar of graphics for trees
# - Create various tree layouts (rectangular, circular, fan)
# - Add metadata and annotations
# - Color and style trees for publications
# - Combine multiple tree visualizations
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

# SECTION 1: INSTALLING AND LOADING GGTREE ===================================

# ggtree is part of Bioconductor, not CRAN
# Install ggtree (only once)
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ggtree")

# Load required libraries
library(ape)        # For reading trees
library(ggtree)     # For beautiful tree plots
library(ggplot2)    # ggplot2 features work with ggtree
library(dplyr)      # For data manipulation

# Check versions
cat("ggtree version:", as.character(packageVersion("ggtree")), "\n")
cat("ggplot2 version:", as.character(packageVersion("ggplot2")), "\n")

# SECTION 2: GGTREE BASICS ===================================================

# ggtree extends ggplot2's "grammar of graphics" to phylogenetic trees
# This means you build plots layer by layer

# Create example mosquito tree
mosquito_tree <- read.tree(text = "((((Aedes_aegypti_USA:0.01,Aedes_aegypti_Brazil:0.01)100:0.02,
                                        (Aedes_albopictus_Japan:0.015,Aedes_albopictus_USA:0.015)98:0.015)95:0.05,
                                      ((Culex_pipiens_USA:0.02,Culex_pipiens_Europe:0.02)100:0.03,
                                       (Culex_quinquefasciatus_Brazil:0.025,Culex_quinquefasciatus_India:0.025)97:0.025)92:0.02)88:0.1,
                                     ((Anopheles_gambiae_Africa:0.03,Anopheles_gambiae_Lab:0.03)100:0.05,
                                      (Anopheles_stephensi_India:0.035,Anopheles_stephensi_Lab:0.035)99:0.045)85:0.05);")

# Basic ggtree plot (compare with ape's plot)
ggtree(mosquito_tree)

# Add tip labels
ggtree(mosquito_tree) +
  geom_tiplab()

# SECTION 3: TREE LAYOUTS ====================================================

# ggtree supports many layouts

# Rectangular/phylogram (default)
p1 <- ggtree(mosquito_tree) +
  geom_tiplab(size = 3) +
  ggtitle("Rectangular Layout")

# Slanted
p2 <- ggtree(mosquito_tree, layout = "slanted") +
  geom_tiplab(size = 3) +
  ggtitle("Slanted Layout")

# Circular
p3 <- ggtree(mosquito_tree, layout = "circular") +
  geom_tiplab(size = 3) +
  ggtitle("Circular Layout")

# Fan (radial)
p4 <- ggtree(mosquito_tree, layout = "fan", open.angle = 120) +
  geom_tiplab(size = 3) +
  ggtitle("Fan Layout")

# Display all layouts (requires patchwork or gridExtra)
# install.packages("patchwork")
library(patchwork)
(p1 | p2) / (p3 | p4)

# SECTION 4: STYLING THE TREE ================================================

# Basic tree with common customizations
basic_plot <- ggtree(mosquito_tree) +
  geom_tiplab(size = 3.5,          # Tip label size
              color = "darkblue",   # Label color
              fontface = "italic") + # Italic (for species names)
  geom_nodelab(size = 2.5,          # Node label size
               color = "red",        # Node label color
               hjust = 1.2,          # Horizontal adjustment
               vjust = -0.5) +       # Vertical adjustment
  theme_tree2() +                    # Theme with scale
  ggtitle("COI Phylogeny - Mosquito Species") +
  xlim(0, 0.25)                      # Extend x-axis for labels

print(basic_plot)

# SECTION 5: BRANCH STYLING ==================================================

# Color and style branches
styled_tree <- ggtree(mosquito_tree,
                      color = "steelblue",     # Branch color
                      size = 1.5,              # Branch width
                      linetype = "solid") +    # Line type
  geom_tiplab(size = 3) +
  theme_tree2()

print(styled_tree)

# Color branches by depth (distance from root)
depth_colored <- ggtree(mosquito_tree, aes(color = branch.length)) +
  geom_tiplab(size = 3, color = "black") +
  scale_color_gradient(low = "blue", high = "red") +
  theme_tree2() +
  theme(legend.position = "right")

print(depth_colored)

# SECTION 6: ADDING BOOTSTRAP/SUPPORT VALUES =================================

# Display node support values (bootstrap, posterior probability)
support_tree <- ggtree(mosquito_tree) +
  geom_tiplab(size = 3, fontface = "italic") +
  geom_nodelab(aes(label = label),         # Show bootstrap values
               size = 3,
               color = "red",
               hjust = 1.3,
               vjust = -0.5) +
  theme_tree2()

print(support_tree)

# Show only high support values (>90)
high_support <- ggtree(mosquito_tree) +
  geom_tiplab(size = 3) +
  geom_nodelab(aes(label = label,
                   subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),
               size = 3,
               color = "darkgreen") +
  theme_tree2()

print(high_support)

# Use symbols for support values
symbol_support <- ggtree(mosquito_tree) +
  geom_tiplab(size = 3) +
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),
                 color = "darkgreen",
                 size = 3) +
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) <= 90 & as.numeric(label) > 70),
                 color = "orange",
                 size = 3) +
  theme_tree2()

print(symbol_support)

# SECTION 7: ADDING METADATA =================================================

# Create metadata for our samples
metadata <- data.frame(
  label = mosquito_tree$tip.label,
  genus = c("Aedes", "Aedes", "Aedes", "Aedes",
            "Culex", "Culex", "Culex", "Culex",
            "Anopheles", "Anopheles", "Anopheles", "Anopheles"),
  location = c("USA", "Brazil", "Japan", "USA",
               "USA", "Europe", "Brazil", "India",
               "Africa", "Lab", "India", "Lab"),
  vector_status = c("Major", "Major", "Major", "Major",
                    "Major", "Major", "Major", "Major",
                    "Major", "Major", "Major", "Major")
)

# Join metadata with tree using %<+% operator
metadata_tree <- ggtree(mosquito_tree) %<+% metadata +
  geom_tiplab(aes(color = genus), size = 3, fontface = "italic") +
  theme_tree2() +
  scale_color_manual(values = c("Aedes" = "red",
                                 "Culex" = "blue",
                                 "Anopheles" = "darkgreen"))

print(metadata_tree)

# SECTION 8: TIP POINTS AND SYMBOLS ==========================================

# Add symbols at tips based on location
symbol_tree <- ggtree(mosquito_tree) %<+% metadata +
  geom_tiplab(size = 3, offset = 0.005) +
  geom_tippoint(aes(color = location, shape = genus), size = 3) +
  theme_tree2() +
  theme(legend.position = "right")

print(symbol_tree)

# SECTION 9: HIGHLIGHTING CLADES ============================================

# Highlight specific clades with background colors
clade_tree <- ggtree(mosquito_tree) +
  geom_tiplab(size = 3, offset = 0.005) +
  theme_tree2()

# Get node numbers for clades
# Use viewClade() to interactively find nodes, or:
# ggtree(mosquito_tree) + geom_text(aes(label=node), hjust=-.3)

# Highlight Aedes clade (adjust node number as needed)
# Note: you may need to adjust node numbers for your specific tree
clade_tree_highlighted <- clade_tree +
  geom_hilight(node = 15, fill = "pink", alpha = 0.3) +      # Aedes
  geom_hilight(node = 18, fill = "lightblue", alpha = 0.3) + # Culex
  geom_hilight(node = 21, fill = "lightgreen", alpha = 0.3)  # Anopheles

print(clade_tree_highlighted)

# SECTION 10: ADDING CLADE LABELS ===========================================

# Add text labels for clades
labeled_clades <- ggtree(mosquito_tree) +
  geom_tiplab(size = 2.5, offset = 0.005) +
  geom_cladelab(node = 15, label = "Aedes", offset = 0.02, fontsize = 4, color = "red") +
  geom_cladelab(node = 18, label = "Culex", offset = 0.02, fontsize = 4, color = "blue") +
  geom_cladelab(node = 21, label = "Anopheles", offset = 0.02, fontsize = 4, color = "darkgreen") +
  theme_tree2()

print(labeled_clades)

# SECTION 11: HEATMAPS AND ASSOCIATED DATA ===================================

# Create some associated data (e.g., gene expression, traits)
trait_data <- data.frame(
  label = mosquito_tree$tip.label,
  wing_length = c(2.5, 2.4, 2.2, 2.3, 3.1, 3.0, 3.2, 3.1, 2.8, 2.9, 2.7, 2.8),
  biting_rate = c(8, 7, 9, 8, 6, 5, 6, 7, 9, 8, 9, 9),
  insecticide_resistance = c(0.8, 0.7, 0.9, 0.85, 0.6, 0.4, 0.5, 0.65, 0.3, 0.2, 0.35, 0.25)
)

# Plot tree with heatmap
library(ggplot2)

# Create main tree plot
p <- ggtree(mosquito_tree) +
  geom_tiplab(size = 3, offset = 0.005) +
  theme_tree2()

# Add heatmap panel
gheatmap(p, trait_data[, c("wing_length", "biting_rate", "insecticide_resistance")],
         offset = 0.03,
         width = 0.3,
         colnames_angle = 90,
         colnames_offset_y = 0.5) +
  scale_fill_gradient(low = "white", high = "darkred")

# SECTION 12: COMPARING TWO TREES ===========================================

# Create a slightly different tree topology
alternative_tree <- read.tree(text = "((((Aedes_aegypti_USA:0.01,Aedes_aegypti_Brazil:0.01):0.02,
                                          (Aedes_albopictus_Japan:0.015,Aedes_albopictus_USA:0.015):0.015):0.04,
                                        (Culex_pipiens_USA:0.02,Culex_pipiens_Europe:0.02):0.06):0.1,
                                       ((Anopheles_gambiae_Africa:0.03,Anopheles_gambiae_Lab:0.03):0.05,
                                        (Anopheles_stephensi_India:0.035,Anopheles_stephensi_Lab:0.035):0.045):0.05);")

# Plot both trees side by side
p1 <- ggtree(mosquito_tree) + geom_tiplab(size = 2.5) + ggtitle("Tree 1")
p2 <- ggtree(alternative_tree) + geom_tiplab(size = 2.5) + ggtitle("Tree 2")

library(patchwork)
p1 | p2

# SECTION 13: SCALE BARS AND AXES ============================================

# Add scale bar
scale_tree <- ggtree(mosquito_tree) +
  geom_tiplab(size = 3) +
  geom_treescale(x = 0, y = 12, offset = 1) +  # Manual scale bar
  theme_tree()

print(scale_tree)

# Use theme_tree2() for automatic axis
axis_tree <- ggtree(mosquito_tree) +
  geom_tiplab(size = 3) +
  theme_tree2() +
  xlab("Substitutions per site")

print(axis_tree)

# SECTION 14: ZOOMING AND SUBSETTING =========================================

# Zoom into a specific clade
zoomed_tree <- ggtree(mosquito_tree) +
  geom_tiplab(size = 3) +
  xlim(0, 0.1) +  # Zoom into early part of tree
  theme_tree2()

print(zoomed_tree)

# SECTION 15: PUBLICATION-QUALITY FIGURE =====================================

# Create a complete publication-ready figure
publication_tree <- ggtree(mosquito_tree) %<+% metadata +
  # Tree with colored branches by genus
  geom_tree(aes(color = genus), size = 1) +
  # Tip labels in italic
  geom_tiplab(aes(color = genus),
              size = 3.5,
              fontface = "italic",
              offset = 0.003) +
  # Bootstrap support values
  geom_nodelab(aes(label = label,
                   subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),
               size = 3,
               color = "black",
               hjust = 1.3,
               vjust = -0.5) +
  # Color scheme
  scale_color_manual(
    values = c("Aedes" = "#E41A1C",
               "Culex" = "#377EB8",
               "Anopheles" = "#4DAF4A"),
    name = "Genus"
  ) +
  # Theme and labels
  theme_tree2(legend.position = c(0.15, 0.85)) +
  labs(title = "Phylogenetic Relationships of Mosquito Species",
       subtitle = "Based on COI Gene Sequences",
       x = "Substitutions per site",
       caption = "Bootstrap values >90 shown at nodes") +
  # Extend x-axis for labels
  xlim(0, 0.25)

print(publication_tree)

# Save the publication figure
ggsave("mosquito_phylogeny_publication.pdf",
       publication_tree,
       width = 10,
       height = 8,
       units = "in",
       dpi = 300)

ggsave("mosquito_phylogeny_publication.png",
       publication_tree,
       width = 10,
       height = 8,
       units = "in",
       dpi = 300)

# SECTION 16: CIRCULAR TREE WITH ANNOTATIONS =================================

# Create an attractive circular tree
circular_pub <- ggtree(mosquito_tree, layout = "circular") %<+% metadata +
  geom_tree(aes(color = genus), size = 1.2) +
  geom_tiplab(aes(color = genus),
              size = 3,
              fontface = "italic") +
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),
                 color = "black",
                 size = 2) +
  scale_color_manual(
    values = c("Aedes" = "#E41A1C",
               "Culex" = "#377EB8",
               "Anopheles" = "#4DAF4A")
  ) +
  theme(legend.position = "right") +
  labs(title = "Mosquito COI Phylogeny - Circular Layout")

print(circular_pub)

# SECTION 17: PRACTICAL FUNCTION - PLOT YOUR TREE ===========================

# Function to quickly create a nice tree plot
plot_phylogeny <- function(tree_file,
                          layout = "rectangular",
                          show_bootstrap = TRUE,
                          bootstrap_cutoff = 70,
                          tip_size = 3,
                          output_file = NULL) {
  # Read tree
  tree <- read.tree(tree_file)

  # Create base plot
  p <- ggtree(tree, layout = layout) +
    geom_tiplab(size = tip_size, fontface = "italic") +
    theme_tree2()

  # Add bootstrap values if requested
  if (show_bootstrap) {
    p <- p + geom_nodelab(
      aes(label = label,
          subset = !is.na(as.numeric(label)) & as.numeric(label) > bootstrap_cutoff),
      size = tip_size - 0.5,
      color = "red",
      hjust = 1.3,
      vjust = -0.5
    )
  }

  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 8, dpi = 300)
    cat("Tree saved to:", output_file, "\n")
  }

  return(p)
}

# Example usage:
# my_tree <- plot_phylogeny("../data/example_tree.newick",
#                           layout = "rectangular",
#                           show_bootstrap = TRUE,
#                           output_file = "my_phylogeny.pdf")

# EXERCISES TO TRY ==========================================================
#
# 1. Load your own tree and create:
#    - A rectangular tree with tip labels
#    - A circular tree
#    - A tree with bootstrap values highlighted
#
# 2. Create metadata for your samples including:
#    - Collection location
#    - Collection date
#    - Collector name
#    Then color tip labels by location
#
# 3. Make a publication-quality figure with:
#    - Colored clades
#    - Bootstrap support values
#    - Informative title and axis labels
#    - Save as PDF and PNG
#
# 4. Create a function that:
#    - Takes a tree file and metadata file
#    - Produces 4 different layouts
#    - Saves all as separate files
#
# 5. Experiment with different color schemes using:
#    - RColorBrewer palettes
#    - Viridis color scales
#    - Custom color vectors
#
################################################################################
# END OF SCRIPT
#
# Summary:
# - ggtree extends ggplot2 for phylogenetic trees
# - Build plots layer by layer with +
# - Many layouts: rectangular, circular, fan, etc.
# - Add metadata with %<+% operator
# - Style branches, tips, and nodes extensively
# - Create publication-quality figures with theme_tree2()
# - Save with ggsave() in multiple formats
#
# Next: Learn to calculate tree statistics and perform phylogenetic analyses!
################################################################################
