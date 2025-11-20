#!/usr/bin/env Rscript
################################################################################
# Exercise 03: Phylogenetic Tree Visualization
#
# Practice creating beautiful tree visualizations with ggtree
# Learn to customize plots for different purposes
#
# Author: Luciano Cosme
# Course: ENTM201L - Molecular Biology Laboratory
# UC Riverside, Fall 2025
################################################################################

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)

# EXERCISE 1: BASIC TREE PLOTTING ===========================================

# 1.1 Read the example tree
# tree <- read.tree("../data/example_tree.newick")

# YOUR CODE HERE:


# 1.2 Create a basic ggtree plot

# YOUR CODE HERE:


# 1.3 Add tip labels to the plot

# YOUR CODE HERE:


# 1.4 Add a theme with axis (theme_tree2)

# YOUR CODE HERE:


# 1.5 Add a title and axis label

# YOUR CODE HERE:


# EXERCISE 2: TREE LAYOUTS ==================================================

# 2.1 Create a rectangular tree layout

# YOUR CODE HERE:


# 2.2 Create a circular tree layout

# YOUR CODE HERE:


# 2.3 Create a fan layout with 120 degree opening

# YOUR CODE HERE:


# 2.4 Create a slanted tree layout

# YOUR CODE HERE:


# 2.5 Which layout is best for your tree? Why?

# YOUR ANSWER:


# EXERCISE 3: STYLING TIPS ==================================================

# 3.1 Change tip label size to 4

# YOUR CODE HERE:


# 3.2 Make tip labels italic (for species names)

# YOUR CODE HERE:


# 3.3 Color tip labels by genus
# Hint: Use metadata or extract genus from tip labels

# YOUR CODE HERE:


# 3.4 Add tip points colored by location

# YOUR CODE HERE:


# 3.5 Offset tip labels so they don't overlap the points

# YOUR CODE HERE:


# EXERCISE 4: BOOTSTRAP VALUES ==============================================

# If your tree has bootstrap values:

# 4.1 Display all bootstrap values at nodes

# YOUR CODE HERE:


# 4.2 Display only bootstrap values >70

# YOUR CODE HERE:


# 4.3 Use colored points for bootstrap support:
# - Green for >90
# - Orange for 70-90
# - Red for <70

# YOUR CODE HERE:


# 4.4 Create a legend explaining the bootstrap colors

# YOUR CODE HERE:


# EXERCISE 5: WORKING WITH METADATA =========================================

# 5.1 Create a metadata data frame for your samples
# Include: genus, species, location, year

# YOUR CODE HERE:


# 5.2 Join the metadata with your tree using %<+%

# YOUR CODE HERE:


# 5.3 Color tip labels by genus using the metadata

# YOUR CODE HERE:


# 5.4 Add symbols at tips showing location

# YOUR CODE HERE:


# 5.5 Create a tree where:
# - Branch colors represent genus
# - Tip labels show species only
# - Tip points show collection year

# YOUR CODE HERE:


# EXERCISE 6: HIGHLIGHTING CLADES ===========================================

# 6.1 Find the node number for the Aedes clade
# Hint: Use ggtree(tree) + geom_text(aes(label=node))

# YOUR CODE HERE:


# 6.2 Highlight the Aedes clade with a colored background

# YOUR CODE HERE:


# 6.3 Highlight all major clades (Aedes, Culex, Anopheles)
# Use different colors

# YOUR CODE HERE:


# 6.4 Add text labels for each clade

# YOUR CODE HERE:


# EXERCISE 7: BRANCH STYLING ================================================

# 7.1 Color all branches by their length
# Hint: Use aes(color = branch.length)

# YOUR CODE HERE:


# 7.2 Make branches thicker (size = 1.5)

# YOUR CODE HERE:


# 7.3 Color branches by genus

# YOUR CODE HERE:


# 7.4 Create a tree where branch color represents:
# - Blue for short branches (<0.01)
# - Red for long branches (>0.05)

# YOUR CODE HERE:


# EXERCISE 8: COLOR SCHEMES =================================================

# 8.1 Create a tree with a colorblind-friendly palette
# Use: c("#E69F00", "#0072B2", "#009E73") for three genera

# YOUR CODE HERE:


# 8.2 Use a viridis color scale for continuous data

# YOUR CODE HERE:


# 8.3 Use an RColorBrewer palette (e.g., "Set1" or "Dark2")

# YOUR CODE HERE:


# 8.4 Create a custom color scheme matching your institution colors

# YOUR CODE HERE:


# EXERCISE 9: COMBINING ANNOTATIONS =========================================

# 9.1 Create a publication-quality tree with:
# - Colored branches by genus
# - Italic tip labels
# - Bootstrap values >90 shown
# - Scale bar
# - Informative title
# - Professional theme

# YOUR CODE HERE:


# 9.2 Add a legend in an appropriate position

# YOUR CODE HERE:


# 9.3 Adjust x-axis limits to make room for labels

# YOUR CODE HERE:


# EXERCISE 10: MULTI-PANEL PLOTS ============================================

# 10.1 Create two versions of your tree:
# - One rectangular
# - One circular
# Display side by side using patchwork

library(patchwork)

# YOUR CODE HERE:


# 10.2 Create a three-panel figure:
# - Tree with metadata
# - Distance distribution plot
# - Bootstrap support histogram

# YOUR CODE HERE:


# EXERCISE 11: HEATMAPS =====================================================

# 11.1 Create trait data for your samples
# Example: wing_length, biting_rate, resistance_level

# YOUR CODE HERE:


# 11.2 Add a heatmap next to the tree using gheatmap()

# YOUR CODE HERE:


# 11.3 Customize the heatmap color scale

# YOUR CODE HERE:


# EXERCISE 12: SAVING FIGURES ===============================================

# 12.1 Save your best tree as a PDF (7" x 8", 300 DPI)

# YOUR CODE HERE:


# 12.2 Save as a high-resolution PNG

# YOUR CODE HERE:


# 12.3 Save as SVG (editable in Illustrator)

# YOUR CODE HERE:


# 12.4 Create a single-column width version (3.5" wide)
# for journal submission

# YOUR CODE HERE:


# EXERCISE 13: CIRCULAR TREE WITH ANNOTATIONS ===============================

# 13.1 Create a circular tree

# YOUR CODE HERE:


# 13.2 Add colored bars around the outside showing genus

# YOUR CODE HERE:


# 13.3 Add a second ring showing location

# YOUR CODE HERE:


# 13.4 Make the circular tree visually appealing with:
# - Good color choices
# - Clear labels
# - Legend
# - Title

# YOUR CODE HERE:


# EXERCISE 14: COMPARISON CHALLENGE =========================================

# 14.1 Create 4 different visualizations of the same tree:
# - Basic rectangular with tip labels
# - Colored by genus with clades highlighted
# - Circular with metadata
# - With heatmap of trait data

# YOUR CODE HERE:


# 14.2 Combine all 4 into one figure using patchwork

# YOUR CODE HERE:


# 14.3 Which visualization best shows:
# - Overall relationships?
# - Genus-level patterns?
# - Sample metadata?
# - Quantitative traits?

# YOUR ANSWERS:


# BONUS CHALLENGE ==========================================================

# Create a function that generates a complete publication figure:
# - Takes tree file, metadata file, trait data file
# - Produces a multi-panel figure showing:
#   * Main phylogeny with colors and annotations
#   * Bootstrap support summary
#   * Trait heatmap
#   * Distance distribution
# - Saves in multiple formats
# - Returns the plot object

create_complete_figure <- function(tree_file,
                                   metadata_file,
                                   trait_file,
                                   output_prefix) {
  # YOUR CODE HERE

}

# Test your function


# ADVANCED CHALLENGE: INTERACTIVE TREES ====================================

# Research how to create an interactive tree with:
# - ggtree + plotly
# - Or use the 'treeio' and 'tidytree' packages
# - Allow users to click on tips to see information

# YOUR CODE HERE:


################################################################################
# REFLECTION QUESTIONS
#
# After completing these exercises, answer these questions:
#
# 1. What are the advantages of ggtree over base R plotting?
#
# 2. When would you use a circular vs rectangular layout?
#
# 3. How do you choose appropriate colors for publication figures?
#
# 4. What information should always be included in a phylogenetic figure?
#
# 5. How do you balance aesthetic appeal with scientific clarity?
#
################################################################################

# Check your answers with:
# source("../solutions/solution_03_visualization.R")
