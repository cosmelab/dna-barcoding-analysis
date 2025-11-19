# R Phylogenetics Cheat Sheet

Quick reference for common R phylogenetic operations.

## Setup

```r
# Load packages
library(ape)        # Phylogenetic analysis
library(ggtree)     # Tree visualization
library(ggplot2)    # Graphics
library(dplyr)      # Data manipulation

# Set working directory
setwd("path/to/03_r_basics")
```

## Reading & Writing Trees

```r
# Read tree
tree <- read.tree("data/example_tree.newick")
tree <- read.nexus("tree.nex")

# Write tree
write.tree(tree, "output.newick")
write.nexus(tree, "output.nex")

# Read multiple trees
trees <- read.tree("multiple_trees.newick")  # Returns list
```

## Basic Tree Information

```r
# Tree properties
Ntip(tree)              # Number of tips
tree$Nnode              # Number of internal nodes
tree$tip.label          # Tip names
tree$edge.length        # Branch lengths
is.rooted(tree)         # Is tree rooted?
is.binary(tree)         # Is fully resolved?
is.ultrametric(tree)    # Are tips equidistant?

# Tree structure
str(tree)               # Detailed structure
print(tree)             # Summary
```

## Tree Manipulation

```r
# Rooting
tree_rooted <- root(tree, outgroup = "Anopheles_gambiae")
tree_midpoint <- midpoint.root(tree)
tree_unrooted <- unroot(tree)

# Pruning
tree_subset <- keep.tip(tree, c("tip1", "tip2", "tip3"))
tree_dropped <- drop.tip(tree, "tip_to_remove")

# Ladderizing (organize for clarity)
tree_ladder <- ladderize(tree, right = TRUE)

# Extract clade
subtree <- extract.clade(tree, node = 20)
```

## Distances

```r
# Pairwise tip-to-tip distances
dist_matrix <- cophenetic.phylo(tree)

# Specific distance
dist_matrix["tip1", "tip2"]

# Distance statistics
mean(dist_matrix[lower.tri(dist_matrix)])  # Mean pairwise
max(dist_matrix)                            # Maximum distance
```

## Tree Statistics

```r
# Branch lengths
sum(tree$edge.length)           # Total tree length
mean(tree$edge.length)          # Mean branch length
range(tree$edge.length)         # Min and max

# Node depths
node_depths <- node.depth.edgelength(tree)
root_to_tip <- node_depths[1:Ntip(tree)]

# Bootstrap values (if present)
boot_vals <- as.numeric(tree$node.label)
boot_vals <- boot_vals[!is.na(boot_vals)]
mean(boot_vals)                 # Mean support
sum(boot_vals > 90)             # Count high support
```

## Basic Plotting with ape

```r
# Simple plot
plot(tree)

# With options
plot(tree,
     cex = 0.8,                  # Label size
     font = 2,                   # Bold
     tip.color = "blue",         # Tip label color
     edge.color = "gray",        # Branch color
     type = "phylogram")         # Tree type

# Add elements
add.scale.bar()                  # Scale bar
nodelabels()                     # Node numbers
tiplabels()                      # Tip numbers
edgelabels()                     # Edge labels

# Tree types
plot(tree, type = "phylogram")   # Default
plot(tree, type = "cladogram")   # Ignore branch lengths
plot(tree, type = "fan")         # Radial
plot(tree, type = "unrooted")    # Unrooted network
```

## ggtree Visualization

```r
# Basic ggtree
ggtree(tree)

# Add components
ggtree(tree) +
  geom_tiplab() +                # Tip labels
  geom_nodelab() +               # Node labels
  theme_tree2()                  # Theme with axis

# Layouts
ggtree(tree, layout = "rectangular")
ggtree(tree, layout = "circular")
ggtree(tree, layout = "fan", open.angle = 120)
ggtree(tree, layout = "slanted")

# Styling
ggtree(tree) +
  geom_tiplab(
    size = 3,                    # Label size
    fontface = "italic",         # Italic
    color = "blue",              # Color
    offset = 0.01                # Distance from tips
  ) +
  geom_tree(
    size = 1,                    # Branch width
    color = "gray50"             # Branch color
  )
```

## Adding Metadata

```r
# Read metadata
metadata <- read.csv("data/tree_metadata.csv")

# Join with tree
ggtree(tree) %<+% metadata +
  geom_tiplab(aes(color = genus))  # Color by genus

# Multiple aesthetics
ggtree(tree) %<+% metadata +
  geom_tiplab(aes(color = genus), size = 3) +
  geom_tippoint(aes(shape = location), size = 3)
```

## Color Schemes

```r
# Manual colors
scale_color_manual(
  values = c("Aedes" = "#E69F00",
             "Culex" = "#0072B2",
             "Anopheles" = "#009E73")
)

# Gradient
scale_color_gradient(low = "blue", high = "red")

# RColorBrewer
scale_color_brewer(palette = "Set1")

# Viridis (colorblind-friendly)
scale_color_viridis_d()
```

## Bootstrap Values

```r
# Show all bootstrap values
ggtree(tree) +
  geom_nodelab(aes(label = label), size = 3)

# Show only high support (>70)
ggtree(tree) +
  geom_nodelab(
    aes(label = label,
        subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
    size = 3
  )

# Color-coded points
ggtree(tree) +
  geom_nodepoint(
    aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),
    color = "darkgreen",
    size = 3
  )
```

## Highlighting Clades

```r
# Find node numbers first
ggtree(tree) + geom_text(aes(label = node), hjust = -0.3)

# Highlight clade
ggtree(tree) +
  geom_hilight(node = 20, fill = "steelblue", alpha = 0.3)

# Add clade label
ggtree(tree) +
  geom_cladelab(node = 20, label = "Aedes", offset = 0.02)

# Multiple clades
ggtree(tree) +
  geom_hilight(node = 20, fill = "red", alpha = 0.2) +
  geom_hilight(node = 25, fill = "blue", alpha = 0.2)
```

## Heatmaps

```r
# Create tree plot
p <- ggtree(tree) + geom_tiplab()

# Add heatmap
gheatmap(p, trait_data,
         offset = 0.03,          # Space from tree
         width = 0.3,            # Heatmap width
         colnames_angle = 90)    # Column name angle
```

## Saving Figures

```r
# Save current plot
ggsave("tree.pdf", width = 8, height = 10, dpi = 300)
ggsave("tree.png", width = 8, height = 10, dpi = 300)
ggsave("tree.svg", width = 8, height = 10)

# Save specific plot
ggsave("tree.pdf", plot = my_plot, width = 8, height = 10)

# Journal formats
ggsave("tree.pdf", width = 3.5, height = 5)  # Single column
ggsave("tree.pdf", width = 7, height = 8)    # Double column
```

## Multi-Panel Figures

```r
library(patchwork)

# Side by side
p1 | p2

# Stacked
p1 / p2

# Complex layouts
(p1 | p2) / p3

# Add annotations
(p1 | p2) +
  plot_annotation(
    title = "Figure 1",
    subtitle = "Phylogenetic analysis"
  )
```

## Common Analyses

```r
# Within vs between group distances
genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
dist_mat <- cophenetic.phylo(tree)

within_dist <- c()
between_dist <- c()
for (i in 1:(Ntip(tree)-1)) {
  for (j in (i+1):Ntip(tree)) {
    if (genera[i] == genera[j]) {
      within_dist <- c(within_dist, dist_mat[i,j])
    } else {
      between_dist <- c(between_dist, dist_mat[i,j])
    }
  }
}

# Compare
t.test(within_dist, between_dist)
```

## Tree Comparison

```r
# Robinson-Foulds distance (topological difference)
rf_dist <- treedist(tree1, tree2, method = "RF")

# Branch score (considers branch lengths)
branch_score <- treedist(tree1, tree2, method = "score")
```

## Useful Functions

```r
# Data manipulation
head(data)              # First 6 rows
tail(data)              # Last 6 rows
str(data)               # Structure
summary(data)           # Summary statistics
dim(data)               # Dimensions
names(data)             # Column names

# Subsetting
data[1:5, ]            # First 5 rows
data[, c("col1", "col2")]  # Select columns
data[data$genus == "Aedes", ]  # Filter rows

# Apply functions
sapply(list, function)  # Returns vector
lapply(list, function)  # Returns list
apply(matrix, 1, mean)  # Row means
apply(matrix, 2, mean)  # Column means
```

## Keyboard Shortcuts (RStudio)

```
Ctrl/Cmd + Enter      Run current line/selection
Ctrl/Cmd + Shift + Enter  Run entire script
Ctrl/Cmd + 1          Move cursor to editor
Ctrl/Cmd + 2          Move cursor to console
Tab                   Auto-complete
Ctrl/Cmd + Shift + C  Comment/uncomment
Alt + -               Insert <-
Ctrl/Cmd + L          Clear console
```

## Getting Help

```r
?function_name          # Function help
??search_term           # Search help
example(function_name)  # Run examples
help.search("topic")    # Search by topic
apropos("pattern")      # Find functions matching pattern
```

## Common Errors

**Error: object not found**
- Check spelling and capitalization
- Make sure you loaded the data
- Check working directory

**Error: could not find function**
- Load required package with `library()`
- Check package is installed

**Error in read.tree: unexpected character**
- Check file format
- Ensure file path is correct
- Look for special characters in tree file

**Plots look weird**
- Resize plot window
- Adjust xlim/ylim
- Change figure dimensions in ggsave()

## Quick Start Template

```r
#!/usr/bin/env Rscript

# Load packages
library(ape)
library(ggtree)
library(ggplot2)

# Read data
tree <- read.tree("data/example_tree.newick")
metadata <- read.csv("data/tree_metadata.csv")

# Basic analysis
dist_mat <- cophenetic.phylo(tree)
cat("Mean pairwise distance:", mean(dist_mat[lower.tri(dist_mat)]), "\n")

# Visualization
p <- ggtree(tree) %<+% metadata +
  geom_tiplab(aes(color = genus), size = 3, fontface = "italic") +
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) &
                      as.numeric(label) > 90),
                 color = "black", size = 2) +
  theme_tree2() +
  labs(title = "Phylogenetic Tree",
       x = "Substitutions per site")

# Save
ggsave("tree.pdf", p, width = 8, height = 10, dpi = 300)

print(p)
```

---

**Print this out and keep it handy while coding!**

**Last Updated**: November 2025
