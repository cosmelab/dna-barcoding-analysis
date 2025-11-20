# Module 03: R Basics for Phylogenetics

**Duration**: 3-4 hours
**Difficulty**: Beginner
**Prerequisites**: Module 01 (Linux Basics)

---

## Overview

Learn R programming for phylogenetic analysis and visualization. R is the gold standard for creating publication-quality phylogenetic trees and performing comparative analyses.

**What You'll Learn**:
- R syntax and data structures
- Read and manipulate phylogenetic trees
- Beautiful tree visualizations with ggtree
- Statistical analysis of phylogenies
- Export publication-ready figures

---

## Module Contents

```
03_r_basics/
├── README.md                    # This file
├── scripts/                     # R scripts
│   ├── 01_r_syntax.R           # Basics
│   ├── 02_ape_basics.R         # ape package
│   ├── 03_ggtree_intro.R       # ggtree visualization
│   ├── 04_tree_statistics.R    # Tree analysis
│   └── 05_publication_figures.R # Export figures
├── exercises/                   # Practice
├── solutions/                   # Answers
└── data/                        # Example trees
    ├── example_tree.newick
    └── tree_metadata.csv
```

---

## Key R Packages

### ape - Analyses of Phylogenetics and Evolution
```r
library(ape)

# Read tree
tree <- read.tree("tree.newick")

# Basic info
Ntip(tree)              # Number of tips
tree$tip.label          # Tip labels

# Plot basic tree
plot(tree)
```

### ggtree - Grammar of Graphics for Trees
```r
library(ggtree)

# Beautiful tree
ggtree(tree) +
  geom_tiplab() +
  theme_tree2()
```

### tidyverse - Data Manipulation
```r
library(tidyverse)

# Read metadata
metadata <- read_csv("metadata.csv")

# Join with tree
ggtree(tree) %<+% metadata +
  geom_tiplab(aes(color=species))
```

---

## Essential Tree Operations

### Read/Write Trees
```r
# Read Newick format
tree <- read.tree("tree.newick")

# Read Nexus format
tree <- read.nexus("tree.nex")

# Write tree
write.tree(tree, "output.newick")
```

### Basic Tree Manipulation
```r
# Root tree
rooted_tree <- root(tree, outgroup="Outgroup_species")

# Drop tips
pruned_tree <- drop.tip(tree, c("tip1", "tip2"))

# Get subtree
subtree <- extract.clade(tree, node=50)
```

### Tree Visualization
```r
# Basic plot
ggtree(tree) +
  geom_tiplab(size=3) +
  geom_nodelab(size=2, color="red") +
  theme_tree2()

# Circular tree
ggtree(tree, layout="circular") +
  geom_tiplab(size=3)

# With bootstrap values
ggtree(tree) +
  geom_tiplab() +
  geom_text(aes(label=label), hjust=-.3)
```

---

## Example: Complete Tree Visualization

```r
#!/usr/bin/env Rscript

library(ggtree)
library(ggplot2)
library(treeio)

# Read tree with bootstrap values
tree <- read.tree("COI_tree.treefile")

# Create beautiful plot
p <- ggtree(tree, layout="rectangular") +
  geom_tiplab(size=3, align=TRUE) +
  geom_nodepoint(aes(subset=as.numeric(label) > 70),
                 color="red", size=2) +
  theme_tree2() +
  ggtitle("COI Phylogenetic Tree") +
  xlim(0, max(tree$edge.length) * 1.3)

# Save publication-quality figure
ggsave("tree_figure.pdf", p, width=10, height=12, units="in", dpi=300)
ggsave("tree_figure.png", p, width=10, height=12, units="in", dpi=300)
```

---

## Tips

1. **RStudio Recommended** - Best IDE for R
2. **Start Simple** - Basic plots first, then add complexity
3. **Use %>% Pipes** - Makes code readable
4. **Save Plots as PDF** - Vector graphics scale perfectly
5. **Check Package Versions** - `sessionInfo()`

---

## Common Tasks

### Add Metadata to Tree
```r
# Metadata file (CSV)
# tip_label,species,location
# Sample1,Aedes_aegypti,California
# Sample2,Aedes_albopictus,Florida

metadata <- read_csv("metadata.csv")

ggtree(tree) %<+% metadata +
  geom_tiplab(aes(color=species)) +
  theme_tree2()
```

### Highlight Clades
```r
ggtree(tree) +
  geom_tiplab() +
  geom_hilight(node=75, fill="steelblue", alpha=0.3) +
  geom_hilight(node=120, fill="darkgreen", alpha=0.3)
```

### Export Multiple Formats
```r
# PDF (vector)
ggsave("tree.pdf", width=10, height=8)

# PNG (raster, high resolution)
ggsave("tree.png", width=10, height=8, dpi=300)

# SVG (editable in Illustrator)
ggsave("tree.svg", width=10, height=8)
```

---

## Resources

- [ggtree Book](https://yulab-smu.top/treedata-book/)
- [ape Documentation](http://ape-package.ird.fr/)
- [R for Data Science](https://r4ds.had.co.nz/)

---

**Next**: Apply these skills in Module 07 (Phylogeny) →
