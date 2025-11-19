# Module 07: Phylogenetic Tree Construction

**Duration**: 3-4 hours
**Prerequisites**: Module 06

---

## Overview

Build phylogenetic trees using maximum likelihood and other methods. Understand evolutionary relationships among sequences.

---

## Tools

### IQ-TREE (Primary Method)
```bash
# Auto model selection + bootstrap
iqtree -s aligned.fasta -m MFP -bb 1000 -nt AUTO

# Output files:
# .treefile      - Best ML tree
# .iqtree        - Analysis log
# .log           - Run log
```

### RAxML-NG (Alternative)
```bash
raxml-ng --all --msa aligned.fasta --model GTR+G --bs-trees 100
```

### FastTree (Quick Trees)
```bash
# Very fast, less accurate
FastTree -nt -gtr < aligned.fasta > tree.newick
```

---

## IQ-TREE Workflow

### Step 1: Model Selection
```bash
# Find best substitution model
iqtree -s aligned.fasta -m MF
```

Common models for COI:
- GTR+G+I
- TN93+G
- HKY+G

### Step 2: Tree Building
```bash
# Build tree with bootstrap
iqtree -s aligned.fasta \
       -m GTR+G+I \
       -bb 1000 \
       -nt AUTO \
       -pre COI_tree

# -m: substitution model
# -bb: ultrafast bootstrap
# -nt: number of threads
# -pre: output prefix
```

### Step 3: Interpret Results

Bootstrap values:
- **>95%**: Strong support
- **70-95%**: Moderate support
- **<70%**: Weak support

---

## Example: Complete Pipeline

```bash
#!/bin/bash
# Build phylogenetic tree from COI sequences

# Input
ALIGNMENT="aligned_sequences.fasta"
PREFIX="COI_phylogeny"

echo "Starting phylogenetic analysis..."

# Run IQ-TREE
iqtree -s $ALIGNMENT \
       -m MFP \
       -bb 1000 \
       -nt AUTO \
       -pre $PREFIX

echo "Analysis complete!"
echo "Best tree: ${PREFIX}.treefile"
echo "View with: figtree ${PREFIX}.treefile"
```

---

## Visualizing Trees

### Using R (ggtree)

```r
library(ggtree)
library(ggplot2)

# Read tree
tree <- read.tree("COI_tree.treefile")

# Plot
p <- ggtree(tree) +
  geom_tiplab(size=3) +
  geom_nodepoint(aes(subset=as.numeric(label) > 70),
                 color="red", size=2) +
  theme_tree2()

ggsave("tree_figure.pdf", p, width=10, height=12)
```

### Using FigTree (GUI)
1. Open .treefile in FigTree
2. Display node labels (bootstrap values)
3. Format tree (colors, fonts)
4. Export as PDF or PNG

---

## Comparing Methods

```bash
# Run multiple methods
mafft --auto sequences.fasta > aligned.fasta

# IQ-TREE
iqtree -s aligned.fasta -m GTR+G -pre iqtree

# RAxML
raxml-ng --msa aligned.fasta --model GTR+G --prefix raxml

# FastTree
FastTree -nt -gtr < aligned.fasta > fasttree.newick

# Compare trees in R
```

---

## Tree Interpretation

### Reading a Phylogenetic Tree

```
                    ┌─ Species_A
              ┌─────┤
        ┌─────┤ 95  └─ Species_B
        │     └─────── Species_C
    ────┤
        │            ┌─ Species_D
        └────────────┤
                 75  └─ Species_E

Numbers = Bootstrap support (%)
Branch length = Evolutionary distance
```

### Key Questions

1. **Which species are most closely related?**
   - Species on adjacent branches

2. **How confident are we?**
   - Bootstrap values >70% are reliable

3. **What's the evolutionary distance?**
   - Branch lengths (substitutions per site)

---

## Common Issues

### Low Bootstrap Support
- Add more sequences
- Use more conserved region
- Try different model

### Long Branch Attraction
- Add intermediate taxa
- Use different tree method
- Check for saturation

### Star Topology (No Resolution)
- Sequences too similar (recent divergence)
- Sequences too divergent (saturation)
- Wrong gene/region

---

## Next Steps

Use tree for **Module 08: Species Identification** →
