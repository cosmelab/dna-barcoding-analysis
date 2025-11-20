# Getting Started with R for Phylogenetics

Welcome to Module 03! This guide will help you navigate the R basics tutorial for DNA barcoding analysis.

## Quick Start

### 1. Setup R and RStudio

**Install R**:
- Download from: https://cran.r-project.org/
- Choose your operating system
- Follow installation instructions

**Install RStudio** (recommended):
- Download from: https://posit.co/download/rstudio-desktop/
- RStudio provides a much better interface for R

**Install Required Packages**:
```r
# In R console, run:
install.packages(c("ape", "ggplot2", "dplyr", "patchwork"))

# For ggtree (from Bioconductor):
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")

# Optional but recommended:
install.packages(c("phangorn", "phytools", "RColorBrewer", "scales"))
```

### 2. Learning Path

Follow scripts in order:

#### Week 1: R Fundamentals
1. **01_r_syntax.R** (2 hours)
   - R basics: variables, vectors, data frames
   - Functions and loops
   - Working with DNA sequence data
   - **Exercise**: exercise_01_r_basics.R

#### Week 2: Phylogenetic Trees
2. **02_ape_basics.R** (2 hours)
   - Reading and writing trees
   - Tree structure and manipulation
   - Distance calculations
   - **Exercise**: exercise_02_tree_reading.R

3. **03_ggtree_intro.R** (2 hours)
   - Beautiful tree visualizations
   - Multiple layouts (rectangular, circular, fan)
   - Adding colors and annotations
   - **Exercise**: exercise_03_visualization.R

#### Week 3: Analysis
4. **04_tree_statistics.R** (2 hours)
   - Branch length analysis
   - Bootstrap support evaluation
   - Within/between group distances
   - **Exercise**: exercise_04_analysis.R

5. **05_publication_figures.R** (2 hours)
   - Creating publication-quality figures
   - Multi-panel layouts
   - Exporting in multiple formats
   - **Apply to your own data**

## How to Use This Module

### For Each Script:

1. **Read through the script first** - Don't just run everything
2. **Run code chunk by chunk** - Select lines and press Ctrl+Enter (Cmd+Enter on Mac)
3. **Modify examples** - Change values and see what happens
4. **Take notes** - Add your own comments
5. **Complete exercises** - Practice makes perfect!

### Working with Scripts in RStudio:

```r
# Open a script
# File > Open File > navigate to script

# Run a single line
# Place cursor on line, press Ctrl+Enter (Cmd+Enter on Mac)

# Run multiple lines
# Select lines, press Ctrl+Enter (Cmd+Enter)

# Run entire script
# Press Ctrl+Shift+Enter (Cmd+Shift+Enter)

# Set working directory to script location
# Session > Set Working Directory > To Source File Location
```

### Example Workflow:

```r
# 1. Set working directory
setwd("/path/to/03_r_basics")

# 2. Load libraries
library(ape)
library(ggtree)

# 3. Read example tree
tree <- read.tree("data/example_tree.newick")

# 4. Read metadata
metadata <- read.csv("data/tree_metadata.csv")

# 5. Create a basic plot
ggtree(tree) + geom_tiplab()

# 6. Add metadata and colors
ggtree(tree) %<+% metadata +
  geom_tiplab(aes(color = genus)) +
  theme_tree2()
```

## File Organization

```
03_r_basics/
├── GETTING_STARTED.md          # This file
├── README.md                    # Module overview
│
├── scripts/                     # Learning scripts
│   ├── 01_r_syntax.R           # Start here!
│   ├── 02_ape_basics.R
│   ├── 03_ggtree_intro.R
│   ├── 04_tree_statistics.R
│   └── 05_publication_figures.R
│
├── exercises/                   # Practice problems
│   ├── exercise_01_r_basics.R
│   ├── exercise_02_tree_reading.R
│   ├── exercise_03_visualization.R
│   └── exercise_04_analysis.R
│
├── solutions/                   # Check your work
│   ├── solution_01_r_basics.R
│   ├── solution_02_tree_reading.R
│   ├── solution_03_visualization.R
│   └── solution_04_analysis.R
│
└── data/                        # Example data
    ├── README.md                # Data documentation
    ├── example_tree.newick      # Main tree file
    ├── bootstrap_tree.treefile  # Tree with bootstrap
    ├── tree_metadata.csv        # Sample information
    └── trait_data.csv           # Phenotypic data
```

## Common Issues and Solutions

### Issue: Package won't install
**Solution**:
```r
# Try different mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("ape")

# Or use RStudio's package installer
# Tools > Install Packages
```

### Issue: Can't read file
**Solution**:
```r
# Check working directory
getwd()

# Set to correct location
setwd("/Users/yourusername/Projects/dna-barcoding-analysis/03_r_basics")

# Or use full file path
tree <- read.tree("/full/path/to/example_tree.newick")
```

### Issue: Tree doesn't plot correctly
**Solution**:
```r
# Make plot window larger in RStudio
# Drag the corner of the plots pane

# Or save to file and view
pdf("my_tree.pdf", width = 10, height = 12)
plot(tree)
dev.off()
```

### Issue: ggtree won't load
**Solution**:
```r
# ggtree is from Bioconductor, not CRAN
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")

# Load it
library(ggtree)
```

## Tips for Success

### 1. Work Through Examples
- Don't just read - actually run the code
- Modify values and see what changes
- Break things - that's how you learn!

### 2. Use R Help
```r
?function_name          # Help for specific function
??phylogenetic          # Search all help
example(plot)           # See examples
```

### 3. Comment Your Code
```r
# Good comments explain WHY, not what
tree <- read.tree("file.newick")  # Load COI tree from IQ-TREE analysis

# Not helpful:
tree <- read.tree("file.newick")  # Read tree
```

### 4. Build a Reference Library
Save useful code snippets:
```r
# My standard tree plot
my_tree_plot <- function(tree_file) {
  tree <- read.tree(tree_file)
  ggtree(tree) +
    geom_tiplab(size = 3, fontface = "italic") +
    theme_tree2() +
    labs(x = "Substitutions per site")
}
```

### 5. Practice with Your Own Data
- Once you complete the exercises, use your own trees
- Apply what you learned to real research questions
- Experiment with different visualizations

## Learning Objectives

By the end of this module, you should be able to:

- [ ] Write basic R code for data analysis
- [ ] Read and manipulate phylogenetic trees
- [ ] Calculate genetic distances and statistics
- [ ] Create beautiful tree visualizations
- [ ] Assess tree quality (bootstrap support)
- [ ] Export publication-ready figures
- [ ] Join metadata with trees
- [ ] Identify potential issues in phylogenies
- [ ] Interpret phylogenetic relationships

## Next Steps

After completing this module:

1. **Apply to your lab data** - Use your own COI sequences
2. **Read documentation** - Explore ape and ggtree manuals
3. **Join communities**:
   - [R phylogenetics Google Group](https://groups.google.com/g/r-sig-phylo)
   - [ggtree discussions](https://github.com/YuLab-SMU/ggtree/discussions)
4. **Explore advanced topics**:
   - Ancestral state reconstruction
   - Molecular clock analysis
   - Comparative phylogenetics
   - Phylogenetic networks

## Additional Resources

### Books
- [ggtree: Data Integration, Manipulation and Visualization of Phylogenetic Trees](https://yulab-smu.top/treedata-book/)
- [Analysis of Phylogenetics and Evolution with R](http://ape-package.ird.fr/APER.html)
- [R for Data Science](https://r4ds.had.co.nz/)

### Online Tutorials
- [R Phylogenetics Workshop](http://www.phytools.org/Cordoba2017/)
- [RStudio Education](https://education.rstudio.com/)
- [Datacamp R Courses](https://www.datacamp.com/courses/free-introduction-to-r)

### Package Documentation
- [ape manual](https://cran.r-project.org/web/packages/ape/ape.pdf)
- [ggtree vignettes](https://bioconductor.org/packages/release/bioc/html/ggtree.html)
- [ggplot2 documentation](https://ggplot2.tidyverse.org/)

### Video Tutorials
- Search YouTube for "R phylogenetics tutorial"
- [RStudio webinars](https://www.rstudio.com/resources/webinars/)

## Getting Help

1. **Read error messages carefully** - They often tell you what's wrong
2. **Check documentation** - Use `?function_name`
3. **Search online** - Google the error message
4. **Ask for help**:
   - Post on Stack Overflow with tag [r] [phylogenetics]
   - Ask in course discussion forum
   - Contact instructor during office hours

## Good Luck!

Remember: Everyone struggles with R at first. It gets easier with practice. Don't be discouraged by errors - they're part of the learning process!

**Start with 01_r_syntax.R and work your way through. You've got this!**

---

**Questions?** Check the main README.md or contact your instructor.

**Last Updated**: November 2025
