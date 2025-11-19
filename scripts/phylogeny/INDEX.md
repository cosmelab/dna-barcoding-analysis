# Phylogenetic Analysis Scripts - Complete Index

## Directory Contents

```
scripts/phylogeny/
├── run_iqtree.sh              [EXECUTABLE] IQ-TREE wrapper (345 lines)
├── visualize_tree.py          [EXECUTABLE] Tree visualization (468 lines)
├── root_tree.py               [EXECUTABLE] Tree rooting methods (525 lines)
├── calculate_distances.py     [EXECUTABLE] Genetic distance calculations (643 lines)
├── README.md                  Comprehensive guide & tutorials (507 lines)
├── QUICK_REFERENCE.md         Command cheatsheet & FAQ
└── INDEX.md                   This file
```

## Choose Your Starting Point

### I'm a Complete Beginner
Start here: **QUICK_REFERENCE.md**
- See one-line summaries
- Find command examples
- Check FAQ section

### I Want to Learn Phylogenetics
Start here: **README.md**
- Read "Phylogenetic Concepts Explained"
- Follow "Quick Start" section
- Try example commands

### I Want to Use These Scripts Now
Start here: **QUICK_REFERENCE.md**
- Find your task in the quick commands
- Copy the example
- Run it!

### I Need Help with a Specific Script
- `run_iqtree.sh`: Read "IQ-TREE Wrapper" section in README.md
- `visualize_tree.py`: Read "Tree Visualization" section in README.md
- `root_tree.py`: Read "Rooting Methods" section in README.md
- `calculate_distances.py`: Read "Genetic Distances" section in README.md

## Script Quick Links

### 1. run_iqtree.sh
**Purpose**: Build phylogenetic trees using maximum likelihood

**Main Concept**: Finds the best-fitting evolutionary tree for your sequences

**Quick Command**:
```bash
./run_iqtree.sh your_alignment.fasta
```

**What You Need**:
- Aligned FASTA file (use MAFFT first if needed)

**What You Get**:
- `.treefile` - Your phylogenetic tree
- `.iqtree` - Detailed analysis report
- `.log` - Run log

**Learn More**: `./run_iqtree.sh -h` or see README.md

---

### 2. visualize_tree.py
**Purpose**: Create beautiful, publication-quality tree figures

**Main Concept**: Draws tree structure with bootstrap values and evolutionary distances

**Quick Command**:
```bash
python visualize_tree.py your_tree.treefile -o tree.png
```

**What You Need**:
- Tree file (from run_iqtree.sh or other phylogenetic tool)

**What You Get**:
- PNG/PDF/JPG image of your tree
- Bootstrap values shown on branches
- Evolutionary distances as branch lengths

**Learn More**: `python visualize_tree.py -h` or see README.md

---

### 3. root_tree.py
**Purpose**: Reroot trees to show evolutionary history

**Main Concept**: Adds directionality to show which species are "ancestral"

**Quick Commands**:
```bash
# Midpoint rooting (use if no outgroup)
python root_tree.py tree.treefile midpoint -o rooted.treefile

# Outgroup rooting (use if you know a distantly related species)
python root_tree.py tree.treefile outgroup --outgroup "Distant_species" -o rooted.treefile
```

**What You Need**:
- Tree file
- Optionally: name of outgroup species

**What You Get**:
- Rooted tree with evolutionary directionality
- Ladderized for better visualization

**Learn More**: `python root_tree.py -h` or see README.md

---

### 4. calculate_distances.py
**Purpose**: Calculate how different sequences are (for species identification)

**Main Concept**: Genetic distance thresholds define species boundaries in DNA barcoding

**Quick Commands**:
```bash
# From alignment
python calculate_distances.py -a alignment.fasta --threshold-analysis

# From tree
python calculate_distances.py -t tree.treefile
```

**What You Need**:
- Either an alignment file OR a tree file

**What You Get**:
- Distance matrix (all pairwise distances)
- Species identification guidance
- CSV export for spreadsheet analysis

**Learn More**: `python calculate_distances.py -h` or see README.md

## Common Workflows

### Workflow 1: Quick Tree Visualization (10 minutes)
```bash
# Start with aligned sequences
./run_iqtree.sh my_alignment.fasta
python visualize_tree.py my_alignment_iqtree.treefile -o tree.png
```

### Workflow 2: Complete Phylogenetic Analysis (30 minutes)
```bash
# Build tree
./run_iqtree.sh alignment.fasta -bb 2000

# Visualize
python visualize_tree.py alignment_iqtree.treefile -o tree.png

# Root with outgroup
python root_tree.py alignment_iqtree.treefile outgroup \
    --outgroup "Outgroup_sp" -o rooted.treefile

# Visualize rooted tree
python visualize_tree.py rooted.treefile -o rooted_tree.png

# Calculate distances for identification
python calculate_distances.py -a alignment.fasta \
    --threshold-analysis -o distances.csv
```

### Workflow 3: Species Identification (5 minutes)
```bash
# Quickly identify unknown species using reference tree
python calculate_distances.py -a alignment.fasta

# Find the closest match distance
# If < 0.03 (3%) = same species
# If > 0.03 (3%) = different species
```

### Workflow 4: Learning Phylogenetics (1 hour)
```bash
# 1. Read concepts
python visualize_tree.py --concepts
python root_tree.py --concepts
python calculate_distances.py --concepts

# 2. Run with example data
./run_iqtree.sh example.fasta

# 3. Experiment with parameters
python visualize_tree.py tree.treefile -w 14 -h 10

# 4. Try different rooting methods
python root_tree.py tree.treefile midpoint -o v1.treefile
python root_tree.py tree.treefile outgroup -g "sp1" -o v2.treefile
```

## Getting Help

### Quick Help
- `run_iqtree.sh -h`
- `python visualize_tree.py -h`
- `python root_tree.py -h`
- `python calculate_distances.py -h`

### Learn Concepts
- `python visualize_tree.py --concepts`
- `python root_tree.py --concepts`
- `python calculate_distances.py --concepts`

### Detailed Guides
- See **README.md** for comprehensive information
- See **QUICK_REFERENCE.md** for command cheatsheet

## Common Tasks

| Task | Command |
|------|---------|
| Build a tree | `./run_iqtree.sh alignment.fasta` |
| View the tree | `python visualize_tree.py tree.treefile` |
| Save tree as image | `python visualize_tree.py tree.treefile -o tree.png` |
| Root with outgroup | `python root_tree.py tree.treefile outgroup --outgroup "species"` |
| Calculate distances | `python calculate_distances.py -a alignment.fasta` |
| Find species threshold | `python calculate_distances.py -a alignment.fasta --threshold-analysis` |
| Save distances to CSV | `python calculate_distances.py -a alignment.fasta -o dist.csv` |
| List all species in tree | `python root_tree.py tree.treefile --list-species` |
| High-quality PDF figure | `python visualize_tree.py tree.treefile -o figure.pdf --dpi 300` |

## What Each Script Teaches

### run_iqtree.sh
Learn about:
- Maximum likelihood phylogenetic inference
- Model selection (how to choose best evolutionary model)
- Bootstrap resampling (testing confidence)
- Tree topology and branch lengths
- Substitution models (GTR, HKY, TN93, etc.)

### visualize_tree.py
Learn about:
- Reading phylogenetic tree structure
- Bootstrap support interpretation
- Branch length scaling
- Scientific visualization with matplotlib
- Publication-quality figure creation

### root_tree.py
Learn about:
- Rooted vs unrooted trees
- Evolutionary directionality
- Outgroup concept
- Midpoint rooting assumptions
- Tree ladderization (visual organization)

### calculate_distances.py
Learn about:
- Genetic distance metrics
- Distance correction methods
- Evolutionary rate models
- DNA barcoding gaps
- Species identification thresholds
- Intraspecific vs interspecific variation

## Educational Philosophy

These scripts are intentionally:
- **Well-commented**: Every section explains what and why
- **Over-documented**: Better too much info than too little
- **Modular**: Each script teaches one concept
- **Practical**: Real-world DNA barcoding applications
- **Progressive**: Start simple, add complexity
- **Self-contained**: Work independently

Remove comments in production code, but keep them for teaching!

## For Instructors

### Creating Assignments

**Easy (20 min)**: "Build a tree and visualize it"
```bash
./run_iqtree.sh provided_alignment.fasta
python visualize_tree.py *_iqtree.treefile -o mytree.png
```

**Medium (1 hour)**: "Complete phylogenetic analysis"
- Build tree with IQ-TREE
- Visualize with bootstrap values
- Try different rooting methods
- Calculate distances
- Write interpretation

**Hard (2+ hours)**: "Reproduce published analysis"
- Download data from BOLD
- Align with MAFFT
- Build tree with these scripts
- Recreate figure from paper
- Compare your results

### Assessment Rubrics

**Tree Quality**:
- Bootstrap values > 70%? (8/10 points)
- Appropriate model selected? (2/10 points)

**Figure Quality**:
- Clear and labeled? (5/10 points)
- Publication-ready resolution? (5/10 points)

**Species Identification**:
- Correct species boundaries? (6/10 points)
- Correct distance interpretation? (4/10 points)

## Requirements

### Software
- Python 3.6+ with:
  - Biopython
  - matplotlib
  - numpy
- IQ-TREE 2.0+
- MAFFT (for alignment)

### Installation
```bash
conda create -n phylogeny -c bioconda \
    python=3.11 biopython matplotlib numpy iqtree mafft

conda activate phylogeny
```

### Example Data
- Aligned FASTA file (< 100 kb typical)
- Or raw sequences (< 10 kb per sequence)

## Tips for Success

1. **Start small**: Test with 5-10 sequences first
2. **Understand each step**: Don't just run commands, read the output
3. **Check the .iqtree file**: Most important output from tree building
4. **Validate rooting**: Make sure outgroup is actually distant
5. **Interpret distances**: Remember 3% = species boundary
6. **Save your work**: Keep figures and trees for reports

## Next Steps

After learning these scripts:
1. Learn about **gene selection** for DNA barcoding
2. Study **multiple sequence alignment** deeper
3. Explore **other phylogenetic methods** (Bayesian, etc.)
4. Practice on **real datasets** from BOLD or GenBank
5. Write up **complete analyses** with figures and interpretation

## Reference Materials

### Key Papers
- Minh et al. (2020) - IQ-TREE 2 paper
- Altschul et al. (1990) - Bootstrap in phylogenetics
- Hebert et al. (2003) - DNA barcoding (original)

### Online Resources
- IQ-TREE website: www.iqtree.org
- Biopython documentation: biopython.org
- BOLD system: www.boldsystems.org
- DNA Barcoding: www.dna-barcoding.info

## File Statistics

| File | Lines | Purpose |
|------|-------|---------|
| run_iqtree.sh | 345 | IQ-TREE automation |
| visualize_tree.py | 468 | Tree graphics |
| root_tree.py | 525 | Tree manipulation |
| calculate_distances.py | 643 | Distance calculations |
| README.md | 507 | Comprehensive guide |
| QUICK_REFERENCE.md | (varies) | Quick lookup |
| INDEX.md | (this file) | Navigation guide |

**Total**: 2,488+ lines of code and documentation

## Version Information

- **Created**: 2025-11-18
- **Version**: 1.0
- **Compatible with**: Python 3.6+, IQ-TREE 2.0+, BioPython 1.7+
- **Platform**: Linux, macOS, Windows (with WSL)

## Contributing & Feedback

These scripts are educational. Feel free to:
- Modify for your research
- Extend with new features
- Share improvements
- Use in teaching
- Adapt for different organisms/genes

---

**Now go build some trees!** Start with QUICK_REFERENCE.md or README.md, then run your first analysis.
