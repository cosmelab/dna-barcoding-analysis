# Phylogenetic Analysis Scripts

A collection of educational scripts for phylogenetic tree construction, visualization, and analysis. These scripts are designed for DNA barcoding analysis, with extensive comments explaining phylogenetic concepts.

## Overview

This directory contains four complementary scripts that form a complete workflow:

```
Raw Sequences (FASTA)
        ↓
[MAFFT Alignment]
        ↓
Aligned FASTA
        ↓
[run_iqtree.sh] ──────→ Phylogenetic Tree
        ↓                      ↓
    Tree File ────→ [visualize_tree.py] ──→ PNG/PDF Figure
        ↓                      ↓
    [root_tree.py] ──→ Rooted Tree with directions
        ↓
    [calculate_distances.py] ──→ Species identification
```

## Quick Start

### 1. Run IQ-TREE Analysis

```bash
# Basic analysis (recommended for learning)
./run_iqtree.sh aligned_sequences.fasta

# With more bootstrap replicates (more robust)
./run_iqtree.sh aligned_sequences.fasta -bb 2000

# Specify model explicitly
./run_iqtree.sh aligned_sequences.fasta -m GTR+G+I
```

**Outputs:**
- `aligned_sequences_iqtree.treefile` - Your phylogenetic tree (most important!)
- `aligned_sequences_iqtree.iqtree` - Detailed analysis report
- `aligned_sequences_iqtree.log` - Run log

### 2. Visualize the Tree

```bash
# Display tree on screen
python visualize_tree.py aligned_sequences_iqtree.treefile

# Save as PNG
python visualize_tree.py aligned_sequences_iqtree.treefile -o tree.png

# Publication-quality PDF (high resolution)
python visualize_tree.py aligned_sequences_iqtree.treefile -o tree.pdf --dpi 300

# Larger figure (for many species)
python visualize_tree.py aligned_sequences_iqtree.treefile -w 14 -h 10
```

### 3. Root the Tree (Understand Evolutionary History)

```bash
# Midpoint rooting (use when no outgroup known)
python root_tree.py aligned_sequences_iqtree.treefile midpoint -o rooted_tree.treefile

# Outgroup rooting (use when you know distantly related species)
python root_tree.py aligned_sequences_iqtree.treefile outgroup \
    --outgroup "Species_A" --outgroup "Species_B" \
    -o rooted_tree.treefile

# List all species in tree
python root_tree.py aligned_sequences_iqtree.treefile --list-species
```

### 4. Calculate Genetic Distances (Species Identification)

```bash
# From tree branch lengths
python calculate_distances.py -t aligned_sequences_iqtree.treefile

# From sequence alignment (more direct)
python calculate_distances.py -a aligned_sequences.fasta

# With species threshold analysis
python calculate_distances.py -a aligned_sequences.fasta --threshold-analysis

# Save to spreadsheet
python calculate_distances.py -a aligned_sequences.fasta -o distances.csv
```

## Scripts in Detail

### run_iqtree.sh

IQ-TREE wrapper for phylogenetic tree construction.

**What it does:**
- Validates input alignment file
- Automatically selects best substitution model (MFP)
- Builds phylogenetic tree using maximum likelihood
- Estimates bootstrap support (confidence values)
- Produces formatted output with explanations

**Why each parameter matters:**
- `-m MFP`: Automatic model selection (tests multiple models, picks best)
- `-bb 1000`: 1000 bootstrap replicates (tests confidence)
- `-nt AUTO`: Uses all available CPU cores (faster analysis)

**Output files:**
- `.treefile`: Best ML tree (Newick format) - **USE THIS FOR DOWNSTREAM**
- `.iqtree`: Full analysis report with statistics
- `.log`: Run log (useful for debugging)

**Learning objectives:**
- Understand what phylogenetic models are
- Interpret bootstrap support values (>95% = strong, >70% = moderate)
- Know how ML estimation works conceptually
- Recognize when to adjust parameters

**Example:**
```bash
./run_iqtree.sh coi_aligned.fasta
```

Produces a tree file ready for visualization or further analysis.

---

### visualize_tree.py

Create beautiful, informative phylogenetic tree figures.

**What it does:**
- Reads tree file in Newick format
- Draws tree structure with branch lengths
- Labels nodes with bootstrap values
- Colors branches by support level
- Produces publication-quality graphics

**Key features:**
- Interactive display or file saving
- Customizable figure size and resolution
- Bootstrap value visualization
- Informative legends and statistics
- Support for multiple output formats (PNG, PDF, etc.)

**Understanding the visualization:**
- **Tree structure**: Shows evolutionary relationships (topology)
- **Branch lengths**: Horizontal distance = evolutionary divergence
- **Bootstrap numbers**: Confidence in that grouping (0-100%)
- **Colors**: Green = well-supported, Red = less supported

**Example usage:**
```bash
# Quick visualization
python visualize_tree.py tree.treefile

# Save as PNG
python visualize_tree.py tree.treefile -o my_tree.png

# Large figure for 50+ species
python visualize_tree.py tree.treefile -w 16 -h 12 -o large_tree.pdf

# Publication quality
python visualize_tree.py tree.treefile -o figure1.pdf --dpi 300
```

**Learning objectives:**
- Interpret tree topology (who's related to whom)
- Understand bootstrap support
- Read branch length scale
- Create scientific figures

---

### root_tree.py

Reroot phylogenetic trees to show evolutionary direction.

**What it does:**
- Takes unrooted tree from IQ-TREE
- Applies rooting method (midpoint or outgroup)
- Reorganizes tree to show ancestor-descendant relationships
- Ladderizes tree for better visualization
- Saves rooted tree for downstream analysis

**Why rooting matters:**
- Unrooted tree: Shows relationships but not direction
- Rooted tree: Shows which species are ancestors/descendants
- Root position affects evolutionary interpretation

**Rooting methods:**

1. **Midpoint rooting:**
   - Root at the middle of longest path
   - Use when you have no outgroup information
   - Assumes roughly equal evolutionary rates

   ```bash
   python root_tree.py tree.treefile midpoint -o rooted_tree.treefile
   ```

2. **Outgroup rooting:**
   - Specify species known to be distantly related
   - Most accurate if you know your outgroups
   - Root placed between outgroup and everything else

   ```bash
   python root_tree.py tree.treefile outgroup \
       --outgroup "Outgroup_species_1" \
       --outgroup "Outgroup_species_2" \
       -o rooted_tree.treefile
   ```

**Example:**
```bash
# List all species to choose outgroup
python root_tree.py tree.treefile --list-species

# Root with outgroup
python root_tree.py tree.treefile outgroup \
    --outgroup "distant_species" \
    -o rooted_tree.treefile

# Visualize rooted tree
python visualize_tree.py rooted_tree.treefile
```

**Learning objectives:**
- Understand rooted vs unrooted trees
- Know when to use different rooting methods
- Recognize outgroup concept
- See how rooting changes tree interpretation

---

### calculate_distances.py

Calculate evolutionary distances between sequences for species identification.

**What it does:**
- Computes all pairwise distances between sequences
- Multiple distance calculation methods
- Analyzes distribution for species threshold
- Saves results to CSV for further analysis

**Why genetic distances matter:**
- Species differ by ~3% (COI gene, DNA barcoding)
- Individuals of same species ~0.3% different
- Creates "barcoding gap" for identification
- Distance < 3% = probably same species

**Distance calculation methods:**

1. **P-distance** (simplest):
   - Count positions where sequences differ
   - Divide by total length
   - Good for closely related sequences

   ```bash
   python calculate_distances.py -a alignment.fasta -m p-distance
   ```

2. **Jukes-Cantor** (corrects for hidden changes):
   - Accounts for multiple substitutions at same site
   - Assumes all changes equally likely
   - Good for moderately divergent sequences

   ```bash
   python calculate_distances.py -a alignment.fasta -m jukes-cantor
   ```

3. **Kimura 2-parameter** (most realistic for DNA):
   - Distinguishes transitions (A↔G, C↔T) from transversions
   - Transitions happen more often in real DNA
   - **Recommended for DNA barcoding**

   ```bash
   python calculate_distances.py -a alignment.fasta -m kimura-2p
   ```

**Example workflow for species identification:**
```bash
# Calculate distances
python calculate_distances.py -a alignment.fasta --threshold-analysis

# Interpret output:
# If your unknown sequence's closest distance is:
#   - < 0.01 (1%): Very likely same species
#   - < 0.03 (3%): Probably same species
#   - > 0.03 (3%): Different species
```

**Learning objectives:**
- Understand genetic distance metrics
- Know when to apply distance corrections
- Interpret distances for species identification
- Recognize the "barcoding gap" concept

---

## Phylogenetic Concepts Explained

### Bootstrap Support

Numbers on tree branches (0-100%) indicate confidence in that grouping.

- **>95%**: Strong support - species definitely distinct
- **70-95%**: Moderate support - probably distinct
- **<70%**: Weak support - uncertain

Why it matters: Low bootstrap doesn't mean wrong, just uncertain.

### Branch Length

Distance from root of branch = evolutionary time (in substitutions per site).

In DNA barcoding:
- Same species: usually 0.01-0.05 changes per site
- Different species: usually 0.05-0.5 changes per site

### Model Selection

IQ-TREE tests multiple substitution models and picks the best.

- **GTR+G**: General Time Reversible + Gamma distribution (usually best for DNA)
- **TN93+G**: Transition/Transversion model
- **HKY+G**: Simpler model, faster calculation

For DNA barcoding, GTR+G is usually optimal.

### Barcoding Gap

The difference between within-species and between-species distances.

```
Within-species:  ~0.003 (0.3%)
Barcoding gap:   0.003 - 0.030
Between-species: ~0.030 (3%)
```

A clear gap = good for species identification!

---

## Workflow Example: Complete Analysis

```bash
#!/bin/bash
# Complete phylogenetic analysis workflow

# Start with aligned sequences (from MAFFT)
ALIGNMENT="coi_aligned.fasta"

# Step 1: Build phylogenetic tree
./run_iqtree.sh $ALIGNMENT -bb 1000
TREE="${ALIGNMENT%.fasta}_iqtree.treefile"

# Step 2: Visualize the tree
python visualize_tree.py $TREE -o coi_tree.png

# Step 3: Root the tree (assuming we know an outgroup)
python root_tree.py $TREE outgroup \
    --outgroup "Outgroup_species" \
    -o coi_rooted.treefile

# Step 4: Visualize rooted tree
python visualize_tree.py coi_rooted.treefile -o coi_rooted_tree.png

# Step 5: Calculate distances for species identification
python calculate_distances.py -a $ALIGNMENT \
    --threshold-analysis \
    -o coi_distances.csv

# Step 6: Identify your unknown species
# Load distances.csv in spreadsheet
# Find closest match to your unknown sequence
# If distance < 0.03 = same species!
```

---

## Tips for Educational Use

### For Students Learning Phylogenetics

1. **Start with small datasets (3-5 species)**
   - Easy to visualize and understand
   - Bootstrap values clear

2. **Run with default parameters first**
   - Don't change model selection or bootstrap
   - Understand standard workflow

3. **Read the .iqtree file**
   - Contains all the important information
   - Model selection results
   - Tree statistics

4. **Compare different rooting methods**
   - See how root position changes interpretation
   - Understand why outgroup matters

5. **Use the --concepts flag**
   ```bash
   python visualize_tree.py --concepts
   python root_tree.py --concepts
   python calculate_distances.py --concepts
   ```

### For Instructors

- These scripts are designed to be self-explanatory
- Heavy commenting explains concepts
- Help text available with `-h` flag
- Can be run without understanding all parameters
- Graduated complexity: basic → advanced features

---

## Troubleshooting

### IQ-TREE not found
```bash
# Install with:
conda install -c bioconda iqtree

# Verify installation:
iqtree --version
```

### Python dependencies missing
```bash
# Install all requirements:
conda install biopython dendropy pandas numpy matplotlib
```

### Alignment file format issue
- Must be FASTA format (lines starting with `>`)
- Sequences must be aligned (same length)
- No special characters in sequence names
- Use MAFFT to align: `mafft sequences.fasta > aligned.fasta`

### Tree visualization not showing
- Try saving to file instead: `-o tree.png`
- Check figure size with `-w` and `-h`
- May need X11 for interactive display on remote servers

### Bootstrap values missing
- Check .iqtree file was created
- IQ-TREE may not have completed (check .log file)
- Tree may have failed - rerun with more verbose output

---

## Further Reading

### Official Documentation
- **IQ-TREE**: http://www.iqtree.org
- **Biopython**: https://biopython.org
- **DNA Barcoding**: http://www.boldinterpretation.org

### Key Concepts
- Maximum Likelihood: The most probable tree given your data
- Bootstrap: Resampling to estimate confidence
- Substitution Models: Mathematical models of DNA evolution
- Genetic Distance: Measure of sequence divergence

### For Publication
- Always cite IQ-TREE: Minh et al. (2020) Molecular Biology and Evolution
- Report model selection results
- Show bootstrap values on figure
- Report branch length scale
- Include distance table in supplementary

---

## Citation

If you use these scripts in publications:

```
Phylogenetic Analysis Scripts for DNA Barcoding
Educational Materials for DNA Barcoding Analysis
http://www.iqtree.org
https://biopython.org
```

---

## Author Notes

These scripts are intentionally over-commented for educational purposes. Each script explains:

1. **What** the code does
2. **Why** it does that (phylogenetic concepts)
3. **How** to use it (examples)
4. **What** the output means

Remove comments in production code, but keep them for teaching!

---

**Last Updated**: 2025-11-18
**Version**: 1.0
**Target Audience**: Beginners learning phylogenetics and DNA barcoding
