# DNA Barcoding - Alignment Scripts

Educational scripts for sequence alignment, visualization, and quality assessment. These scripts demonstrate best practices in bioinformatics workflows using industry-standard tools like MAFFT, BioPython, and Matplotlib.

## Scripts Overview

### 1. `run_mafft.sh` - MAFFT Alignment Wrapper

A comprehensive bash wrapper for MAFFT (Multiple Alignment using Fast Fourier Transform), a widely-used multiple sequence alignment program.

**Features:**
- Multiple alignment strategies (fast, balanced, default, accurate)
- Parallel processing with configurable thread count
- Automatic strategy selection based on sequence count
- Progress reporting and detailed logging
- Input validation and error handling

**Usage:**
```bash
# Basic usage - default strategy
./run_mafft.sh -i sequences.fasta -o aligned.fasta

# Fast alignment with 4 threads
./run_mafft.sh -i sequences.fasta -o aligned.fasta -m fast -t 4

# Accurate alignment (slow, high quality)
./run_mafft.sh -i sequences.fasta -o aligned.fasta -m accurate

# Quiet mode (no progress output)
./run_mafft.sh -i sequences.fasta -o aligned.fasta -q
```

**Requirements:**
- MAFFT installed (install via: `brew install mafft` on macOS)
- Bash 4.0+

**MAFFT Strategies Explained:**
- **Fast (FFT-NS-1):** Fastest, least accurate. Use for very large datasets (>5000 sequences)
- **Balanced (FFT-NS-2):** Default balance of speed and accuracy. Good for 100-5000 sequences
- **Default (FFT-NS-2 + refinement):** Good starting point for most alignments
- **Accurate (L-INS-i):** Highest quality, slowest. Best for <200 sequences

**References:**
- MAFFT Manual: https://mafft.cbrc.jp/alignment/software/manual/manual.html
- Katoh et al. (2018) MSB

---

### 2. `alignment_stats.py` - Alignment Statistics Calculator

Calculate comprehensive quality metrics for sequence alignments including conservation scores, gap content, and sequence variation.

**Features:**
- Gap percentage analysis (overall and per-position)
- Shannon entropy for measuring sequence variation
- Conservation scores (0 = variable, 1 = conserved)
- Overall sequence identity
- Per-sequence statistics (gap content, length)
- CSV export of position-by-position metrics
- Support for DNA and protein sequences

**Usage:**
```bash
# Basic statistics summary
./alignment_stats.py -i alignment.fasta

# Print detailed per-sequence statistics
./alignment_stats.py -i alignment.fasta --verbose

# Export position statistics to CSV
./alignment_stats.py -i alignment.fasta -o position_stats.csv

# Process Stockholm format alignment
./alignment_stats.py -i alignment.sto -f stockholm
```

**Requirements:**
- Python 3.6+
- BioPython: `pip install biopython`
- NumPy: `pip install numpy`

**Understanding the Metrics:**

1. **Gap Percentage:** Percentage of positions containing gaps. Lower is better for phylogenetic analysis.

2. **Shannon Entropy:** Measures sequence variation at each position:
   - Formula: H = -Σ(p_i × log₂(p_i))
   - Range: 0 (highly conserved) to ~2 for DNA, ~4.3 for proteins (highly variable)

3. **Conservation Score:** Normalized entropy (1 - entropy/max_entropy):
   - Range: 0 (variable) to 1 (fully conserved)
   - Used to identify functionally important regions

4. **Sequence Identity:** Percentage of positions where all sequences have the same character.

**References:**
- BioPython AlignIO: https://biopython.org/wiki/AlignIO
- Conservation scoring: https://avrilomics.blogspot.com/2016/07/calculating-conservation-score-for.html

---

### 3. `visualize_alignment.py` - Alignment Visualization

Create publication-quality visualizations of sequence alignments with multiple plot types.

**Features:**
- Sequence alignment matrix (color-coded by nucleotide/amino acid)
- Conservation profile plot (line and fill)
- Gap distribution histogram
- Multiple color schemes (DNA and protein)
- Region-specific visualization
- Export to PNG, PDF, or other formats
- Interactive display option

**Usage:**
```bash
# Basic visualization to PNG
./visualize_alignment.py -i alignment.fasta -o alignment.png

# High-resolution PDF output
./visualize_alignment.py -i alignment.fasta -o alignment.pdf --dpi 300

# Visualize specific region (positions 1-500)
./visualize_alignment.py -i alignment.fasta -o region.png --region 1 500

# Interactive display
./visualize_alignment.py -i alignment.fasta --show
```

**Requirements:**
- Python 3.6+
- BioPython: `pip install biopython`
- NumPy: `pip install numpy`
- Matplotlib: `pip install matplotlib`
- Seaborn (optional): `pip install seaborn`

**Color Schemes:**
- **DNA:** A=Green, C=Blue, G=Orange, T=Red, Gap=Gray
- **Protein:** Hydrophobic=Blue, Aromatic=Orange, Polar=Green, etc.

**References:**
- BioPython: https://biopython.org/
- Matplotlib for bioinformatics: https://www.numberanalytics.com/blog/matplotlib-bioinformatics-pipelines

---

### 4. `trim_alignment.py` - Alignment Trimming and Cleaning

Remove poorly aligned and ambiguous regions from alignments to improve phylogenetic analysis.

**Features:**
- Gap-based trimming (remove high-gap positions and sequences)
- Entropy-based trimming (remove low-conservation regions)
- Multiple preset strategies (strict, moderate, lenient)
- Custom parameter combinations
- Contiguous block detection
- Detailed trimming reports
- Support for various alignment formats

**Usage:**
```bash
# Apply strict trimming preset
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --strict

# Apply moderate trimming (default)
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --moderate

# Custom: remove positions with >50% gaps
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --max-gap 50

# Remove sequences with >30% gaps
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --seq-gap 30

# Keep only well-conserved regions
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --min-conservation 0.7

# Combine custom parameters
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --max-gap 30 --seq-gap 50
```

**Trimming Strategies:**

| Strategy | Position Gap Threshold | Sequence Gap Threshold | Use Case |
|----------|----------------------|----------------------|----------|
| **Strict** | >10% | >30% | Small, high-quality alignments; phylogenetically critical |
| **Moderate** | >30% | >50% | Most phylogenetic analyses (default) |
| **Lenient** | >50% | >70% | Large alignments with many indels; exploratory analysis |

**Requirements:**
- Python 3.6+
- BioPython: `pip install biopython`
- NumPy: `pip install numpy`

**Important Considerations:**

1. **Information Loss:** Trimming removes data. Consider downstream analyses when choosing thresholds.

2. **Alignment Quality:** Poorly aligned regions can bias phylogenetic inference more than missing data.

3. **Comparative Studies:** trimAl studies show it outperforms Gblocks in ~90% of scenarios.

**References:**
- trimAl paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC2712344/
- Comparison: https://pmc.ncbi.nlm.nih.gov/articles/PMC4538881/
- pytrimal: https://github.com/althonos/pytrimal

---

## Typical Workflow

```bash
# 1. Align sequences
./run_mafft.sh -i unaligned.fasta -o aligned.fasta -m accurate

# 2. Visualize alignment quality
./visualize_alignment.py -i aligned.fasta -o alignment_viz.png

# 3. Calculate statistics
./alignment_stats.py -i aligned.fasta -o stats.csv --verbose

# 4. Trim poorly aligned regions
./trim_alignment.py -i aligned.fasta -o trimmed.fasta --moderate

# 5. Visualize trimmed alignment
./visualize_alignment.py -i trimmed.fasta -o trimmed_viz.png

# 6. Verify quality of trimmed alignment
./alignment_stats.py -i trimmed.fasta
```

---

## Installation and Setup

### Prerequisites

Install BioPython and required dependencies:
```bash
# Using pip
pip install biopython numpy matplotlib

# Using conda (recommended for scientific computing)
conda install -c bioconda biopython numpy matplotlib

# Install MAFFT
# macOS
brew install mafft

# Linux (Ubuntu/Debian)
sudo apt-get install mafft

# Linux (Fedora)
sudo dnf install mafft
```

### Running Scripts

1. **Make executable:**
```bash
chmod +x *.py *.sh
```

2. **Run with Python interpreter:**
```bash
python3 alignment_stats.py -i alignment.fasta
```

3. **Or directly (if shebang is present):**
```bash
./alignment_stats.py -i alignment.fasta
```

---

## Common Issues and Troubleshooting

### MAFFT Not Found
```bash
# Check if MAFFT is installed
which mafft

# Install MAFFT
brew install mafft  # macOS
sudo apt-get install mafft  # Ubuntu
```

### BioPython Import Error
```bash
# Install BioPython
pip install --upgrade biopython
```

### Permission Denied
```bash
# Make scripts executable
chmod +x *.py *.sh
```

### Matplotlib Display Issues
For remote servers, use non-interactive backend:
```bash
export MPLBACKEND=Agg
./visualize_alignment.py -i alignment.fasta -o output.png
```

---

## Educational Resources

- **Multiple Sequence Alignment:** https://www.ebi.ac.uk/training/online/courses/multiple-sequence-alignment
- **BioPython Tutorial:** https://biopython.org/wiki/Documentation
- **Sequence Analysis:** https://www.ncbi.nlm.nih.gov/pubmed/25402006
- **Phylogenetics:** https://crsl.soe.ucsc.edu/

---

## License and Attribution

These scripts are educational examples for the DNA Barcoding Analysis course. They demonstrate best practices for scientific computing and alignment quality assessment.

**Sources Referenced:**
- MAFFT (Katoh et al.)
- BioPython (Cock et al.)
- trimAl (Capella-Gutiérrez et al.)
- Matplotlib and NumPy developers

---

## Script Maintenance Notes

- Scripts are designed to be self-documenting with extensive comments
- Each script includes usage examples in the header
- Error handling and input validation included
- Compatible with Python 3.6+ and current tool versions (as of 2024)

For questions or improvements, refer to the documentation and research papers linked in each script.
