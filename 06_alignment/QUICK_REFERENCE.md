# Alignment Scripts - Quick Reference Guide

Quick lookup for common commands and parameters.

## Scripts Overview

| Script | Purpose | Language | Main Use |
|--------|---------|----------|----------|
| `run_mafft.sh` | Multiple sequence alignment | Bash | Align unaligned sequences |
| `alignment_stats.py` | Quality metrics calculation | Python | Assess alignment quality |
| `visualize_alignment.py` | Create visualizations | Python | Publication figures |
| `trim_alignment.py` | Remove poorly aligned regions | Python | Clean alignments |

---

## Command Quick Reference

### 1. MAFFT Alignment - `run_mafft.sh`

```bash
# Basic usage
./run_mafft.sh -i input.fasta -o output.fasta

# With custom strategy
./run_mafft.sh -i input.fasta -o output.fasta -m [fast|balanced|default|accurate]

# With parallel processing
./run_mafft.sh -i input.fasta -o output.fasta -t 4

# Quiet mode
./run_mafft.sh -i input.fasta -o output.fasta -q

# Complete example
./run_mafft.sh -i sequences.fasta -o aligned.fasta -m accurate -t 8
```

**Strategy Selection:**
- `fast`: >5000 sequences, speed priority
- `balanced`: 500-5000 sequences
- `default`: Best general purpose
- `accurate`: <200 sequences, quality priority

---

### 2. Alignment Statistics - `alignment_stats.py`

```bash
# Basic statistics
./alignment_stats.py -i alignment.fasta

# Export to CSV
./alignment_stats.py -i alignment.fasta -o stats.csv

# Verbose (per-sequence details)
./alignment_stats.py -i alignment.fasta --verbose

# Different formats
./alignment_stats.py -i alignment.sto -f stockholm
./alignment_stats.py -i alignment.aln -f clustal

# Complete example
./alignment_stats.py -i alignment.fasta -o stats.csv -f fasta --verbose
```

**Output Metrics:**
- `Gap percentage`: Positions with gaps
- `Sequence identity`: Same character at position
- `Shannon entropy`: Sequence variation measure
- `Conservation score`: 0 (variable) to 1 (conserved)

---

### 3. Alignment Visualization - `visualize_alignment.py`

```bash
# Basic PNG output
./visualize_alignment.py -i alignment.fasta -o output.png

# High-resolution PDF
./visualize_alignment.py -i alignment.fasta -o output.pdf --dpi 300

# Specific region (1-based)
./visualize_alignment.py -i alignment.fasta -o output.png --region 1 500

# Interactive display
./visualize_alignment.py -i alignment.fasta --show

# Complete example
./visualize_alignment.py -i alignment.fasta -o figure.pdf --dpi 300 --region 100 600
```

**Output Components:**
- Alignment matrix (color-coded nucleotides/amino acids)
- Conservation profile (line plot)
- Gap distribution (bar chart)

---

### 4. Alignment Trimming - `trim_alignment.py`

```bash
# Preset strategies
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --strict    # >10% gaps
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --moderate  # >30% gaps
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --lenient   # >50% gaps

# Custom parameters
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --max-gap 25
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --seq-gap 40
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --min-conservation 0.6

# Combined
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --max-gap 30 --seq-gap 50

# Verbose report
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --moderate --verbose
```

**Strategy Comparison:**
```
Strategy  | Position Gap | Sequence Gap | Typical Use
strict    | >10%        | >30%        | Small, critical alignments
moderate  | >30%        | >50%        | Most phylogenetic analyses
lenient   | >50%        | >70%        | Large alignments, exploratory
```

---

## Typical Workflows

### Workflow 1: Full Analysis Pipeline
```bash
# Align → Assess → Visualize → Trim → Verify
./run_mafft.sh -i raw.fasta -o aligned.fasta
./alignment_stats.py -i aligned.fasta
./visualize_alignment.py -i aligned.fasta -o before_trim.png
./trim_alignment.py -i aligned.fasta -o trimmed.fasta --moderate
./alignment_stats.py -i trimmed.fasta
./visualize_alignment.py -i trimmed.fasta -o after_trim.png
```

### Workflow 2: Quick Quality Check
```bash
# Align → Stats → Visualize
./run_mafft.sh -i sequences.fasta -o aligned.fasta -q
./alignment_stats.py -i aligned.fasta
./visualize_alignment.py -i aligned.fasta -o viz.png
```

### Workflow 3: Production Trimming
```bash
# Align with accuracy → Trim conservatively
./run_mafft.sh -i sequences.fasta -o aligned.fasta -m accurate
./trim_alignment.py -i aligned.fasta -o trimmed.fasta --strict
```

---

## Parameter Guide

### MAFFT Strategy Selection Table

| Parameter | Sequences | Speed | Accuracy | Memory |
|-----------|-----------|-------|----------|--------|
| fast | >5000 | Fast | Low | Low |
| balanced | 500-5000 | Medium | Medium | Medium |
| default | 100-1000 | Medium | Good | Medium |
| accurate | <200 | Slow | High | High |

### Trimming Thresholds

**When to use strict:**
- Phylogenetic analysis where tree accuracy is critical
- Small, well-curated datasets
- Publication-quality analysis

**When to use moderate:**
- Standard phylogenetic analyses
- Balanced information retention
- Most general use cases

**When to use lenient:**
- Large datasets (>1000 sequences)
- Exploratory analyses
- When preserving data diversity is important

---

## File Format Support

### Input Formats
```bash
-f fasta      # FASTA format (default)
-f stockholm  # Stockholm format
-f clustal    # Clustal format
-f phylip     # PHYLIP format
```

### Output Formats (visualization)
```bash
-o file.png   # PNG (256 DPI raster)
-o file.pdf   # PDF (vector, high quality)
-o file.jpg   # JPEG (raster)
-o file.svg   # SVG (vector)
```

---

## Common Parameter Combinations

### DNA Barcoding Analysis
```bash
./run_mafft.sh -i barcodes.fasta -o aligned.fasta -m accurate -t 4
./trim_alignment.py -i aligned.fasta -o trimmed.fasta --moderate
./visualization_alignment.py -i trimmed.fasta -o barcode_viz.pdf
```

### Large-Scale Metabarcoding
```bash
./run_mafft.sh -i reads.fasta -o aligned.fasta -m fast -t 16 -q
./trim_alignment.py -i aligned.fasta -o trimmed.fasta --lenient
./alignment_stats.py -i trimmed.fasta -o quality.csv
```

### Phylogenetic Analysis
```bash
./run_mafft.sh -i genes.fasta -o aligned.fasta -m accurate
./alignment_stats.py -i aligned.fasta -o metrics.csv --verbose
./trim_alignment.py -i aligned.fasta -o trimmed.fasta --strict
./visualize_alignment.py -i trimmed.fasta -o phylo.pdf --dpi 300
```

### Protein Analysis
```bash
./run_mafft.sh -i proteins.fasta -o aligned.fasta
./alignment_stats.py -i aligned.fasta --verbose
./visualize_alignment.py -i aligned.fasta -o protein_viz.png
./trim_alignment.py -i aligned.fasta -o trimmed.fasta --moderate
```

---

## Troubleshooting Quick Fixes

| Error | Solution |
|-------|----------|
| MAFFT not found | `brew install mafft` or `apt-get install mafft` |
| Bio module not found | `pip install biopython` |
| Permission denied | `chmod +x *.sh *.py` |
| Cannot display plot | `export MPLBACKEND=Agg` |
| Memory error | Use `-m fast` strategy or reduce dataset size |
| Slow alignment | Use `-t 4` for parallel processing |

---

## Performance Tips

### Speed Up Alignment
```bash
# Use more threads
./run_mafft.sh -i input.fasta -o output.fasta -t 16

# Use faster strategy
./run_mafft.sh -i input.fasta -o output.fasta -m fast
```

### Reduce Memory Usage
```bash
# Process in chunks and realign
# Cluster sequences first, then align clusters

# Use faster strategies (less iterative refinement)
./run_mafft.sh -i input.fasta -o output.fasta -m balanced
```

### Improve Visualization
```bash
# Higher DPI for publications
./visualize_alignment.py -i aligned.fasta -o output.pdf --dpi 300

# Specific regions for clarity
./visualize_alignment.py -i aligned.fasta -o output.png --region 1 200
```

---

## Dependencies Summary

### Required
```bash
pip install biopython numpy matplotlib
```

### Optional
```bash
pip install seaborn  # Better visualizations
```

### System Requirements
```bash
# Install MAFFT
brew install mafft     # macOS
apt-get install mafft  # Ubuntu/Debian
dnf install mafft      # Fedora
```

---

## File Size Guidelines

| Parameter | Typical Size | Max Size |
|-----------|-------------|----------|
| Small alignment | <100 seqs | 1 MB |
| Medium | 100-1000 seqs | 10 MB |
| Large | 1000-10000 seqs | 100 MB |
| Very large | >10000 seqs | >1 GB |

---

## Quick Help

```bash
# Get help for any script
./run_mafft.sh -h
./alignment_stats.py -h
./visualize_alignment.py -h
./trim_alignment.py -h

# Check installed versions
python3 --version
mafft --version
```

---

## Example Data Locations

Test scripts with sequences from:
- **Local examples:** `06_alignment/examples/` (if present)
- **Create test data:**
  ```bash
  cat > test.fasta << 'EOF'
  >seq1
  ATCGATCGATCGATCGATCGATCG
  >seq2
  ATCGATCGATCGATCGATCGATCG
  EOF
  ```

---

## Reference Materials

- Full documentation: See `README.md`
- Detailed examples: See `EXAMPLES.md`
- Script headers contain research citations

Last Updated: 2024
Part of DNA Barcoding Analysis Course
