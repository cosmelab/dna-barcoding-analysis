# Alignment Scripts - Usage Examples

Complete examples demonstrating the alignment scripts in real-world scenarios.

## Example 1: Complete Alignment and Analysis Workflow

### Scenario: Analyzing cytochrome oxidase (COX1) barcoding sequences

```bash
cd 06_alignment/scripts

# 1. Align unaligned COX1 sequences using high-accuracy method
./run_mafft.sh -i cox1_unaligned.fasta -o cox1_aligned.fasta -m accurate -t 4
# Output: cox1_aligned.fasta (aligned sequences)

# 2. Assess alignment quality
./alignment_stats.py -i cox1_aligned.fasta --verbose -o cox1_stats.csv
# Shows:
#   - Overall alignment statistics
#   - Per-sequence gap content
#   - Position-by-position conservation
#   - Saves metrics to CSV for further analysis

# 3. Visualize the alignment
./visualize_alignment.py -i cox1_aligned.fasta -o cox1_alignment.png --dpi 300
# Creates: Alignment matrix, conservation profile, gap distribution

# 4. Trim poorly aligned regions
./trim_alignment.py -i cox1_aligned.fasta -o cox1_trimmed.fasta --moderate
# Reports positions and sequences removed

# 5. Verify trimmed alignment quality
./alignment_stats.py -i cox1_trimmed.fasta
# Compare stats before/after trimming
```

**Expected Output:**
```
Loaded alignment: 42 sequences, 658 bp
Applied moderate trimming strategy...

ALIGNMENT TRIMMING REPORT
===============================================================================
Original alignment:
  Sequences: 42
  Positions: 658

Trimmed alignment:
  Sequences: 42 (removed 0)
  Positions: 612 (removed 46)

Retention rates:
  Positions: 92.9%
  Sequences: 100.0%
===============================================================================
```

---

## Example 2: Quick Assessment of Pre-aligned Sequences

### Scenario: Checking quality of sequences from GenBank

```bash
# Single command to get alignment statistics
./alignment_stats.py -i genbank_alignment.fasta -o quality_metrics.csv

# View summary statistics
cat genbank_alignment.fasta | ./alignment_stats.py -i /dev/stdin

# Find highly conserved regions (visualization of region 1-200)
./visualize_alignment.py -i genbank_alignment.fasta -o region_viz.pdf --region 1 200
```

---

## Example 3: Batch Processing Multiple Alignments

### Scenario: Processing multiple gene alignments

```bash
#!/bin/bash
# batch_align.sh - Process multiple FASTA files

for gene in COX1 16S 18S matK rbcL; do
    echo "Processing ${gene}..."

    # Align with automatic strategy selection
    ./run_mafft.sh -i ${gene}_unaligned.fasta -o ${gene}_aligned.fasta

    # Calculate statistics
    ./alignment_stats.py -i ${gene}_aligned.fasta -o ${gene}_stats.csv

    # Trim alignment
    ./trim_alignment.py -i ${gene}_aligned.fasta -o ${gene}_trimmed.fasta --moderate

    echo "${gene} complete!"
done
```

---

## Example 4: Comparing Trimming Strategies

### Scenario: Determining best trimming approach for your data

```bash
# Original alignment
cp original_alignment.fasta test_alignment.fasta

# Test different strategies
for strategy in strict moderate lenient; do
    echo "Testing $strategy trimming..."
    ./trim_alignment.py -i test_alignment.fasta \
        -o trimmed_${strategy}.fasta \
        --${strategy}

    # Calculate statistics for comparison
    ./alignment_stats.py -i trimmed_${strategy}.fasta > stats_${strategy}.txt
done

# Compare retention rates
echo "=== Retention Comparison ==="
grep "Retention rates" stats_*.txt
```

**Example Output:**
```
strict:    Positions: 45.2%, Sequences: 95.2%
moderate:  Positions: 76.3%, Sequences: 100.0%
lenient:   Positions: 88.1%, Sequences: 100.0%
```

---

## Example 5: Fine-tuned Custom Trimming

### Scenario: Removing sequences with many gaps but keeping well-aligned regions

```bash
# Remove positions with >25% gaps AND sequences with >40% gaps
./trim_alignment.py \
    -i alignment.fasta \
    -o trimmed_custom.fasta \
    --max-gap 25 \
    --seq-gap 40 \
    --verbose

# Only keep highly conserved positions
./trim_alignment.py \
    -i alignment.fasta \
    -o trimmed_conserved.fasta \
    --min-conservation 0.6
```

---

## Example 6: Large-Scale Alignment (Many Sequences)

### Scenario: Aligning 5000+ environmental barcode sequences

```bash
# Use fast strategy with maximum parallelization
./run_mafft.sh \
    -i large_barcode_set.fasta \
    -o aligned_barcodes.fasta \
    -m fast \
    -t 16 \
    --quiet

# Check results without verbose output
./alignment_stats.py -i aligned_barcodes.fasta

# Use lenient trimming (fast and reasonable)
./trim_alignment.py \
    -i aligned_barcodes.fasta \
    -o trimmed_barcodes.fasta \
    --lenient
```

---

## Example 7: Creating Publication-Ready Figures

### Scenario: Preparing figures for a research paper

```bash
# High-resolution visualization for main figure
./visualize_alignment.py \
    -i alignment.fasta \
    -o figure1_alignment.pdf \
    --dpi 300

# Detailed view of specific gene region (supplementary figure)
./visualize_alignment.py \
    -i alignment.fasta \
    -o supp_figure_region.pdf \
    --dpi 300 \
    --region 100 300

# Export detailed statistics table
./alignment_stats.py \
    -i alignment.fasta \
    -o supplementary_stats.csv \
    --verbose
```

---

## Example 8: Troubleshooting Alignment Issues

### Scenario: Alignment seems off, investigating quality

```bash
# 1. Get detailed statistics
./alignment_stats.py -i questionable_alignment.fasta --verbose

# Examine specific positions:
# - High gaps: Positions with many insertions/deletions
# - Low conservation: Highly variable regions
# - Identity: Overall sequence similarity

# 2. Visualize to identify problem regions
./visualize_alignment.py -i questionable_alignment.fasta -o diagnostic.png

# 3. Check individual sequence quality
grep -c "^>" questionable_alignment.fasta  # Count sequences
awk '/^>/ {n++} END {print n}' questionable_alignment.fasta  # Alternative

# 4. If needed, re-align with different strategy
./run_mafft.sh \
    -i original_unaligned.fasta \
    -o realigned.fasta \
    -m accurate  # Try more accurate strategy
```

---

## Example 9: Integration with Phylogenetic Analysis

### Scenario: Preparing alignment for RAxML/IQ-TREE analysis

```bash
# 1. Align sequences
./run_mafft.sh -i sequences.fasta -o aligned.fasta -m accurate

# 2. Export statistics for methods section
./alignment_stats.py -i aligned.fasta -o methods_stats.csv

# 3. Trim for phylogenetic signal
./trim_alignment.py -i aligned.fasta -o trimmed.fasta --moderate

# 4. Convert to phylip format for RAxML
# (Note: Convert using your preferred tool or BioPython script)
python3 << 'EOF'
from Bio import AlignIO
alignment = AlignIO.read("trimmed.fasta", "fasta")
AlignIO.write(alignment, "trimmed.phy", "phylip-relaxed")
EOF

# 5. Document alignment statistics for paper
echo "Alignment Statistics:" > alignment_report.txt
./alignment_stats.py -i trimmed.fasta >> alignment_report.txt
```

---

## Example 10: Working with Different Sequence Types

### Scenario: Analyzing protein sequences instead of DNA

```bash
# For protein sequences, scripts auto-detect type
./alignment_stats.py -i protein_alignment.fasta --verbose

# Alignment statistics for proteins show:
# - Gap percentages (same as DNA)
# - Conservation scores (calculated from 20 amino acids)
# - Shannon entropy (log₂(20) = 4.32 bits maximum)

# Visualization uses Clustal color scheme for proteins
./visualize_alignment.py -i protein_alignment.fasta -o protein_viz.png

# Trimming works identically
./trim_alignment.py -i protein_alignment.fasta -o protein_trimmed.fasta
```

---

## Tips and Best Practices

### Alignment Strategy Selection

```
Number of Sequences | Recommended Strategy | Typical Runtime
< 50               | accurate             | seconds to minutes
50-500             | balanced/default     | minutes
500-5000           | default              | minutes to hours
> 5000             | fast                 | hours or less
```

### Memory Considerations

- MAFFT uses significant memory for large alignments
- Use multiple threads for better resource utilization
- For very large datasets (>10,000 seqs), consider clustering first

### Quality Assessment Workflow

```
1. Calculate basic stats (gap %, identity)
           ↓
2. Visualize alignment (identify problem regions)
           ↓
3. Consider trimming (remove outliers)
           ↓
4. Recalculate stats (compare before/after)
           ↓
5. Decide on downstream analyses
```

### Common Parameters for Different Use Cases

**Strict phylogenetics (best tree inference):**
```bash
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --strict
```

**Exploratory analysis (keep more data):**
```bash
./trim_alignment.py -i alignment.fasta -o trimmed.fasta --lenient
```

**Automated pipelines (balance speed/quality):**
```bash
./run_mafft.sh -i input.fasta -o output.fasta -m default -t 4 -q
./trim_alignment.py -i output.fasta -o trimmed.fasta --moderate -q
```

---

## Example Data

To test these scripts, you'll need:
- FASTA files with unaligned or aligned sequences
- Minimum: 2 sequences
- Recommended: 10-100 sequences for testing

Create a test file:
```bash
# test_sequences.fasta
cat > test_sequences.fasta << 'EOF'
>seq1
ATCGATCGATCGATCG
>seq2
ATCGAT-GATCGATCG
>seq3
ATCGATCGATC-ATCG
EOF
```

---

## Troubleshooting Common Errors

### "MAFFT not found"
```bash
# Install MAFFT
brew install mafft  # macOS
sudo apt-get install mafft  # Ubuntu
```

### "No module named Bio"
```bash
pip install biopython numpy matplotlib
```

### "Permission denied"
```bash
chmod +x *.py *.sh
```

### "Matplotlib: Could not connect to display"
```bash
export MPLBACKEND=Agg
./visualize_alignment.py -i alignment.fasta -o output.png
```

---

## Further Reading

- **MAFFT Paper:** https://doi.org/10.1093/molbev/msr121
- **Multiple Sequence Alignment Review:** https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3848042/
- **Phylogenetic Best Practices:** https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2872364/
- **BioPython Documentation:** https://biopython.org/wiki/Documentation

---

Last Updated: 2024
These examples are for educational purposes in the DNA Barcoding Analysis course.
