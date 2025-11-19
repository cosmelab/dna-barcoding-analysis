# Module 06: Multiple Sequence Alignment

**Duration**: 2 hours
**Prerequisites**: Module 05

---

## Overview

Align DNA sequences to identify homologous positions. Essential step before phylogenetic analysis.

---

## Why Alignment?

Sequences evolve through:
- **Substitutions** (mutations)
- **Insertions** (additions)
- **Deletions** (losses)

Alignment identifies which bases correspond across sequences despite evolutionary changes.

---

## Tools

### MAFFT (Recommended)
```bash
# Auto mode (best for most cases)
mafft --auto sequences.fasta > aligned.fasta

# More accurate (slower)
mafft --maxiterate 1000 --localpair sequences.fasta > aligned.fasta

# Very fast (less accurate)
mafft --retree 1 sequences.fasta > aligned.fasta
```

### MUSCLE
```bash
muscle -in sequences.fasta -out aligned.fasta
```

### ClustalW
```bash
clustalw -INFILE=sequences.fasta -OUTFILE=aligned.fasta
```

---

## Workflow

```
Input: Unaligned FASTA
        ↓
   Run MAFFT
        ↓
Output: Aligned FASTA
        ↓
Manual Inspection (AliView)
        ↓
Clean Alignment
        ↓
Ready for Phylogeny
```

---

## Example: Align COI Sequences

```bash
#!/bin/bash
# Align COI sequences using MAFFT

INPUT="sequences.fasta"
OUTPUT="aligned_sequences.fasta"

echo "Aligning sequences with MAFFT..."

mafft --auto \
      --thread 4 \
      --adjustdirection \
      $INPUT > $OUTPUT

echo "Alignment complete!"
echo "Output: $OUTPUT"

# Quick stats
echo "Number of sequences: $(grep -c '^>' $OUTPUT)"
echo "Alignment length: $(grep -v '^>' $OUTPUT | head -1 | tr -d '\n' | wc -c)"
```

---

## Checking Alignment Quality

### Visual Inspection with AliView

1. Open aligned FASTA in AliView
2. Check for:
   - Conserved regions (should be aligned)
   - Variable regions (expected)
   - Large gaps (might need adjustment)
   - Incorrect reverse complements

### Automated Checks

```python
from Bio import AlignIO

# Read alignment
alignment = AlignIO.read("aligned.fasta", "fasta")

# Check alignment length
print(f"Alignment length: {alignment.get_alignment_length()}")

# Check for gaps
gap_counts = []
for record in alignment:
    gaps = record.seq.count("-")
    gap_counts.append(gaps)

avg_gaps = sum(gap_counts) / len(gap_counts)
print(f"Average gaps per sequence: {avg_gaps:.1f}")
```

---

## Common Issues

### Problem: Sequences in Different Orientations
```bash
# Solution: Use --adjustdirection flag
mafft --auto --adjustdirection sequences.fasta > aligned.fasta
# Creates _R_ files for reversed sequences
```

### Problem: Large Gaps
```bash
# Solution: Trim alignment
# Use trimAl or manual curation
trimal -in aligned.fasta -out trimmed.fasta -automated1
```

### Problem: Very Divergent Sequences
```bash
# Solution: Use more accurate algorithm
mafft --maxiterate 1000 --localpair sequences.fasta > aligned.fasta
```

---

## Next Steps

Aligned sequences proceed to **Module 07: Phylogeny** →
