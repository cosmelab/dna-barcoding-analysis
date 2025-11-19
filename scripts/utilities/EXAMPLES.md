# Utility Scripts - Practical Examples

This document provides real-world examples of how to use the bioinformatics utility scripts.

## Table of Contents
1. [Basic File Operations](#basic-file-operations)
2. [Sequence Analysis](#sequence-analysis)
3. [Data Preparation](#data-preparation)
4. [Advanced Workflows](#advanced-workflows)
5. [Command-Line Usage](#command-line-usage)

---

## Basic File Operations

### Example 1.1: Count Sequences in a File

**Scenario**: You have a FASTA file and want to know how many sequences it contains.

```python
#!/usr/bin/env python3
"""Count sequences in FASTA files."""

from utilities.fasta_tools import count_sequences
import os

# Single file
fasta_file = "sequences.fasta"
num_seqs = count_sequences(fasta_file)
print(f"{fasta_file}: {num_seqs} sequences")

# Multiple files in a directory
data_dir = "data"
total_seqs = 0
for filename in os.listdir(data_dir):
    if filename.endswith('.fasta'):
        filepath = os.path.join(data_dir, filename)
        count = count_sequences(filepath)
        total_seqs += count
        print(f"  {filename}: {count}")

print(f"Total: {total_seqs} sequences")
```

### Example 1.2: Read and Display Sequences

**Scenario**: You want to preview sequences from a file.

```python
#!/usr/bin/env python3
"""Display sequences from a FASTA file."""

from utilities.fasta_tools import read_fasta

# Display first 5 sequences
fasta_file = "sequences.fasta"
count = 0

for header, sequence in read_fasta(fasta_file):
    print(f">Header: {header}")
    print(f"Length: {len(sequence)} bp")
    print(f"Sequence: {sequence[:80]}...")  # First 80 characters
    print()

    count += 1
    if count >= 5:
        break
```

### Example 1.3: Extract Specific Sequence by ID

**Scenario**: You need to get a single sequence by its ID.

```python
#!/usr/bin/env python3
"""Extract specific sequence."""

from utilities.fasta_tools import get_sequence_by_id

# Find and display a specific sequence
fasta_file = "sequences.fasta"
target_id = "seq_001"

sequence = get_sequence_by_id(fasta_file, target_id)

if sequence:
    print(f"Found {target_id}:")
    print(f"Length: {len(sequence)} bp")
    print(f"Sequence: {sequence}")
else:
    print(f"Sequence '{target_id}' not found")
```

### Example 1.4: Split Large FASTA File

**Scenario**: Your FASTA file is too large and you want to split it for parallel processing.

```python
#!/usr/bin/env python3
"""Split large FASTA file into manageable chunks."""

from utilities.fasta_tools import split_fasta

# Split into files with 500 sequences each
input_file = "large_database.fasta"
output_prefix = "chunk"

results = split_fasta(input_file, output_prefix, sequences_per_file=500)

print(f"Split {input_file} into {len(results)} files:")
for filename, count in sorted(results.items()):
    print(f"  {filename}: {count} sequences")
```

### Example 1.5: Merge Multiple FASTA Files

**Scenario**: You have multiple FASTA files from different sources and want to combine them.

```python
#!/usr/bin/env python3
"""Merge multiple FASTA files."""

from utilities.fasta_tools import merge_fasta
import glob

# Find all FASTA files in directory
fasta_files = glob.glob("datasets/*.fasta")
output_file = "merged_all.fasta"

# Merge them
total = merge_fasta(fasta_files, output_file)
print(f"Merged {len(fasta_files)} files → {total} total sequences")
```

---

## Sequence Analysis

### Example 2.1: Calculate GC Content

**Scenario**: You want to analyze the GC content of all sequences in a file.

```python
#!/usr/bin/env python3
"""Calculate GC content for all sequences."""

from utilities.fasta_tools import read_fasta
from utilities.sequence_stats import SequenceStats
import statistics

fasta_file = "sequences.fasta"
gc_contents = []
gc_outliers = []

for header, sequence in read_fasta(fasta_file):
    stats = SequenceStats(sequence)
    gc = stats.gc_content()
    gc_contents.append(gc)

    # Flag sequences with unusual GC content
    if gc < 30 or gc > 60:
        gc_outliers.append((header, gc))

# Statistics
mean_gc = statistics.mean(gc_contents)
stdev_gc = statistics.stdev(gc_contents)

print(f"GC Content Analysis ({len(gc_contents)} sequences)")
print(f"  Mean: {mean_gc:.1f}%")
print(f"  Std Dev: {stdev_gc:.1f}%")
print(f"  Range: {min(gc_contents):.1f}% - {max(gc_contents):.1f}%")

if gc_outliers:
    print(f"\nSequences with unusual GC content:")
    for header, gc in gc_outliers[:10]:
        print(f"  {header}: {gc:.1f}%")
```

### Example 2.2: Analyze Sequence Composition

**Scenario**: You want to understand nucleotide composition of sequences.

```python
#!/usr/bin/env python3
"""Analyze nucleotide composition."""

from utilities.fasta_tools import read_fasta
from utilities.sequence_stats import SequenceStats

fasta_file = "sequences.fasta"

# Aggregate composition across all sequences
total_composition = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
total_length = 0

for header, sequence in read_fasta(fasta_file):
    stats = SequenceStats(sequence)
    comp = stats.nucleotide_composition()

    for base, count in comp.items():
        if base in total_composition:
            total_composition[base] += count

    total_length += len(sequence)

# Calculate percentages
print(f"Overall Nucleotide Composition ({total_length:,} bp)")
print("-" * 40)

for base in ['A', 'T', 'G', 'C', 'N']:
    count = total_composition[base]
    percentage = (count / total_length) * 100
    print(f"{base}: {count:>10,} ({percentage:>5.1f}%)")

# Verify complementary bases
print(f"\nVerification:")
at_total = total_composition['A'] + total_composition['T']
gc_total = total_composition['G'] + total_composition['C']
print(f"A+T: {at_total:,} ({100*at_total/total_length:.1f}%)")
print(f"G+C: {gc_total:,} ({100*gc_total/total_length:.1f}%)")
```

### Example 2.3: Find Sequences with Specific Properties

**Scenario**: You want to find sequences with specific characteristics.

```python
#!/usr/bin/env python3
"""Find sequences with specific properties."""

from utilities.fasta_tools import read_fasta, write_fasta
from utilities.sequence_stats import SequenceStats

fasta_file = "sequences.fasta"
output_file = "filtered_sequences.fasta"

# Criteria
min_length = 400
max_length = 700
max_gc_content = 55
min_at_content = 35

filtered_sequences = []

for header, sequence in read_fasta(fasta_file):
    stats = SequenceStats(sequence)

    # Check all criteria
    if not (min_length <= len(sequence) <= max_length):
        continue

    gc = stats.gc_content()
    if gc > max_gc_content:
        continue

    at = stats.at_content()
    if at < min_at_content:
        continue

    filtered_sequences.append((header, sequence))

# Write results
write_fasta(output_file, filtered_sequences)

print(f"Filtered Results:")
print(f"  Original: {len(list(read_fasta(fasta_file)))} sequences")
print(f"  Filtered: {len(filtered_sequences)} sequences")
print(f"  Retention: {100*len(filtered_sequences)/len(list(read_fasta(fasta_file))):.1f}%")
```

### Example 2.4: Check Sequence Validity

**Scenario**: You want to validate sequence data before analysis.

```python
#!/usr/bin/env python3
"""Validate sequences and report issues."""

from utilities.fasta_tools import read_fasta
from utilities.sequence_stats import SequenceStats

fasta_file = "sequences.fasta"

valid_count = 0
invalid_sequences = []

for header, sequence in read_fasta(fasta_file):
    stats = SequenceStats(sequence)

    # Check validity
    is_valid, message = stats.is_valid_dna(), "Valid"

    if not is_valid:
        ambiguous = stats.ambiguous_nucleotides()
        invalid_sequences.append((header, ambiguous))
        continue

    # Check for excessive ambiguous bases
    ambiguous = stats.ambiguous_nucleotides()
    if len(ambiguous) > len(sequence) * 0.05:  # >5% ambiguous
        invalid_sequences.append((header, ambiguous))
        continue

    valid_count += 1

# Report
print(f"Validation Report")
print(f"  Valid sequences: {valid_count}")
print(f"  Invalid sequences: {len(invalid_sequences)}")

if invalid_sequences:
    print(f"\nInvalid sequences (first 10):")
    for header, ambiguous in invalid_sequences[:10]:
        print(f"  {header}")
        print(f"    Issues: {len(ambiguous)} ambiguous bases")
```

---

## Data Preparation

### Example 3.1: Convert FASTQ to FASTA

**Scenario**: You have quality-controlled reads from a sequencer (FASTQ) and want to convert to FASTA for database building.

```python
#!/usr/bin/env python3
"""Convert FASTQ sequences to FASTA format."""

from utilities.format_converter import FormatConverter

# Convert FASTQ to FASTA
converter = FormatConverter()
converter.fastq_to_fasta('reads.fastq', 'reads.fasta')

print("Conversion complete: reads.fastq → reads.fasta")
```

### Example 3.2: Prepare Sequences for Alignment

**Scenario**: You have sequences that need to be aligned for phylogenetic analysis.

```python
#!/usr/bin/env python3
"""Prepare sequences for alignment."""

from utilities.fasta_tools import filter_fasta
from utilities.sequence_stats import SequenceStats, analyze_sequence_file

input_file = "sequences.fasta"
filtered_file = "sequences_for_alignment.fasta"

# Step 1: Filter by length (must be similar for good alignment)
# Most tools require >80% sequence length similarity
filter_fasta(input_file, filtered_file, min_length=500)

# Step 2: Check what we have
analyze_sequence_file(filtered_file)

# Step 3: Convert to PHYLIP for phylogenetics software
from utilities.format_converter import auto_convert
auto_convert(filtered_file, 'aligned.phy', 'PHYLIP')

print("\nNext steps:")
print("1. Align sequences with MAFFT or Clustal")
print("2. Convert aligned FASTA to PHYLIP")
print("3. Run phylogenetic analysis (RAxML, PhyML, etc.)")
```

### Example 3.3: Create Training and Test Sets

**Scenario**: You want to create training/test data for machine learning.

```python
#!/usr/bin/env python3
"""Create training and test datasets."""

from utilities.subsample_sequences import SequenceSubsampler

# Load sequences
subsampler = SequenceSubsampler(seed=42)  # Fixed seed for reproducibility
all_sequences = subsampler.read_fasta_with_metadata('sequences.fasta')

print(f"Total sequences: {len(all_sequences)}")

# Step 1: Remove duplicates
unique_sequences, n_dupes = subsampler.remove_duplicates(all_sequences)
print(f"After deduplication: {len(unique_sequences)} "
      f"(removed {n_dupes} duplicates)")

# Step 2: Split into training (80%) and test (20%) sets
train_seqs = subsampler.random_sample(unique_sequences, fraction=0.8)
test_indices = set(id(s) for s in train_seqs)
test_seqs = [s for s in unique_sequences if id(s) not in test_indices]

print(f"Training set: {len(train_seqs)}")
print(f"Test set: {len(test_seqs)}")

# Step 3: Write datasets
subsampler.write_fasta('train_data.fasta', train_seqs)
subsampler.write_fasta('test_data.fasta', test_seqs)

# Step 4: Verify
train_stats = subsampler.get_sequence_statistics(train_seqs)
test_stats = subsampler.get_sequence_statistics(test_seqs)

print(f"\nTraining set statistics:")
print(f"  Length: {train_stats['min_length']}-{train_stats['max_length']} bp")
print(f"  Mean: {train_stats['mean_length']:.0f} bp")

print(f"\nTest set statistics:")
print(f"  Length: {test_stats['min_length']}-{test_stats['max_length']} bp")
print(f"  Mean: {test_stats['mean_length']:.0f} bp")
```

### Example 3.4: Create Non-Redundant Database

**Scenario**: You want to reduce sequence redundancy before building a reference database.

```python
#!/usr/bin/env python3
"""Create a non-redundant sequence database."""

from utilities.fasta_tools import merge_fasta
from utilities.subsample_sequences import SequenceSubsampler

# Step 1: Merge all source files
source_files = [
    'ncbi_sequences.fasta',
    'custom_sequences.fasta',
    'additional_data.fasta'
]
merged_file = 'merged_all.fasta'
merge_fasta(source_files, merged_file)

# Step 2: Deduplicate
subsampler = SequenceSubsampler()
sequences = subsampler.read_fasta_with_metadata(merged_file)
unique_seqs, n_dupes = subsampler.remove_duplicates(sequences)

print(f"Merged {len(sequences)} sequences → {len(unique_seqs)} unique")

# Step 3: Reduce redundancy (keep sequences <95% similar)
nr_seqs = subsampler.subsample_by_identity(unique_seqs, max_identity=0.95)

print(f"After clustering: {len(nr_seqs)} non-redundant sequences")

# Step 4: Write final database
subsampler.write_fasta('reference_database_nr.fasta', nr_seqs)

print(f"\nFinal database: reference_database_nr.fasta")
print(f"Compression: {100*(1 - len(nr_seqs)/len(sequences)):.1f}%")
```

---

## Advanced Workflows

### Example 4.1: Complete DNA Barcoding Pipeline

**Scenario**: Analyze barcode sequences from raw reads to final alignment.

```python
#!/usr/bin/env python3
"""Complete DNA barcoding analysis pipeline."""

import sys
from pathlib import Path
from utilities.fasta_tools import read_fasta, write_fasta, filter_fasta
from utilities.sequence_stats import SequenceStats
from utilities.format_converter import auto_convert
from utilities.subsample_sequences import SequenceSubsampler

# Configuration
INPUT_FASTQ = "raw_reads.fastq"
MIN_SEQ_LENGTH = 600
MAX_SEQ_LENGTH = 700
MAX_AMBIGUOUS = 0.05  # Max 5% N's

def run_pipeline():
    """Execute complete pipeline."""

    print("DNA Barcoding Analysis Pipeline")
    print("=" * 50)

    # Step 1: Convert FASTQ to FASTA
    print("\n[1/5] Converting FASTQ to FASTA...")
    auto_convert(INPUT_FASTQ, 'sequences.fasta')

    # Step 2: Quality control filtering
    print("\n[2/5] Filtering sequences...")
    filtered_seqs = []

    for header, sequence in read_fasta('sequences.fasta'):
        stats = SequenceStats(sequence)

        # Check length
        if not (MIN_SEQ_LENGTH <= len(sequence) <= MAX_SEQ_LENGTH):
            continue

        # Check ambiguous bases
        ambiguous = stats.ambiguous_nucleotides()
        if len(ambiguous) > len(sequence) * MAX_AMBIGUOUS:
            continue

        # Check GC content (flag unusual)
        gc = stats.gc_content()
        if gc < 30 or gc > 60:
            print(f"  Warning: {header} has unusual GC content ({gc:.1f}%)")

        filtered_seqs.append((header, sequence))

    write_fasta('filtered.fasta', filtered_seqs)
    print(f"  Passed QC: {len(filtered_seqs)} sequences")

    # Step 3: Deduplicate
    print("\n[3/5] Removing duplicates...")
    subsampler = SequenceSubsampler()
    sequences = subsampler.read_fasta_with_metadata('filtered.fasta')
    unique_seqs, n_dupes = subsampler.remove_duplicates(sequences)
    subsampler.write_fasta('unique.fasta', unique_seqs)
    print(f"  Removed {n_dupes} duplicates")

    # Step 4: Create representative sample
    print("\n[4/5] Creating stratified sample for alignment...")
    sample = subsampler.stratified_sample(unique_seqs, fraction=0.1)
    subsampler.write_fasta('sample_for_alignment.fasta', sample)
    print(f"  Selected {len(sample)} representatives")

    # Step 5: Prepare for phylogenetics
    print("\n[5/5] Preparing for phylogenetic analysis...")
    auto_convert('aligned_sample.fasta', 'aligned.phy', 'PHYLIP')

    # Summary
    print("\n" + "=" * 50)
    print("Pipeline Complete!")
    print("\nOutput files:")
    print("  - unique.fasta: Non-redundant sequences")
    print("  - sample_for_alignment.fasta: Representatives for alignment")
    print("  - aligned.phy: Ready for phylogenetics software")
    print("\nNext steps:")
    print("  1. Align sample_for_alignment.fasta with MAFFT")
    print("  2. Convert aligned FASTA to PHYLIP")
    print("  3. Run phylogenetic analysis (RAxML, PhyML, IQTree)")

if __name__ == "__main__":
    run_pipeline()
```

### Example 4.2: Comparative Analysis of Multiple Datasets

**Scenario**: Compare GC content and composition across different datasets.

```python
#!/usr/bin/env python3
"""Compare multiple sequence datasets."""

import os
import statistics
from utilities.fasta_tools import read_fasta
from utilities.sequence_stats import SequenceStats

# Datasets to compare
datasets = {
    'bacteria': 'bacteria_sequences.fasta',
    'archaea': 'archaea_sequences.fasta',
    'mitochondria': 'mitochondrial_sequences.fasta'
}

print("Comparative Dataset Analysis")
print("=" * 60)

results = {}

for dataset_name, fasta_file in datasets.items():
    if not os.path.exists(fasta_file):
        print(f"Warning: {fasta_file} not found")
        continue

    gc_contents = []
    lengths = []
    count = 0

    for header, sequence in read_fasta(fasta_file):
        count += 1
        stats = SequenceStats(sequence)
        gc_contents.append(stats.gc_content())
        lengths.append(len(sequence))

    results[dataset_name] = {
        'count': count,
        'avg_length': statistics.mean(lengths),
        'avg_gc': statistics.mean(gc_contents),
        'gc_stdev': statistics.stdev(gc_contents) if len(gc_contents) > 1 else 0,
    }

# Print comparison table
print(f"\n{'Dataset':<20} {'Count':>10} {'Avg Length':>12} {'Avg GC':>10} {'GC StDev':>10}")
print("-" * 60)

for dataset, stats in results.items():
    print(f"{dataset:<20} {stats['count']:>10} "
          f"{stats['avg_length']:>12.0f} {stats['avg_gc']:>10.1f}% "
          f"{stats['gc_stdev']:>10.1f}%")

# Statistical analysis
print("\n" + "=" * 60)
print("Observations:")

if len(results) >= 2:
    gc_values = list(results.values())
    max_gc = max(gc_values, key=lambda x: x['avg_gc'])
    min_gc = min(gc_values, key=lambda x: x['avg_gc'])

    for name, dataset in results.items():
        if dataset == max_gc:
            print(f"  - {name}: Highest GC content ({dataset['avg_gc']:.1f}%)")
        if dataset == min_gc:
            print(f"  - {name}: Lowest GC content ({dataset['avg_gc']:.1f}%)")
```

---

## Command-Line Usage

You can also use these utilities from the command line by creating simple wrapper scripts:

### Example: count_seqs.py

```python
#!/usr/bin/env python3
"""Command-line tool to count sequences."""

import sys
from utilities.fasta_tools import count_sequences

if len(sys.argv) < 2:
    print("Usage: python3 count_seqs.py <fasta_file>")
    sys.exit(1)

fasta_file = sys.argv[1]
count = count_sequences(fasta_file)
print(f"{fasta_file}: {count} sequences")
```

Usage:
```bash
python3 count_seqs.py sequences.fasta
```

### Example: convert_format.py

```python
#!/usr/bin/env python3
"""Command-line tool for format conversion."""

import sys
from utilities.format_converter import auto_convert

if len(sys.argv) < 3:
    print("Usage: python3 convert_format.py <input> <output> [format]")
    print("Format options: FASTA, FASTQ, PHYLIP")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]
output_format = sys.argv[3] if len(sys.argv) > 3 else None

auto_convert(input_file, output_file, output_format)
```

Usage:
```bash
python3 convert_format.py reads.fastq reads.fasta
python3 convert_format.py aligned.fasta aligned.phy PHYLIP
```

### Example: subsample.py

```python
#!/usr/bin/env python3
"""Command-line tool for subsampling sequences."""

import sys
from utilities.subsample_sequences import subsample_fasta

if len(sys.argv) < 4:
    print("Usage: python3 subsample.py <input> <output> <n> [seed]")
    print("       python3 subsample.py <input> <output> 0.<fraction> [seed]")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    # Handle both integer and fraction
    if '.' in sys.argv[3]:
        fraction = float(sys.argv[3])
        seed = int(sys.argv[4]) if len(sys.argv) > 4 else None
        from utilities.subsample_sequences import subsample_fasta_stratified
        subsample_fasta_stratified(input_file, output_file, fraction=fraction, seed=seed)
    else:
        n = int(sys.argv[3])
        seed = int(sys.argv[4]) if len(sys.argv) > 4 else None
        subsample_fasta(input_file, output_file, n=n, seed=seed)
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
```

Usage:
```bash
python3 subsample.py large.fasta sample.fasta 100 42
python3 subsample.py large.fasta sample.fasta 0.1 42
```

---

## Tips for Using Examples

1. **Save scripts**: Save these examples as Python files and use them in your own work
2. **Modify as needed**: Adjust file paths, thresholds, and parameters for your data
3. **Add logging**: Add print statements to track progress on large datasets
4. **Error handling**: Add try/except blocks for production use
5. **Test first**: Test on small subsets before running on full datasets

## Common Errors and Solutions

| Error | Solution |
|-------|----------|
| `FileNotFoundError` | Use absolute paths or check file exists |
| `MemoryError` | Use iterators instead of loading all sequences |
| `ValueError: Fraction must be between 0 and 1` | Use `0.1` not `10` for 10% |
| `ImportError: No module named 'utilities'` | Add path: `sys.path.insert(0, '/path/to/scripts')` |

