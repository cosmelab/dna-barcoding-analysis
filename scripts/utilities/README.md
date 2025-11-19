# Bioinformatics Utility Scripts

A collection of general-purpose Python tools for DNA sequence analysis. These utilities are designed for students learning bioinformatics and can be reused in your own analysis workflows.

## Overview

This package provides five essential utility modules for working with biological sequence data:

| Module | Purpose | Key Functions |
|--------|---------|---------------|
| **fasta_tools.py** | FASTA file operations | Parse, split, merge, filter |
| **sequence_stats.py** | Sequence analysis | GC content, composition, codons |
| **reverse_complement.py** | Sequence manipulation | Reverse complement, palindromes |
| **format_converter.py** | Format conversion | FASTA ↔ FASTQ ↔ PHYLIP ↔ GenBank |
| **subsample_sequences.py** | Subsampling & filtering | Random, stratified, grouped sampling |

## Installation

The utilities are already installed in your DNA Barcoding Analysis environment. To use them:

```python
# Option 1: Import from utilities package
from utilities.fasta_tools import read_fasta, write_fasta
from utilities.sequence_stats import SequenceStats

# Option 2: Add scripts directory to Python path
import sys
sys.path.insert(0, '/path/to/scripts/utilities')
from fasta_tools import read_fasta
```

## Quick Start Examples

### 1. Reading and Parsing FASTA Files

```python
from utilities.fasta_tools import read_fasta, count_sequences

# Count sequences in a file
num_seqs = count_sequences('sequences.fasta')
print(f"File contains {num_seqs} sequences")

# Parse FASTA file and process each sequence
for header, sequence in read_fasta('sequences.fasta'):
    print(f"{header}: {len(sequence)} bp")
```

### 2. Calculating Sequence Statistics

```python
from utilities.sequence_stats import SequenceStats

# Calculate statistics for a sequence
seq = "ATCGATCGATCGATCG"
stats = SequenceStats(seq)

print(f"Length: {stats.length} bp")
print(f"GC Content: {stats.gc_content():.1f}%")
print(f"Composition: {stats.nucleotide_composition()}")

# Print formatted summary
stats.print_summary()

# Analyze entire FASTA file
from utilities.sequence_stats import analyze_sequence_file
analyze_sequence_file('sequences.fasta')
```

### 3. Working with Reverse Complements

```python
from utilities.reverse_complement import reverse_complement, is_palindrome

# Get reverse complement
original = "ATCGATCG"
rc = reverse_complement(original)
print(f"Original: {original}")
print(f"Rev Comp: {rc}")

# Verify reverse complement is same as opposite strand
assert reverse_complement(reverse_complement(original)) == original

# Check for restriction sites (palindromes)
restriction_sites = ["GAATTC", "GGATCC", "GCTAGC"]
for site in restriction_sites:
    if is_palindrome(site):
        print(f"{site} is a palindrome (restriction site)")
```

### 4. Converting Between Sequence Formats

```python
from utilities.format_converter import FormatConverter, auto_convert

# Convert FASTQ to FASTA
converter = FormatConverter()
converter.fastq_to_fasta('reads.fastq', 'reads.fasta')

# Auto-detect format and convert
auto_convert('input.fastq', 'output.fasta')

# Convert to PHYLIP for phylogenetics analysis
auto_convert('aligned.fasta', 'aligned.phy', 'PHYLIP')
```

### 5. Subsampling Large Sequence Files

```python
from utilities.subsample_sequences import SequenceSubsampler, subsample_fasta

# Random sampling - get 100 sequences from a large file
subsample_fasta('large.fasta', 'sample.fasta', n=100, seed=42)

# Stratified sampling - ensure length diversity
subsample_fasta_stratified('large.fasta', 'sample.fasta',
                           fraction=0.1, seed=42)

# Custom subsampling
subsampler = SequenceSubsampler(seed=42)
sequences = subsampler.read_fasta_with_metadata('large.fasta')

# Get 50% of sequences, removing duplicates
unique_seqs, num_dupes = subsampler.remove_duplicates(sequences)
sample = subsampler.random_sample(unique_seqs, fraction=0.5)
subsampler.write_fasta('deduped_sample.fasta', sample)
```

## Detailed Module Documentation

### fasta_tools.py

Functions for reading, writing, and manipulating FASTA format files.

**Key Functions:**

- `read_fasta(filename)` - Iterator over (header, sequence) tuples
- `write_fasta(filename, sequences)` - Write sequences to file
- `split_fasta(input_file, output_prefix, sequences_per_file)` - Split large file
- `merge_fasta(input_files, output_file)` - Combine multiple files
- `filter_fasta(input_file, output_file, min_length, max_length)` - Filter by length
- `count_sequences(filename)` - Count total sequences
- `get_sequence_by_id(filename, seq_id)` - Retrieve specific sequence

**Example: Split and Filter**

```python
from utilities.fasta_tools import split_fasta, filter_fasta

# Split large file into chunks
results = split_fasta('large_db.fasta', 'chunks', sequences_per_file=1000)
for file, count in results.items():
    print(f"{file}: {count} sequences")

# Keep only sequences longer than 100 bp
filter_fasta('sequences.fasta', 'filtered.fasta', min_length=100)
```

### sequence_stats.py

Calculate statistical properties of DNA/RNA sequences.

**Key Functions:**

- `SequenceStats.gc_content()` - Percentage of G and C nucleotides
- `SequenceStats.nucleotide_composition()` - Count of each base
- `SequenceStats.at_content()` - Percentage of A and T nucleotides
- `SequenceStats.codon_usage()` - Frequency of codons
- `SequenceStats.is_valid_dna()` / `is_valid_rna()` - Validate sequence
- `SequenceStats.summary()` - Complete statistics dictionary
- `SequenceStats.print_summary()` - Formatted output

**Example: Analyze Entire File**

```python
from utilities.sequence_stats import SequenceStats, analyze_sequence_file
from utilities.fasta_tools import read_fasta
import statistics

# Analyze all sequences and calculate summary statistics
gc_contents = []
for header, sequence in read_fasta('sequences.fasta'):
    stats = SequenceStats(sequence)
    gc_contents.append(stats.gc_content())

print(f"Average GC content: {statistics.mean(gc_contents):.1f}%")
print(f"GC content range: {min(gc_contents):.1f}% - {max(gc_contents):.1f}%")
```

### reverse_complement.py

DNA sequence manipulation functions.

**Key Functions:**

- `reverse_complement(sequence)` - Get reverse complement
- `complement(sequence)` - Get complement (3' to 5' direction)
- `reverse(sequence)` - Reverse the sequence
- `is_palindrome(sequence)` - Check if sequence equals reverse complement
- `find_palindromes(sequence, min_length)` - Find restriction sites
- `validate_sequence(sequence)` - Check for valid bases
- `process_fasta_reverse_complement()` - Generate RC for file

**Important Concepts:**

- **DNA is double-stranded**: The same gene appears on both strands
- **Reverse Complement**: The sequence from the opposite strand read 5' to 3'
- **Why it matters**: BLAST searches find matches on both strands

```
   Forward strand:     5'-ATCGATCG-3'
   Complement strand:  3'-TAGCTAGC-5'

   Reverse complement: 5'-CGATCGAT-3' (opposite strand, read backwards)
```

**Example: Find Restriction Sites**

```python
from utilities.reverse_complement import find_palindromes

# Find common restriction enzyme recognition sites
sequence = "ATGATGATGAATTCATCGATGCTAGC"
palindromes = find_palindromes(sequence, min_length=4)

for position, site in palindromes:
    print(f"Position {position}: {site} (potential restriction site)")
```

### format_converter.py

Convert between sequence formats: FASTA, FASTQ, PHYLIP, and GenBank.

**Format Handlers:**

- `FASTAHandler` - Simple sequence format
- `FASTQHandler` - Sequence with quality scores
- `PHYLIPHandler` - Aligned sequences for phylogenetics
- `FormatConverter` - Convert between formats

**Key Functions:**

- `auto_convert(input_file, output_file, output_format)` - Smart conversion
- `FormatConverter.fastq_to_fasta()` - Drop quality scores
- `FormatConverter.fasta_to_fastq()` - Add quality scores
- `FormatConverter.fasta_to_phylip()` - Format for phylogenetics
- `FormatConverter.detect_format()` - Identify format

**Example: Prepare for Phylogenetics**

```python
from utilities.format_converter import FormatConverter

converter = FormatConverter()

# Read aligned FASTA, validate, convert to PHYLIP
converter.fasta_to_phylip('aligned_sequences.fasta', 'aligned.phy')

# Now use PHYLIP file with phylogenetics software (RAxML, MrBayes, etc.)
```

### subsample_sequences.py

Create random or stratified subsamples of sequence files.

**Key Functions:**

- `random_sample()` - Select n random sequences or fraction
- `stratified_sample()` - Sample proportionally by length ranges
- `length_based_sample()` - Sample within length range
- `grouped_sample()` - Sample from custom groups
- `remove_duplicates()` - Remove identical sequences
- `subsample_by_identity()` - Create non-redundant set
- `subsample_fasta()` - Convenience function

**Example: Prepare Training/Test Data**

```python
from utilities.subsample_sequences import SequenceSubsampler

subsampler = SequenceSubsampler(seed=42)  # Set seed for reproducibility

# Read sequences
sequences = subsampler.read_fasta_with_metadata('dataset.fasta')

# Remove duplicates
unique_seqs, num_dupes = subsampler.remove_duplicates(sequences)
print(f"Removed {num_dupes} duplicate sequences")

# Create training set (80%) and test set (20%)
train_size = int(len(unique_seqs) * 0.8)
train_set = subsampler.random_sample(unique_seqs, n=train_size)
test_indices = set(id(s) for s in train_set)
test_set = [s for s in unique_seqs if id(s) not in test_indices]

# Write subsets
subsampler.write_fasta('train.fasta', train_set)
subsampler.write_fasta('test.fasta', test_set)
```

## Common Workflows

### Workflow 1: Clean and Prepare Sequences

```python
from utilities.fasta_tools import filter_fasta, count_sequences
from utilities.sequence_stats import SequenceStats
from utilities.reverse_complement import reverse_complement

# Count original sequences
original_count = count_sequences('raw_sequences.fasta')
print(f"Original: {original_count} sequences")

# Filter by length (remove very short sequences)
filter_fasta('raw_sequences.fasta', 'filtered.fasta', min_length=200)

# Check statistics
from utilities.sequence_stats import analyze_sequence_file
analyze_sequence_file('filtered.fasta')

# Verify quality
filtered_count = count_sequences('filtered.fasta')
print(f"After filtering: {filtered_count} sequences "
      f"({100*filtered_count/original_count:.1f}% retained)")
```

### Workflow 2: Create Non-Redundant Database

```python
from utilities.fasta_tools import merge_fasta, filter_fasta
from utilities.subsample_sequences import SequenceSubsampler

# Merge all sequence files
merge_fasta(['fasta_dir/file1.fasta', 'fasta_dir/file2.fasta'], 'merged.fasta')

# Remove duplicates
subsampler = SequenceSubsampler()
seqs = subsampler.read_fasta_with_metadata('merged.fasta')
unique_seqs, num_dupes = subsampler.remove_duplicates(seqs)

# Reduce redundancy (keep sequences <95% similar)
nr_seqs = subsampler.subsample_by_identity(unique_seqs, max_identity=0.95)
subsampler.write_fasta('non_redundant.fasta', nr_seqs)

print(f"Original: {len(seqs)} → Unique: {len(unique_seqs)} "
      f"→ Non-redundant: {len(nr_seqs)}")
```

### Workflow 3: Phylogenetic Analysis Prep

```python
from utilities.format_converter import auto_convert
from utilities.subsample_sequences import subsample_fasta_stratified

# Read raw sequences (could be FASTQ from sequencer)
# auto_convert('raw_reads.fastq', 'sequences.fasta')

# Create a representative sample for alignment
subsample_fasta_stratified('sequences.fasta', 'sample_stratified.fasta',
                          fraction=0.1, seed=42)

# Convert to PHYLIP for phylogenetics software
auto_convert('aligned_sample.fasta', 'aligned.phy', 'PHYLIP')

# Now use with RAxML, MrBayes, BEAST, etc.
```

### Workflow 4: Sequence Quality Control

```python
from utilities.fasta_tools import read_fasta, write_fasta
from utilities.sequence_stats import SequenceStats
from utilities.reverse_complement import validate_sequence

filtered_sequences = []

for header, sequence in read_fasta('raw_sequences.fasta'):
    stats = SequenceStats(sequence)

    # Check criteria
    if stats.length < 100:
        continue  # Too short

    if stats.gc_content() > 60 or stats.gc_content() < 40:
        continue  # Unusual GC content

    is_valid, msg = validate_sequence(sequence)
    if not is_valid:
        print(f"Skipping {header}: {msg}")
        continue

    # Check for excessive ambiguous bases
    ambiguous = stats.ambiguous_nucleotides()
    if len(ambiguous) > stats.length * 0.05:  # >5% N's
        continue

    filtered_sequences.append((header, sequence))

write_fasta('qc_passed.fasta', filtered_sequences)
print(f"Retained {len(filtered_sequences)} sequences after QC")
```

## DNA/RNA Basics for Bioinformatics

### Nucleotide Bases

| Base | Full Name | Type | DNA/RNA |
|------|-----------|------|---------|
| A | Adenine | Purine | Both |
| T | Thymine | Pyrimidine | DNA only |
| U | Uracil | Pyrimidine | RNA only |
| G | Guanine | Purine | Both |
| C | Cytosine | Pyrimidine | Both |

### IUPAC Ambiguity Codes

| Code | Meaning | Bases |
|------|---------|-------|
| R | Purine | A or G |
| Y | Pyrimidine | C or T |
| W | Weak | A or T |
| S | Strong | G or C |
| K | Keto | G or T |
| M | Amino | A or C |
| N | Any | A, T, G, C |

### Base Pairing

In double-stranded DNA:
- A pairs with T (2 hydrogen bonds)
- G pairs with C (3 hydrogen bonds)

This is why:
- `complement(ATCG)` = `TAGC`
- G-C bonds are stronger than A-T bonds
- GC content affects DNA stability

### Reverse Complement

The same genetic information appears on both DNA strands. To read the opposite strand in the 5' to 3' direction, take the reverse complement:

```
Original (5' to 3'):    5'-ATCGATCG-3'
Complement (3' to 5'):  3'-TAGCTAGC-5'
Reverse complement:     5'-CGATCGAT-3'
```

## Tips and Best Practices

1. **Use seeds for reproducibility**: When subsampling, set `seed=42` (or any number) to get consistent results
   ```python
   subsampler = SequenceSubsampler(seed=42)
   ```

2. **Check for duplicates**: Before analysis, remove duplicate sequences
   ```python
   unique_seqs, num_dupes = subsampler.remove_duplicates(sequences)
   ```

3. **Validate before analysis**: Always check sequence validity
   ```python
   stats = SequenceStats(seq)
   assert stats.is_valid_dna(), f"Invalid DNA: {seq}"
   ```

4. **Use stratified sampling**: For machine learning, use stratified sampling to preserve data distribution
   ```python
   train = subsampler.stratified_sample(seqs, fraction=0.8)
   ```

5. **Create backups**: Always keep original files
   ```bash
   cp sequences.fasta sequences.fasta.backup
   ```

6. **Log your operations**: Document what you did for reproducibility
   ```python
   import sys
   print(f"Processed {len(seqs)} sequences", file=sys.stderr)
   ```

## File Format Examples

### FASTA Example
```
>seq_001 mitochondrial DNA from species A
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
>seq_002 mitochondrial DNA from species B
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
```

### FASTQ Example
```
@seq_001 first read
ATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@seq_002 second read
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

### PHYLIP Example
```
4 20
seq_001  ATCGATCGATCGATCGATCG
seq_002  GCTAGCTAGCTAGCTAGCTA
seq_003  ATCGNCGATCGATCGATCGA
seq_004  GCTAGCTAGCTAGCTAGCTA
```

## Troubleshooting

### ImportError: No module named 'utilities'

Make sure you're in the correct directory or add the path:
```python
import sys
sys.path.insert(0, '/Users/lucianocosme/Projects/dna-barcoding-analysis/scripts')
```

### FileNotFoundError

Use absolute paths:
```python
import os
fasta_file = os.path.join(os.path.expanduser('~'), 'data', 'sequences.fasta')
```

### Memory Issues with Large Files

Use iterators instead of loading everything into memory:
```python
# Good - processes one sequence at a time
for header, seq in read_fasta('huge_file.fasta'):
    process(seq)

# Avoid - loads everything into memory
sequences = list(read_fasta('huge_file.fasta'))
```

## Contributing and Extending

These tools are designed to be extended. To add your own functions:

1. Edit the appropriate module file
2. Add your function with detailed docstrings
3. Test your code
4. Update the README with examples

Example of adding a custom function:
```python
# In sequence_stats.py
def calculate_entropy(sequence: str) -> float:
    """Calculate Shannon entropy of sequence."""
    from collections import Counter
    counts = Counter(sequence.upper())
    total = len(sequence)
    entropy = 0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    return entropy
```

## References and Further Reading

- **NCBI Bioinformatics**: https://www.ncbi.nlm.nih.gov/
- **IUPAC Nucleotide Codes**: https://www.bioinformatics.org/sms/iupac.html
- **Sequence File Formats**: http://www.sanger.ac.uk/resources/software/seqtk/
- **BioPython Documentation**: https://biopython.org/

## License

MIT License - These tools are free to use and modify for educational purposes.

## Questions?

Refer to the docstrings in each module for detailed function documentation:
```python
from utilities.fasta_tools import read_fasta
help(read_fasta)
```

Or run any module directly to see usage examples:
```bash
python3 /path/to/utilities/fasta_tools.py
```
