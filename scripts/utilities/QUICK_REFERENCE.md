# Bioinformatics Utilities - Quick Reference Card

## Module Imports

```python
# Individual imports
from utilities.fasta_tools import read_fasta, write_fasta
from utilities.sequence_stats import SequenceStats
from utilities.reverse_complement import reverse_complement
from utilities.format_converter import auto_convert
from utilities.subsample_sequences import SequenceSubsampler

# Or from utilities package
from utilities import read_fasta, SequenceStats, reverse_complement
```

---

## FASTA Tools

### Reading FASTA Files
```python
# Single iteration (memory-efficient)
for header, sequence in read_fasta('data.fasta'):
    print(f"{header}: {len(sequence)} bp")

# Convert to list
sequences = list(read_fasta('data.fasta'))

# Count sequences
n = count_sequences('data.fasta')

# Get specific sequence
seq = get_sequence_by_id('data.fasta', 'seq_001')
```

### Writing FASTA Files
```python
# Write sequences
sequences = [('seq1', 'ATCGATCG'), ('seq2', 'GCTAGCTA')]
write_fasta('output.fasta', sequences, line_width=80)

# Write with custom line width
write_fasta('output.fasta', sequences, line_width=100)
```

### Manipulating FASTA Files
```python
# Split large file
results = split_fasta('large.fasta', 'chunk', sequences_per_file=1000)

# Merge files
merge_fasta(['file1.fasta', 'file2.fasta'], 'merged.fasta')

# Filter by length
filter_fasta('input.fasta', 'output.fasta', min_length=100, max_length=1000)

# Filter by header pattern
filter_fasta('input.fasta', 'output.fasta', header_pattern='mitochondrial')
```

---

## Sequence Statistics

### Basic Statistics
```python
seq = "ATCGATCGATCG"
stats = SequenceStats(seq)

# Basic properties
print(stats.length)                      # 12
print(stats.gc_content())                # 50.0
print(stats.at_content())                # 50.0

# Composition
comp = stats.nucleotide_composition()    # {'A': 3, 'T': 3, 'G': 3, 'C': 3, 'N': 0}
pct = stats.nucleotide_percentages()     # {'A': 25.0, 'T': 25.0, 'G': 25.0, 'C': 25.0}

# Codons
codons = stats.codon_usage()             # Codon frequencies

# Validation
is_dna = stats.is_valid_dna()            # True
is_rna = stats.is_valid_rna()            # True
```

### File Analysis
```python
# Analyze entire file
analyze_sequence_file('sequences.fasta')

# Calculate statistics
gc_contents = []
for header, seq in read_fasta('data.fasta'):
    stats = SequenceStats(seq)
    gc_contents.append(stats.gc_content())

print(f"Mean GC: {statistics.mean(gc_contents):.1f}%")
print(f"Range: {min(gc_contents):.1f}% - {max(gc_contents):.1f}%")
```

### Summary
```python
stats = SequenceStats(seq)

# Get dictionary
summary = stats.summary()
# {'length': 16, 'gc_content': 50.0, 'at_content': 50.0, ...}

# Print formatted
stats.print_summary()
```

---

## Reverse Complement

### Basic Operations
```python
# Reverse complement (most important!)
rc = reverse_complement("ATCGATCG")      # "CGATCGAT"

# Just complement
comp = complement("ATCGATCG")            # "TAGCTAGC"

# Just reverse
rev = reverse("ATCGATCG")                # "GCTAGCTA"

# Verify roundtrip
assert reverse_complement(reverse_complement(seq)) == seq
```

### Detection
```python
# Check if palindrome (restriction site)
if is_palindrome("GAATTC"):              # True (EcoRI)
    print("Restriction site!")

# Find all palindromes
pals = find_palindromes(sequence, min_length=4)
# [(position, palindrome), ...]

# Validate sequence
is_valid, msg = validate_sequence("ATCGATCG")
```

### Orientations
```python
# Get all 4 orientations
orientations = generate_all_orientations("ATCG")
# {
#   'forward': 'ATCG',
#   'reverse': 'GCTA',
#   'complement': 'TAGC',
#   'reverse_complement': 'CGAT'
# }

# Process FASTA file
process_fasta_reverse_complement('input.fasta', 'output_rc.fasta')
```

---

## Format Converter

### Auto Convert
```python
# Smart conversion (detects input format, infers output from extension)
auto_convert('reads.fastq', 'reads.fasta')
auto_convert('aligned.fasta', 'aligned.phy', 'PHYLIP')

# Specify format
auto_convert('input.fasta', 'output.fastq', 'FASTQ')
```

### FASTQ Operations
```python
# FASTQ to FASTA
converter = FormatConverter()
converter.fastq_to_fasta('reads.fastq', 'reads.fasta')

# FASTA to FASTQ (with quality scores)
converter.fasta_to_fastq('reads.fasta', 'reads.fastq', quality_score=30)
```

### PHYLIP Operations
```python
# FASTA to PHYLIP (for phylogenetics)
converter.fasta_to_phylip('aligned.fasta', 'aligned.phy')

# PHYLIP to FASTA
converter.phylip_to_fasta('aligned.phy', 'aligned.fasta')

# Read PHYLIP
num_species, num_sites, sequences = PHYLIPHandler.read('aligned.phy')
```

### Format Detection
```python
detected = FormatConverter.detect_format('unknown_file')
# Returns 'FASTA', 'FASTQ', or 'UNKNOWN'
```

---

## Subsampling

### Quick Sampling
```python
# Random sample (n sequences)
subsample_fasta('large.fasta', 'sample.fasta', n=100, seed=42)

# Fraction (10% of file)
subsample_fasta('large.fasta', 'sample.fasta', fraction=0.1, seed=42)

# Stratified (preserves length distribution)
subsample_fasta_stratified('large.fasta', 'sample.fasta', fraction=0.1, seed=42)
```

### Advanced Sampling
```python
subsampler = SequenceSubsampler(seed=42)

# Load sequences
seqs = subsampler.read_fasta_with_metadata('data.fasta')

# Random sample
sample = subsampler.random_sample(seqs, n=100)

# Stratified by length
sample = subsampler.stratified_sample(seqs, n=100, num_strata=5)

# By length range
sample = subsampler.length_based_sample(seqs, min_length=400, max_length=600, n=100)

# Custom groups
group_func = lambda h: h.split('_')[0]  # Group by first part
sample = subsampler.grouped_sample(seqs, group_func, n=100)

# Write results
subsampler.write_fasta('output.fasta', sample)
```

### Deduplication
```python
subsampler = SequenceSubsampler()
seqs = subsampler.read_fasta_with_metadata('data.fasta')

# Remove duplicates
unique_seqs, n_dupes = subsampler.remove_duplicates(seqs)

# Reduce by identity (non-redundant)
nr_seqs = subsampler.subsample_by_identity(seqs, max_identity=0.95)

# Statistics
stats = subsampler.get_sequence_statistics(seqs)
# {'count': 1000, 'total_length': 500000, 'min_length': 400, ...}
```

---

## Common Workflows

### Quality Control
```python
from utilities.fasta_tools import read_fasta, write_fasta
from utilities.sequence_stats import SequenceStats

qc_passed = []
for header, seq in read_fasta('raw.fasta'):
    stats = SequenceStats(seq)

    if stats.length < 100: continue       # Too short
    if stats.gc_content() > 60: continue  # Unusual GC
    if not stats.is_valid_dna(): continue # Invalid

    qc_passed.append((header, seq))

write_fasta('qc_passed.fasta', qc_passed)
```

### Train/Test Split
```python
from utilities.subsample_sequences import SequenceSubsampler

subsampler = SequenceSubsampler(seed=42)
seqs = subsampler.read_fasta_with_metadata('all.fasta')

# Remove duplicates
unique, _ = subsampler.remove_duplicates(seqs)

# Split
train = subsampler.random_sample(unique, fraction=0.8)
test_ids = set(id(s) for s in train)
test = [s for s in unique if id(s) not in test_ids]

subsampler.write_fasta('train.fasta', train)
subsampler.write_fasta('test.fasta', test)
```

### Prepare for Phylogenetics
```python
from utilities.fasta_tools import filter_fasta
from utilities.format_converter import auto_convert

# Filter similar lengths
filter_fasta('seqs.fasta', 'filtered.fasta', min_length=600, max_length=700)

# Convert to PHYLIP
auto_convert('filtered.fasta', 'aligned.phy', 'PHYLIP')

# Use with RAxML, MrBayes, IQTree, etc.
```

### Build Reference Database
```python
from utilities.fasta_tools import merge_fasta
from utilities.subsample_sequences import SequenceSubsampler

# Merge sources
merge_fasta(['source1.fasta', 'source2.fasta'], 'merged.fasta')

# Deduplicate and reduce redundancy
subsampler = SequenceSubsampler()
seqs = subsampler.read_fasta_with_metadata('merged.fasta')
unique, _ = subsampler.remove_duplicates(seqs)
nr_seqs = subsampler.subsample_by_identity(unique, max_identity=0.95)

# Export
subsampler.write_fasta('reference_nr.fasta', nr_seqs)
```

---

## Common Issues and Solutions

| Problem | Solution |
|---------|----------|
| `ImportError: No module named 'utilities'` | Add path: `sys.path.insert(0, '/path/to/scripts')` |
| `FileNotFoundError` | Use absolute paths or verify file exists |
| Out of memory | Use iterators: `for h, s in read_fasta(file)` not `list(read_fasta(file))` |
| Fraction not working | Use `0.1` not `10` for 10% |
| Non-reproducible results | Set seed: `SequenceSubsampler(seed=42)` |
| Alignment fails | Check all sequences have same length in PHYLIP |

---

## Function Quick Index

### fasta_tools
- `read_fasta(f)` - Read FASTA file
- `write_fasta(f, seqs)` - Write FASTA file
- `count_sequences(f)` - Count sequences
- `get_sequence_by_id(f, id)` - Get one sequence
- `split_fasta(in, prefix, n)` - Split file
- `merge_fasta(files, out)` - Merge files
- `filter_fasta(in, out, min, max)` - Filter sequences

### sequence_stats
- `SequenceStats(seq)` - Create stats object
- `.gc_content()` - Get GC%
- `.nucleotide_composition()` - Count bases
- `.is_valid_dna()` - Validate
- `.summary()` - Get stats dict
- `.print_summary()` - Print formatted

### reverse_complement
- `reverse_complement(seq)` - RC (most used!)
- `complement(seq)` - Get complement
- `reverse(seq)` - Reverse string
- `is_palindrome(seq)` - Check palindrome
- `find_palindromes(seq)` - Find restriction sites
- `validate_sequence(seq)` - Check validity

### format_converter
- `auto_convert(in, out, fmt)` - Convert formats
- `FASTQHandler.read(f)` - Read FASTQ
- `PHYLIPHandler.write(f, seqs)` - Write PHYLIP
- `FormatConverter.detect_format(f)` - Detect format

### subsample_sequences
- `subsample_fasta(in, out, n)` - Quick sample
- `SequenceSubsampler()` - Create sampler
- `.random_sample(seqs, n)` - Random sample
- `.stratified_sample(seqs, n)` - Stratified sample
- `.remove_duplicates(seqs)` - Remove dupes
- `.read_fasta_with_metadata(f)` - Read with stats

---

## Key Concepts

### Reverse Complement
```
Original 5'-3':    ATCGATCG
Opposite strand 3'-5': TAGCTAGC
Reverse complement 5'-3': CGATCGAT
```
Use when: Searching for sequences (BLAST), building databases

### GC Content
- % of G and C nucleotides
- Ranges 0-100%
- Affects: DNA stability, melting temperature, complexity
- Normal range: 30-70%

### Stratified Sampling
- Sample proportionally from each group
- Preserves distribution (e.g., length ranges)
- Better for machine learning than random sampling

### IUPAC Codes
- N = Any (ATGC)
- R = puRine (AG)
- Y = pYrimidine (CT)
- W = Weak (AT, 2 H-bonds)
- S = Strong (GC, 3 H-bonds)

---

## Tips

1. Always use **absolute paths**: `/path/to/file.fasta` not `./file.fasta`
2. Set **seed** for reproducibility: `SequenceSubsampler(seed=42)`
3. Use **iterators** for large files: `read_fasta()` not list
4. **Test on small files** first before processing large datasets
5. **Check line endings**: Some utilities prefer Unix (LF) over Windows (CRLF)
6. **Verify output**: Always spot-check first few sequences
7. **Keep backups**: Never overwrite original files
8. **Use try/except**: Add error handling for production use

---

## Learning Path

1. Start with `fasta_tools` - basic file operations
2. Learn `sequence_stats` - understand sequence properties
3. Master `reverse_complement` - grasp DNA biology
4. Use `format_converter` - work with different data types
5. Advanced `subsample_sequences` - prepare data for analysis

---

## Resources

- **Full documentation**: See README.md
- **Code examples**: See EXAMPLES.md
- **Docstrings**: `help(function_name)` in Python
- **Inline comments**: Read source code
- **Module info**: Run module directly: `python3 fasta_tools.py`

---

**Remember**: These utilities are designed for learning AND production use. Use them, modify them, improve them!
