# Utilities Directory Index

## Quick Navigation

### Python Modules (Core Tools)
1. **fasta_tools.py** (307 lines)
   - Read, write, split, merge, filter FASTA files
   - Memory-efficient iterators
   - Sequence retrieval and counting

2. **sequence_stats.py** (359 lines)
   - GC content, composition, statistics
   - Sequence validation
   - Codon analysis

3. **reverse_complement.py** (392 lines)
   - Reverse, complement, reverse-complement operations
   - Palindrome detection (restriction sites)
   - Sequence validation

4. **format_converter.py** (508 lines)
   - FASTA ↔ FASTQ ↔ PHYLIP conversion
   - Format detection
   - Quality score handling

5. **subsample_sequences.py** (506 lines)
   - Random and stratified sampling
   - Deduplication
   - Identity-based clustering

6. **__init__.py** (102 lines)
   - Package initialization
   - Convenient imports

### Documentation Files
1. **README.md** (Full manual with everything)
   - Overview of all modules
   - Detailed function documentation
   - Common workflows
   - File format specs
   - Troubleshooting

2. **EXAMPLES.md** (Practical examples)
   - 30+ code examples
   - Real-world workflows
   - Command-line tools
   - Complete pipelines

3. **QUICK_REFERENCE.md** (Cheat sheet)
   - Function quick reference
   - Common code snippets
   - Tips and tricks
   - Quick index

4. **SUMMARY.md** (Overview)
   - Project summary
   - Feature highlights
   - Integration with course
   - Statistics

5. **INDEX.md** (This file)
   - Directory navigation
   - File organization

## By Use Case

### Reading/Writing Sequences
→ Use **fasta_tools.py**
- `read_fasta()` - Parse files
- `write_fasta()` - Save sequences
- `count_sequences()` - Get total count

### Analyzing Sequences
→ Use **sequence_stats.py**
- `SequenceStats.gc_content()` - Get GC%
- `SequenceStats.nucleotide_composition()` - Count bases
- `analyze_sequence_file()` - Full file analysis

### DNA Manipulation
→ Use **reverse_complement.py**
- `reverse_complement()` - Get RC (most important!)
- `is_palindrome()` - Find restriction sites
- `validate_sequence()` - Check validity

### Converting Formats
→ Use **format_converter.py**
- `auto_convert()` - Smart conversion
- `FASTQHandler` - Work with FASTQ
- `PHYLIPHandler` - Prepare for phylogenetics

### Subsampling Data
→ Use **subsample_sequences.py**
- `subsample_fasta()` - Quick sampling
- `SequenceSubsampler` - Advanced sampling
- `remove_duplicates()` - Deduplication

## By Module

### fasta_tools.py

**Class**: None (all functions)

**Main Functions**:
```
read_fasta(filename)                    - Iterator over (header, seq) tuples
write_fasta(filename, sequences)        - Write FASTA file
split_fasta(input, prefix, n)           - Split into chunks
merge_fasta(inputs, output)             - Combine files
filter_fasta(input, output, ...)        - Filter by length/pattern
count_sequences(filename)               - Count total sequences
get_sequence_by_id(filename, id)        - Get specific sequence
```

**Examples**: See README.md § FASTA Tools

---

### sequence_stats.py

**Main Class**: SequenceStats(sequence)

**Key Methods**:
```
.gc_content()                   - GC percentage (0-100)
.nucleotide_composition()       - Count each base
.nucleotide_percentages()       - Percentage of each base
.at_content()                   - AT percentage
.codon_usage()                  - Codon frequencies
.is_valid_dna()                 - Check validity
.is_valid_rna()                 - Check if RNA
.ambiguous_nucleotides()        - Find N positions
.summary()                      - Complete stats dict
.print_summary()                - Formatted output
```

**Functions**:
```
analyze_sequence_file(filename)  - Analyze all sequences in file
get_sequence_statistics(seqs)    - Stats for sequence list
```

**Examples**: See README.md § Sequence Statistics

---

### reverse_complement.py

**Functions**:
```
reverse_complement(seq)                 - Get reverse complement
complement(seq)                         - Get complement only
reverse(seq)                            - Reverse string
is_palindrome(seq)                      - Check if RC == seq
find_palindromes(seq, min_length)       - Find restriction sites
validate_sequence(seq)                  - Check if valid DNA/RNA
rotate_sequence(seq, positions)         - Rotate sequence
generate_all_orientations(seq)          - Get all 4 forms
process_fasta_reverse_complement(in,out)- Process entire file
```

**Key Constants**:
```
COMPLEMENT_DNA      - A↔T, G↔C mapping
COMPLEMENT_RNA      - A↔U, G↔C mapping
IUPAC_CODES         - Ambiguity code meanings
```

**Examples**: See README.md § Reverse Complement

---

### format_converter.py

**Classes**:
```
Sequence                - Data class (header, sequence, quality, features)
FASTAHandler           - Read/write FASTA files
FASTQHandler           - Read/write FASTQ with quality scores
PHYLIPHandler          - Read/write aligned PHYLIP format
FormatConverter        - High-level conversion functions
```

**Main Functions**:
```
auto_convert(input, output, format)     - Smart conversion
FASTAHandler.read(filename)             - Parse FASTA
FASTAHandler.write(filename, seqs)      - Write FASTA
FASTQHandler.read(filename)             - Parse FASTQ
FASTQHandler.write(filename, seqs)      - Write FASTQ
PHYLIPHandler.read(filename)            - Parse PHYLIP
PHYLIPHandler.write(filename, seqs)     - Write PHYLIP
FormatConverter.detect_format(file)     - Identify format
```

**Examples**: See README.md § Format Conversion

---

### subsample_sequences.py

**Main Class**: SequenceSubsampler(seed=None)

**Key Methods**:
```
.read_fasta_with_metadata(file)         - Load with lengths
.write_fasta(file, seqs)                - Save sequences
.random_sample(seqs, n, fraction)       - Random selection
.stratified_sample(seqs, n, num_strata) - Stratified sampling
.length_based_sample(seqs, min, max)    - Filter by length
.grouped_sample(seqs, group_func, n)    - Custom grouping
.remove_duplicates(seqs)                - Deduplicate
.subsample_by_identity(seqs, threshold) - Reduce redundancy
```

**Functions**:
```
subsample_fasta(input, output, n, seed)                 - Quick sample
subsample_fasta_stratified(input, output, fraction)     - Stratified sample
get_sequence_statistics(seqs)                           - Summary stats
```

**Examples**: See README.md § Subsampling

---

## Statistics

```
Total Lines of Code:        2,174
Total Documentation:        4,000+ words
Example Code Snippets:      30+
Functions:                  40+
Classes:                    5
Python Modules:             5
Documentation Files:        5
```

## Import Methods

### Method 1: Direct Imports
```python
import sys
sys.path.insert(0, '/Users/lucianocosme/Projects/dna-barcoding-analysis/scripts')
from utilities.fasta_tools import read_fasta
```

### Method 2: Package Imports
```python
from utilities import read_fasta, SequenceStats
```

### Method 3: Run as Scripts
```bash
python3 /path/to/utilities/fasta_tools.py
python3 /path/to/utilities/sequence_stats.py
```

## Documentation Reading Order

For **beginners**:
1. QUICK_REFERENCE.md - Get overview
2. README.md - Learn module by module
3. EXAMPLES.md - See practical code
4. Source code - Study implementation

For **experienced programmers**:
1. QUICK_REFERENCE.md - Quick lookup
2. Source code - Understand implementation
3. Docstrings - Function details
4. EXAMPLES.md - Integration patterns

## File Organization

```
/scripts/utilities/
├── Python Modules
│   ├── fasta_tools.py           (FASTA operations)
│   ├── sequence_stats.py        (Sequence analysis)
│   ├── reverse_complement.py    (DNA manipulation)
│   ├── format_converter.py      (Format conversion)
│   ├── subsample_sequences.py   (Subsampling)
│   └── __init__.py              (Package init)
│
├── Documentation
│   ├── README.md                (Full manual)
│   ├── EXAMPLES.md              (Code examples)
│   ├── QUICK_REFERENCE.md       (Cheat sheet)
│   ├── SUMMARY.md               (Overview)
│   └── INDEX.md                 (This file)
│
└── Auto-generated
    └── __pycache__/             (Compiled Python)
```

## Quality Assurance

All modules:
- ✓ Python 3.6+ compatible
- ✓ Comprehensive docstrings
- ✓ Type hints provided
- ✓ Error handling included
- ✓ Syntax validated
- ✓ Functions tested
- ✓ Memory efficient
- ✓ Production ready

## Typical Workflows

1. **Load sequences**: fasta_tools.read_fasta()
2. **Analyze**: sequence_stats.SequenceStats()
3. **Manipulate**: reverse_complement operations
4. **Convert**: format_converter.auto_convert()
5. **Sample**: subsample_sequences methods

## Getting Help

Inside Python:
```python
help(read_fasta)                    # Function help
help(SequenceStats)                 # Class help
help(SequenceStats.gc_content)      # Method help
```

Command line:
```bash
python3 -c "from utilities.fasta_tools import read_fasta; help(read_fasta)"
python3 /path/to/utilities/fasta_tools.py  # Shows usage
```

In editor:
```python
# Docstrings appear in most editors
# Ctrl+Click to jump to definition
# Hover to see function signature
```

## Integration Points

These utilities work well with:
- **BioPython** - Advanced sequence analysis
- **Pandas** - Sequence data as DataFrames
- **NumPy** - Numerical analysis
- **Matplotlib** - Visualization
- **Alignment tools** - MAFFT, Clustal, MUSCLE
- **Phylogenetics** - RAxML, MrBayes, IQTree

## Next Steps After Using Utilities

1. **Extend** - Add custom functions
2. **Optimize** - Improve for your use case
3. **Integrate** - Combine with other tools
4. **Deploy** - Use in pipelines
5. **Share** - Contribute improvements

---

**Last Updated**: November 18, 2024
**Version**: 1.0.0
**Status**: Production Ready
