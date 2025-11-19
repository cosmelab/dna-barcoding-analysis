# Quality Control Scripts Summary

## Overview

This directory contains educational quality control scripts for DNA barcoding analysis. The scripts demonstrate core bioinformatics concepts and can be used as a complete QC pipeline for Sanger sequencing data.

**Total Code:** 1,174 lines of well-commented Python and Bash

## Contents

### Python Scripts (3 scripts, 882 lines)

#### 1. parse_ab1.py (227 lines)
**Purpose:** Parse ABI Sanger sequencer files

**What it teaches:**
- Binary file format parsing
- Using BioPython SeqIO module
- Working with metadata
- Quality score concepts
- Multiple output formats

**Key functions:**
- `parse_ab1_file()` - Read .ab1 file using BioPython
- `print_sequence_info()` - Display extracted sequences
- `print_metadata()` - Show sequencer and sample information
- `save_fasta()` - Export sequence in FASTA format
- `save_qual()` - Export quality scores in QUAL format

**Usage:**
```bash
python3 parse_ab1.py sample.ab1
python3 parse_ab1.py sample.ab1 --save-fasta output.fasta --save-qual output.qual
```

---

#### 2. trim_quality.py (298 lines)
**Purpose:** Trim low-quality bases from sequence ends

**What it teaches:**
- Quality score statistics
- Sequence trimming algorithms
- File I/O with different formats
- Command-line argument parsing
- Statistics tracking and reporting

**Key functions:**
- `trim_sequence()` - Find and extract high-quality region
- `process_sequences()` - Apply trimming to all sequences
- `save_sequences()` - Output trimmed sequences
- `print_statistics()` - Report before/after metrics

**Algorithm:**
1. For each sequence, extract quality scores
2. Find first base with quality >= threshold
3. Find last base with quality >= threshold
4. Keep only bases between these positions
5. Discard if trimmed sequence < 50 bp

**Usage:**
```bash
python3 trim_quality.py input.fastq --quality 20 --output trimmed.fastq
python3 trim_quality.py input.fasta --qual input.qual --quality 25 --format fasta
```

---

#### 3. filter_sequences.py (357 lines)
**Purpose:** Filter sequences by quality metrics

**What it teaches:**
- Object-oriented Python (classes)
- Multiple filtering criteria
- Ambiguous base detection
- Complex statistics tracking
- Informative error reporting

**Key classes:**
- `SequenceFilter` - Encapsulates filtering logic
  - `check_length()` - Verify minimum length
  - `check_quality()` - Verify minimum average quality
  - `check_ambiguity()` - Check N character percentage
  - `filter()` - Apply all criteria to a sequence
  - `get_statistics()` - Retrieve statistics

**Filtering criteria:**
1. Minimum length (default: 0 bp)
2. Minimum average quality (default: 0)
3. Maximum N percentage (default: 100%)

**Usage:**
```bash
python3 filter_sequences.py input.fastq --min-length 300 --min-quality 20 --max-n 5
```

---

### Bash Script (1 script, 292 lines)

#### batch_qc.sh (292 lines)
**Purpose:** Automate QC workflow across multiple files

**What it teaches:**
- Bash scripting fundamentals
- File operations and loops
- Error handling with `set -e` and conditionals
- Function definitions and reuse
- Command-line parsing
- Logging and reporting
- Workflow automation

**Key functions:**
- `log_message()` - Timestamped logging to console and file
- `print_header()` - Formatted section headers
- `check_script()` - Verify dependencies exist
- `convert_ab1_to_fastq()` - Parse single .ab1 file
- `trim_sequence_file()` - Trim individual file
- `filter_sequence_file()` - Filter individual file
- `process_sequence()` - Complete pipeline for one file
- `generate_report()` - Summary report generation

**Workflow:**
```
For each .ab1 file:
  1. Parse .ab1 → FASTA + QUAL
  2. Trim low-quality bases
  3. Filter by quality criteria
  4. Output final QC sequence
```

**Usage:**
```bash
bash batch_qc.sh /path/to/ab1_files
QUALITY_THRESHOLD=25 bash batch_qc.sh /path/to/ab1_files
```

---

### Documentation Files

#### README.md (9.9 KB)
Comprehensive guide covering:
- What each script does and why
- How to install and configure
- Complete workflow examples
- Quality score concepts
- Common file formats
- Troubleshooting guide
- Educational learning points
- Professional tool references

#### EXAMPLES.md (11 KB)
Practical workflow scenarios:
- Quick start examples
- DNA barcoding project setup
- Trimming impact analysis
- Filter criteria comparison
- Batch processing of large datasets
- Quality distribution analysis
- Pipeline with error logging

#### QUICKSTART.txt (7.3 KB)
Quick reference guide:
- Installation instructions
- Command syntax for all scripts
- Workflow example
- Quality score reference table
- Common filter parameters
- Troubleshooting quick reference
- File format summary

#### SUMMARY.md (this file)
Overview of all components and learning objectives

---

## Learning Path

### Beginner Level
1. **Start with QUICKSTART.txt** - Get scripts running
2. **Run parse_ab1.py** - View sequence structure
3. **Try trim_quality.py** - Understand trimming concept
4. **Explore README.md** - Learn background concepts

### Intermediate Level
1. **Read parse_ab1.py source** - Understand file parsing
2. **Read trim_quality.py source** - Learn trimming algorithm
3. **Run through EXAMPLES.md scenarios** - Practice workflows
4. **Experiment with parameters** - See what they do

### Advanced Level
1. **Study filter_sequences.py** - Object-oriented design
2. **Analyze batch_qc.sh** - Bash scripting patterns
3. **Modify scripts** - Add your own features
4. **Integrate into pipelines** - Use with other tools

---

## Concepts Demonstrated

### Python Concepts
- **Modules and imports:** BioPython, file I/O, argparse
- **Data structures:** Lists, dictionaries, tuples
- **Functions:** Parameter passing, return values, documentation
- **Classes:** Object-oriented design (SequenceFilter class)
- **File handling:** Reading/writing different formats
- **Error handling:** Try/except blocks, error messages
- **Command-line parsing:** argparse module
- **String manipulation:** Slicing, formatting, regex patterns

### Bioinformatics Concepts
- **File formats:** FASTA, QUAL, FASTQ, AB1
- **Quality scores:** PHRED scale, scoring basics
- **Quality control:** Trimming and filtering rationale
- **Sequence analysis:** Length, composition, quality metrics
- **Metadata:** Sequencer information, sample tracking
- **DNA barcoding:** Workflow for specimen identification

### Bash Concepts
- **Variables:** Setting, expanding, environment variables
- **Control flow:** for loops, if statements, conditionals
- **Functions:** Definition, parameters, return values
- **File operations:** find, mkdir, ls, grep
- **Piping and redirection:** chaining commands
- **Error handling:** set -e, return codes, conditional execution
- **Text processing:** awk, sed basics

---

## Quality Control Workflow

### Complete Pipeline
```
Raw .ab1 files
    ↓
[parse_ab1.py] → Extract sequence + quality
    ↓
FASTA + QUAL files
    ↓
[trim_quality.py] → Remove low-quality ends
    ↓
Trimmed sequences
    ↓
[filter_sequences.py] → Remove bad sequences
    ↓
QC-passed sequences
    ↓
Ready for: alignment, phylogeny, identification
```

### Decision Points in QC

**1. Quality Trimming**
- Default: Q20 (1% error rate)
- Choose based on downstream analysis
- Balance sensitivity vs. specificity

**2. Sequence Filtering**
- Minimum length: 300-400 bp (species-specific)
- Minimum quality: 20-25 (application-specific)
- Maximum Ns: 0-5% (depends on data)

**3. Parameter Selection**
- Permissive: Maximize sequences kept
- Standard: Balance quality and quantity
- Strict: Maximum confidence in data

---

## Standard Parameter Sets

### For DNA Barcoding (COI, rbcL, etc.)
```bash
# Quick QC
python3 trim_quality.py input.fastq --quality 20
python3 filter_sequences.py input.fastq --min-length 300 --min-quality 20

# Strict QC
python3 trim_quality.py input.fastq --quality 25
python3 filter_sequences.py input.fastq --min-length 400 --min-quality 25 --max-n 1
```

### For High-Throughput Identification
```bash
# Lenient QC
python3 trim_quality.py input.fastq --quality 15
python3 filter_sequences.py input.fastq --min-length 200 --min-quality 15

# Standard QC
python3 trim_quality.py input.fastq --quality 20
python3 filter_sequences.py input.fastq --min-length 300 --min-quality 20
```

### For Assembly Projects
```bash
# Strict QC
python3 trim_quality.py input.fastq --quality 30
python3 filter_sequences.py input.fastq --min-length 500 --min-quality 30 --max-n 0
```

---

## Testing & Validation

### How to Validate Scripts
```bash
# 1. Check syntax
python3 -m py_compile *.py

# 2. Test help functions
python3 parse_ab1.py -h
python3 trim_quality.py -h
python3 filter_sequences.py -h

# 3. Test with small inputs
python3 parse_ab1.py sample.ab1

# 4. Compare before/after
wc -l input.fastq output_filtered.fastq
```

### Performance Characteristics
- **parse_ab1.py:** ~10-100 ms per file (depends on file size)
- **trim_quality.py:** ~100-500 ms per file
- **filter_sequences.py:** ~100-500 ms per file
- **batch_qc.sh:** ~1-5 seconds per file (overhead)

---

## Integration with Other Tools

### Downstream Tools
- **Alignment:** BLAST, Clustal, MAFFT
- **Phylogeny:** RAxML, PAUP*, IQ-TREE
- **Identification:** BOLD, Barcode of Life
- **Assembly:** SPADES, VELVET, MIRA
- **Visualization:** MEGA, FigTree, Jalview

### Upstream Tools
- **Sequencing:** Applied Biosystems, Sanger, NGS
- **Base calling:** KB Basecaller, GATK
- **Demultiplexing:** Saute, demultiplex_fq

---

## Troubleshooting & FAQs

### Common Issues

**Q: My quality threshold doesn't seem to work**
A: Make sure your input file has quality scores. FASTA files don't have quality; you need FASTQ or FASTA+QUAL.

**Q: How do I know what threshold to use?**
A: Start with Q20 for trimming and --min-quality 20 for filtering. Adjust based on your project needs.

**Q: Can I process very large files?**
A: Yes, these scripts use generators for memory efficiency. Tested with files up to several GB.

**Q: How do I parallelize processing?**
A: Use GNU parallel or shell loops. See EXAMPLES.md for parallel processing examples.

**Q: Can I modify the scripts?**
A: Absolutely! That's part of the learning. Add features like:
- Adapter trimming
- Chimera detection
- Base composition analysis
- Custom quality cutoffs

---

## Next Steps

After mastering QC scripts:

1. **Learn alignment:** Sequence comparison and multiple alignment
2. **Learn phylogeny:** Tree building and evolutionary analysis
3. **Learn identification:** Species assignment methods
4. **Learn visualization:** Creating publication-quality figures
5. **Learn automation:** Creating production pipelines
6. **Learn best practices:** Quality assurance and validation

---

## References

- BioPython: https://biopython.org/
- FASTQ format: https://en.wikipedia.org/wiki/FASTQ_format
- PHRED quality: https://en.wikipedia.org/wiki/Phred_quality_score
- Sanger sequencing: https://en.wikipedia.org/wiki/Sanger_sequencing
- DNA barcoding: https://en.wikipedia.org/wiki/DNA_barcoding

---

## Credits

Created as educational material for DNA barcoding analysis course.

**Key Learning Objectives Achieved:**
✓ Parse binary sequencing files
✓ Calculate and interpret quality metrics
✓ Implement filtering algorithms
✓ Design automated workflows
✓ Write production-ready code
✓ Document scientific software

---

**Ready to start?** See QUICKSTART.txt for immediate next steps!
