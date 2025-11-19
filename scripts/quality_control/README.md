# DNA Sequence Quality Control Scripts

Educational examples of quality control workflows for DNA barcoding sequences. These scripts are designed to be easy to understand and learn from, with extensive comments explaining each step.

## Scripts Overview

### 1. `parse_ab1.py` - Parse ABI Sequencing Files

**What it does:**
- Reads .ab1 binary files from DNA sequencers
- Extracts DNA sequences and quality scores
- Displays metadata (machine info, sample name, etc.)
- Saves sequence in FASTA format
- Saves quality scores in QUAL format

**What's an AB1 file?**
- Binary file format used by Applied Biosystems DNA sequencers
- Contains: basecalls, quality scores, trace data, metadata
- One file = one DNA sequence reaction (one sample)

**Usage examples:**
```bash
# View sequence information from an AB1 file
python3 parse_ab1.py sample.ab1

# Extract sequence to FASTA
python3 parse_ab1.py sample.ab1 --save-fasta output.fasta

# Extract quality scores to QUAL format
python3 parse_ab1.py sample.ab1 --save-qual output.qual

# Save both
python3 parse_ab1.py sample.ab1 --save-fasta seq.fasta --save-qual seq.qual
```

**Learning points:**
- How to work with binary file formats
- BioPython SeqIO for parsing different sequence formats
- Accessing metadata from sequence records
- Working with quality scores (PHRED scale)

---

### 2. `trim_quality.py` - Trim Low-Quality Bases

**What it does:**
- Removes low-quality bases from sequence ends
- Finds first/last position above quality threshold
- Removes sequences that become too short after trimming
- Reports trimming statistics

**Why trim?**
- DNA sequencers typically have worse accuracy at read ends
- Low-quality bases cause errors in downstream analysis
- Trimming improves alignment, assembly, and identification accuracy

**Quality Score Scale (PHRED):**
- Q20 = 1% error rate (99% accuracy) ← common minimum
- Q25 = 0.3% error rate
- Q30 = 0.1% error rate (very high quality)

**Usage examples:**
```bash
# Trim with default Q20 threshold
python3 trim_quality.py input.fastq --output trimmed.fastq

# Stricter trimming with Q25
python3 trim_quality.py input.fastq --quality 25 --output trimmed.fastq

# Trim FASTA file with separate QUAL file
python3 trim_quality.py input.fasta --qual input.qual --format fasta

# Very strict trimming (Q30, minimum 400 bp)
python3 trim_quality.py input.fastq --quality 30 --min-length 400
```

**Output statistics:**
- Total/passed/failed sequences
- Bases before and after trimming
- Average sequence length
- Quality distribution

**Learning points:**
- How to calculate quality statistics
- Sequence slicing and manipulation
- Handling file I/O in Python
- Writing informative error messages

---

### 3. `filter_sequences.py` - Filter by Quality Metrics

**What it does:**
- Removes entire sequences that don't meet quality standards
- Filters by: minimum length, average quality, ambiguous bases (N)
- Keeps track of why sequences were filtered out
- Reports detailed filtering statistics

**Different from trimming:**
- **Trimming:** Removes low-quality bases from ends
- **Filtering:** Removes entire bad sequences

**Common filtering parameters:**
- Minimum length: 300-400 bp
- Minimum average quality: Q20-Q30
- Maximum N characters: 0-5%

**Usage examples:**
```bash
# Basic filtering with defaults
python3 filter_sequences.py input.fastq --output filtered.fastq

# Filter with specific criteria
python3 filter_sequences.py input.fastq \
    --min-length 300 \
    --min-quality 20 \
    --max-n 5 \
    --output filtered.fastq

# Strict filtering for high-quality assemblies
python3 filter_sequences.py input.fastq \
    --min-length 400 \
    --min-quality 25 \
    --max-n 1 \
    --output filtered.fastq

# Filter FASTA files
python3 filter_sequences.py input.fasta \
    --format fasta \
    --min-length 300
```

**Output statistics:**
- Sequences passed/failed by each criterion
- Percentage of sequences retained
- Total bases before/after filtering
- Failure breakdown by reason

**Learning points:**
- Object-oriented programming with Python classes
- Multiple filtering criteria
- Statistics tracking and reporting
- Handling edge cases (no quality scores, empty sequences)

---

### 4. `batch_qc.sh` - Batch Processing Workflow

**What it does:**
- Automates QC workflow across multiple files
- Converts AB1 → trimmed → filtered
- Generates summary reports
- Logs all processing steps

**Workflow steps:**
1. Find all .ab1 files in input directory
2. For each file:
   - Parse AB1 and extract FASTA/QUAL
   - Trim low-quality bases
   - Filter by quality criteria
3. Generate summary report

**Usage examples:**
```bash
# Process AB1 files in current directory
bash batch_qc.sh

# Process AB1 files in specific directory
bash batch_qc.sh /path/to/ab1_files

# Use custom quality thresholds
QUALITY_THRESHOLD=25 MIN_LENGTH=400 bash batch_qc.sh

# Set multiple parameters
export QUALITY_THRESHOLD=25
export MIN_LENGTH=400
export MAX_N=1
bash batch_qc.sh /data/sequences
```

**Output structure:**
```
qc_results/
├── sample1_trimmed.fasta
├── sample1_qc.fasta
├── sample2_trimmed.fasta
├── sample2_qc.fasta
└── qc_log.txt          # Detailed processing log
```

**Learning points:**
- Bash scripting fundamentals
- Loops and conditional statements
- File operations in bash
- Error handling and exit codes
- Logging and reporting
- Function definitions
- Command-line argument handling

---

## Installation & Requirements

### Requirements
- Python 3.6+
- BioPython library

### Install BioPython
```bash
pip install biopython
```

### Make scripts executable
```bash
chmod +x parse_ab1.py
chmod +x trim_quality.py
chmod +x filter_sequences.py
chmod +x batch_qc.sh
```

---

## Complete Workflow Example

Process a set of AB1 files from Sanger sequencing:

```bash
# 1. Create directory for raw AB1 files
mkdir raw_sequences
cd raw_sequences

# 2. Copy your AB1 files here
cp /path/to/*.ab1 .

# 3. Run batch QC with custom parameters
cd ..
QUALITY_THRESHOLD=20 MIN_LENGTH=300 MAX_N=5 bash scripts/quality_control/batch_qc.sh raw_sequences

# 4. Check results
ls -lh raw_sequences/qc_results/
cat raw_sequences/qc_results/qc_log.txt

# 5. Use QC-passed sequences for downstream analysis
# (alignment, phylogeny, species identification, etc.)
```

---

## Understanding Quality Scores

**PHRED Quality Scale** (standard in bioinformatics):
```
Q Score | Error Rate | Accuracy
--------|-----------|----------
Q10     | 10%       | 90%
Q20     | 1%        | 99%        ← minimum for most applications
Q30     | 0.1%      | 99.9%      ← high quality
Q40     | 0.01%     | 99.99%     ← very high (rare)
```

**Where do quality scores come from?**
- DNA sequencers estimate base-calling confidence
- Based on signal strength and peak resolution
- Reported as PHRED quality scores
- Different sequencing methods have different typical ranges

**Quality issues in Sanger sequencing:**
- Weak signal at sequence start (ramp-up)
- Weak signal at sequence end (signal decay)
- Template secondary structure
- Non-template peaks (secondary sequences)

---

## Common Bioinformatics Concepts

### Sequence Formats

**FASTA** - Sequence only
```
>sequence_id description
ATCGATCGATCGATCG...
```

**QUAL** - Quality scores
```
>sequence_id description
30 28 32 31 29 30...
```

**FASTQ** - Sequence + quality (combined)
```
@sequence_id description
ATCGATCGATCGATCG...
+
IIIIIIIIII9A)3/5
```

### Quality Encoding

FASTQ files use ASCII characters to represent quality:
```
ASCII 33:  ! (Q0)
ASCII 48:  0 (Q15)
ASCII 73:  I (Q40)
```

This allows compact storage of quality information.

---

## Troubleshooting

### "No module named Bio"
```bash
# Install BioPython
pip install biopython
```

### "File not found"
- Check file path is correct
- Check file actually exists: `ls -l filename`
- Use absolute paths if having issues

### "Permission denied" on .sh file
```bash
# Make shell script executable
chmod +x batch_qc.sh
```

### Quality threshold not making a difference
- Check input file actually has quality scores
- Verify PHRED scale (not Solexa or other)
- Try a more reasonable threshold (try Q20 first)

### No sequences passed filtering
- Thresholds might be too strict
- Start with defaults: min-length=300, min-quality=20, max-n=5
- Check input file quality distribution first

---

## Educational Notes for Students

These scripts are designed to teach bioinformatics concepts. Each one demonstrates:

1. **parse_ab1.py** - File I/O, data structures, metadata handling
2. **trim_quality.py** - Quality statistics, sequence manipulation, filtering algorithms
3. **filter_sequences.py** - Object-oriented programming, class design, multiple criteria
4. **batch_qc.sh** - Shell scripting, automation, workflow design

As you learn more bioinformatics, you'll recognize these patterns in professional tools:
- **fastqc** - Quality assessment (similar to our analysis)
- **trim_galore** - Quality trimming (similar to our trimming)
- **fastp** - Combined QC pipeline (like our batch_qc)
- **QIIME2** - Complete microbiome QC pipeline (uses similar concepts)

---

## References

- BioPython Tutorial: https://biopython.org/wiki/Documentation
- FASTQ Format Specification: https://en.wikipedia.org/wiki/FASTQ_format
- PHRED Quality Scores: https://en.wikipedia.org/wiki/Phred_quality_score
- Sanger Sequencing: https://en.wikipedia.org/wiki/Sanger_sequencing

---

## License

Educational use - feel free to modify and learn from these scripts!
