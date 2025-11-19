# Quality Control Scripts - Usage Examples

Quick reference guide with practical examples for each QC script.

## Quick Start

### 1. Parse AB1 Files
```bash
# View sequence information
python3 parse_ab1.py sample.ab1

# Extract sequence (FASTA)
python3 parse_ab1.py sample.ab1 --save-fasta sample.fasta

# Extract quality scores
python3 parse_ab1.py sample.ab1 --save-qual sample.qual
```

### 2. Trim Low Quality Bases
```bash
# Trim with Q20 threshold (default)
python3 trim_quality.py input.fastq --output trimmed.fastq

# Stricter trimming (Q25)
python3 trim_quality.py input.fastq --quality 25 --output trimmed.fastq
```

### 3. Filter Low Quality Sequences
```bash
# Filter with default parameters
python3 filter_sequences.py input.fastq --output filtered.fastq

# Filter with custom criteria
python3 filter_sequences.py input.fastq \
    --min-length 300 \
    --min-quality 20 \
    --max-n 5 \
    --output filtered.fastq
```

### 4. Batch QC Processing
```bash
# Process all AB1 files in directory
bash batch_qc.sh /path/to/ab1_files

# Use custom parameters
QUALITY_THRESHOLD=25 MIN_LENGTH=400 bash batch_qc.sh /path/to/ab1_files
```

---

## Scenario: DNA Barcoding Project

### Project Setup
```bash
# Create project structure
mkdir dna_barcoding_project
cd dna_barcoding_project

# Create subdirectories
mkdir raw_sequences
mkdir qc_sequences
mkdir aligned_sequences
mkdir results

# Copy your AB1 files
cp /mnt/sequencer/sample_*.ab1 raw_sequences/
```

### Step 1: Convert AB1 to FASTQ
```bash
# Parse individual AB1 files
for file in raw_sequences/*.ab1; do
    python3 ../scripts/quality_control/parse_ab1.py "$file" \
        --save-fasta "qc_sequences/$(basename $file .ab1).fasta" \
        --save-qual "qc_sequences/$(basename $file .ab1).qual"
done
```

### Step 2: Quality Trimming
```bash
# Trim all sequences
for file in qc_sequences/*.fasta; do
    python3 ../scripts/quality_control/trim_quality.py "$file" \
        --qual "${file%.fasta}.qual" \
        --quality 20 \
        --format fasta \
        --output "${file%.fasta}_trimmed.fasta"
done
```

### Step 3: Quality Filtering
```bash
# Filter sequences
for file in qc_sequences/*_trimmed.fasta; do
    python3 ../scripts/quality_control/filter_sequences.py "$file" \
        --format fasta \
        --min-length 400 \
        --min-quality 20 \
        --max-n 2 \
        --output "qc_sequences/$(basename $file _trimmed.fasta)_final.fasta"
done
```

### Step 4: Generate QC Report
```bash
# Count sequences before and after QC
echo "=== QC Summary ===" > results/qc_summary.txt
echo "Raw AB1 files: $(ls raw_sequences/*.ab1 | wc -l)" >> results/qc_summary.txt
echo "Final sequences: $(ls qc_sequences/*_final.fasta | wc -l)" >> results/qc_summary.txt
echo "" >> results/qc_summary.txt

# Show sequence statistics
echo "Sequence lengths:" >> results/qc_summary.txt
for file in qc_sequences/*_final.fasta; do
    grep -v "^>" "$file" | awk '{print length}' | sort -n | \
    awk '{count++; sum+=$1} END {print "  Min: " min ", Max: " max ", Avg: " sum/count}' >> results/qc_summary.txt
done
```

---

## Scenario: Analyzing Trimming Impact

### Before & After Comparison
```bash
# Original sequences
INPUT_FILE="sequences.fastq"

# Run trimming
python3 trim_quality.py "$INPUT_FILE" --quality 20 --output trimmed_q20.fastq
python3 trim_quality.py "$INPUT_FILE" --quality 25 --output trimmed_q25.fastq
python3 trim_quality.py "$INPUT_FILE" --quality 30 --output trimmed_q30.fastq

# Compare results
echo "=== Trimming Impact ==="
echo "Original sequences:"
grep "^@" "$INPUT_FILE" | wc -l
echo "After Q20 trimming:"
grep "^@" trimmed_q20.fastq | wc -l
echo "After Q25 trimming:"
grep "^@" trimmed_q25.fastq | wc -l
echo "After Q30 trimming:"
grep "^@" trimmed_q30.fastq | wc -l
```

---

## Scenario: Comparing Filter Criteria

### Test Different Thresholds
```bash
INPUT="trimmed.fastq"

# Test different minimum quality scores
for QUAL in 15 20 25 30; do
    echo "Testing Q${QUAL}..."
    python3 filter_sequences.py "$INPUT" \
        --min-quality $QUAL \
        --min-length 300 \
        --output "filtered_q${QUAL}.fastq"
done

# Test different minimum lengths
for LEN in 200 300 400 500; do
    echo "Testing min length ${LEN}..."
    python3 filter_sequences.py "$INPUT" \
        --min-quality 20 \
        --min-length $LEN \
        --output "filtered_len${LEN}.fastq"
done

# Test different N percentage thresholds
for MAXN in 0 2 5 10; do
    echo "Testing max ${MAXN}% Ns..."
    python3 filter_sequences.py "$INPUT" \
        --min-quality 20 \
        --max-n $MAXN \
        --output "filtered_maxn${MAXN}.fastq"
done
```

---

## Scenario: Batch Processing Many Samples

### Process Large Dataset
```bash
#!/bin/bash
# Process all AB1 files from a sequencing run

INPUT_DIR="raw_ab1_files"
OUTPUT_DIR="qc_results"

mkdir -p "$OUTPUT_DIR"

# Process each AB1 file
for ab1_file in "$INPUT_DIR"/*.ab1; do
    sample_name=$(basename "$ab1_file" .ab1)
    echo "Processing $sample_name..."

    # Parse AB1
    python3 parse_ab1.py "$ab1_file" \
        --save-fasta "$OUTPUT_DIR/${sample_name}.fasta" \
        --save-qual "$OUTPUT_DIR/${sample_name}.qual"

    # Trim
    python3 trim_quality.py "$OUTPUT_DIR/${sample_name}.fasta" \
        --qual "$OUTPUT_DIR/${sample_name}.qual" \
        --quality 20 \
        --format fasta \
        --output "$OUTPUT_DIR/${sample_name}_trimmed.fasta"

    # Filter
    python3 filter_sequences.py "$OUTPUT_DIR/${sample_name}_trimmed.fasta" \
        --format fasta \
        --min-length 300 \
        --min-quality 20 \
        --max-n 5 \
        --output "$OUTPUT_DIR/${sample_name}_final.fasta"

    # Clean up intermediate files
    rm "$OUTPUT_DIR/${sample_name}.fasta" "$OUTPUT_DIR/${sample_name}.qual"
    rm "$OUTPUT_DIR/${sample_name}_trimmed.fasta"

    echo "  Complete!"
done

echo "All samples processed."
```

---

## Scenario: Analyzing Sequence Quality Distribution

### Quality Analysis Script
```bash
#!/bin/bash
# Analyze quality distribution in FASTQ files

FASTQ_FILE="$1"

echo "=== Quality Analysis for $FASTQ_FILE ==="

# Extract quality scores (4th line of FASTQ repeats)
awk 'NR % 4 == 0' "$FASTQ_FILE" | head -1000 | while read qual; do
    # Convert ASCII to numbers (this is simplified)
    echo "$qual"
done | fold -w1 | while read char; do
    printf "%d\n" "'$char"
done | sort -n | uniq -c | sort -rn | head -20

echo ""
echo "Note: Quality scores are ASCII-encoded."
echo "! = 0, \" = 1, ... # = 2, etc. (Phred+33 encoding)"
```

---

## Scenario: Pipeline with Logging

### Automated QC with Error Tracking
```bash
#!/bin/bash
# Complete QC pipeline with error handling and logging

LOG_FILE="qc_processing.log"
ERROR_LOG="qc_errors.log"

log_message() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

process_sample() {
    local ab1_file="$1"
    local base_name=$(basename "$ab1_file" .ab1)

    log_message "Processing $base_name..."

    # Parse AB1
    if ! python3 parse_ab1.py "$ab1_file" \
        --save-fasta "${base_name}.fasta" \
        --save-qual "${base_name}.qual"; then
        echo "$base_name: Failed to parse AB1" >> "$ERROR_LOG"
        return 1
    fi

    # Trim
    if ! python3 trim_quality.py "${base_name}.fasta" \
        --qual "${base_name}.qual" \
        --quality 20 \
        --format fasta \
        --output "${base_name}_trimmed.fasta"; then
        echo "$base_name: Failed trimming" >> "$ERROR_LOG"
        return 1
    fi

    # Filter
    if ! python3 filter_sequences.py "${base_name}_trimmed.fasta" \
        --format fasta \
        --min-length 300 \
        --min-quality 20 \
        --max-n 5 \
        --output "${base_name}_final.fasta"; then
        echo "$base_name: Failed filtering" >> "$ERROR_LOG"
        return 1
    fi

    log_message "  Success!"
    return 0
}

# Initialize logs
> "$LOG_FILE"
> "$ERROR_LOG"

log_message "Starting QC pipeline"

# Process all AB1 files
for ab1_file in *.ab1; do
    process_sample "$ab1_file" || true
done

log_message "Pipeline complete"
echo "Errors logged to: $ERROR_LOG"
```

---

## Quick Reference: Command Options

### parse_ab1.py
```
--save-fasta FILE    Output FASTA file
--save-qual FILE     Output QUAL file (quality scores)
```

### trim_quality.py
```
--quality Q          Minimum quality threshold (default: 20)
--format FORMAT      Input format: fastq or fasta (default: fastq)
--qual FILE          QUAL file (for FASTA format)
--output FILE        Output filename
```

### filter_sequences.py
```
--min-length L       Minimum sequence length (default: 0)
--min-quality Q      Minimum average quality (default: 0)
--max-n N            Maximum N percentage (default: 100)
--format FORMAT      Input format: fastq or fasta (default: fastq)
--output FILE        Output filename
```

### batch_qc.sh
```
QUALITY_THRESHOLD    Environment variable for quality threshold
MIN_LENGTH           Environment variable for minimum length
MAX_N                Environment variable for maximum N percentage
```

---

## Tips & Tricks

### Monitor Progress
```bash
# Watch file sizes grow as processing completes
watch -n 2 'ls -lh output_dir/'
```

### Parallel Processing
```bash
# Use GNU parallel to process multiple files
find . -name "*.ab1" | parallel python3 parse_ab1.py {}
```

### Check Quality Before Processing
```bash
# Sample first few sequences to check quality
head -40 input.fastq | python3 filter_sequences.py /dev/stdin --format fastq
```

### Generate Statistics
```bash
# Count total bases in all FASTA files
find . -name "*.fasta" -exec cat {} \; | grep -v "^>" | wc -c
```

---

## Troubleshooting Quick Reference

| Problem | Solution |
|---------|----------|
| "No such file" | Use `ls -l filename` to verify file exists |
| "No module named Bio" | Run `pip install biopython` |
| "Permission denied" | Run `chmod +x script.py` |
| "Quality threshold has no effect" | Check input file has quality data |
| "No sequences passed filter" | Thresholds might be too strict, try defaults |
| "OutOfMemory" | Process files in smaller batches |

---

For more details, see README.md in this directory.
