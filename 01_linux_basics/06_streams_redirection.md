# Lesson 6: Streams and Redirection - The Power of Pipes

**Duration**: 45 minutes
**Topics**: `>`, `>>`, `|`, `<`, `tee`, stdin, stdout, stderr

---

## Introduction

This is where Linux becomes **powerful** for bioinformatics. You'll learn to:
- Redirect command output to files
- Chain commands together with pipes
- Build complex data processing pipelines

**Key Concept**: In Linux, everything is a **stream** of data. You can redirect these streams wherever you want.

---

## The Three Streams

Every Linux command has three standard streams:

1. **stdin** (Standard Input) - Channel 0
   - Where command gets its input
   - Usually keyboard, or piped data

2. **stdout** (Standard Output) - Channel 1
   - Where command sends its normal output
   - Usually your terminal screen

3. **stderr** (Standard Error) - Channel 2
   - Where command sends error messages
   - Usually your terminal screen (same as stdout)

```
           ┌─────────────────┐
  stdin ──→│                 │──→ stdout
    0      │    COMMAND      │     1
           │                 │──→ stderr
           └─────────────────┘     2
```

---

## Operator 1: `>` (Redirect Output, Overwrite)

**Purpose**: Send output to a file instead of screen. **OVERWRITES** existing file.

**Syntax**:
```bash
command > output_file.txt
```

**Examples**:

```bash
# Save directory listing to file
ls -l > directory_listing.txt

# Count sequences and save result
grep -c "^>" sequences.fasta > sequence_count.txt

# Extract sequence IDs to file
grep "^>" sequences.fasta > sequence_ids.txt
```

**DANGER - This Overwrites**:
```bash
# First run
echo "Hello" > file.txt
cat file.txt
# Output: Hello

# Second run - FILE IS REPLACED
echo "World" > file.txt
cat file.txt
# Output: World (Hello is gone!)
```

**Try It**:
```bash
# Create a file with directory listing
ls -l > my_listing.txt

# View it
cat my_listing.txt

# Overwrite it
echo "This is new content" > my_listing.txt
cat my_listing.txt
```

---

## Operator 2: `>>` (Redirect Output, Append)

**Purpose**: Add output to end of file. **APPENDS** to existing file.

**Syntax**:
```bash
command >> output_file.txt
```

**Examples**:

```bash
# First entry
echo "Sample1" >> samples.txt

# Add more entries
echo "Sample2" >> samples.txt
echo "Sample3" >> samples.txt

# View accumulated data
cat samples.txt
# Output:
# Sample1
# Sample2
# Sample3
```

**Practical Use - Combining FASTA Files**:
```bash
# Combine multiple FASTA files into one
cat sample1.fasta > all_samples.fasta
cat sample2.fasta >> all_samples.fasta
cat sample3.fasta >> all_samples.fasta

# Or more elegantly (we'll cover this later)
cat sample*.fasta > all_samples.fasta
```

**Try It**:
```bash
# Create log file
echo "Analysis started" > analysis_log.txt

# Append entries
echo "Step 1: Quality control complete" >> analysis_log.txt
echo "Step 2: Alignment complete" >> analysis_log.txt
echo "Step 3: Tree building complete" >> analysis_log.txt

# View log
cat analysis_log.txt
```

---

## Operator 3: `|` (Pipe - Chain Commands)

**Purpose**: Send output of one command as input to another command.

**This is the most powerful operator in Linux.**

**Syntax**:
```bash
command1 | command2 | command3
```

**Simple Examples**:

```bash
# Count number of directories
ls -l | grep "^d" | wc -l

# Find and count sequence IDs
grep "^>" sequences.fasta | wc -l

# Show only first 10 lines of output
ls -l | head -10

# Sort directory listing by size
ls -lh | sort -k5 -h
```

**Bioinformatics Examples**:

```bash
# Extract sequence IDs and save to file
grep "^>" sequences.fasta | sed 's/>//' | sort > sorted_ids.txt

# Count unique species in BLAST output
cut -f2 blast_results.txt | sort | uniq | wc -l

# Get top 10 longest sequences
awk '/^>/ {if (seqlen){print seqlen, header}; header=$0; seqlen=0; next}
     {seqlen+=length($0)}
     END {print seqlen, header}' sequences.fasta |
     sort -rn |
     head -10

# Find sequences with specific pattern and count them
grep -v "^>" sequences.fasta | grep "ATGATG" | wc -l
```

**Complex Pipeline Example**:
```bash
# Find all FASTA files, count sequences in each, sort by count
find . -name "*.fasta" |
while read file; do
    count=$(grep -c "^>" "$file")
    echo -e "$count\t$file"
done |
sort -rn
```

**Try It**:
```bash
# Count files in current directory
ls | wc -l

# Show only README files
ls | grep README

# Get first 5 files
ls | head -5

# Sort files and show last 5
ls | sort | tail -5

# Chain multiple operations
ls -l | grep "^-" | awk '{print $9}' | sort
```

---

## Operator 4: `<` (Redirect Input)

**Purpose**: Read input from file instead of keyboard.

**Syntax**:
```bash
command < input_file.txt
```

**Examples**:

```bash
# Send file contents to command
wc -l < sequences.fasta

# Sort file contents
sort < unsorted_data.txt > sorted_data.txt

# Count unique lines
sort < data.txt | uniq | wc -l
```

**Note**: This is less common because most commands accept filenames:
```bash
# These are equivalent:
wc -l < sequences.fasta
wc -l sequences.fasta
```

---

## Command: `tee` (Write to File AND Screen)

**Purpose**: Send output to file AND display on screen simultaneously.

**Syntax**:
```bash
command | tee output_file.txt
command | tee -a output_file.txt  # Append mode
```

**Examples**:

```bash
# Save and display directory listing
ls -l | tee directory_contents.txt

# Save log while watching progress
./long_running_script.sh | tee analysis_log.txt

# Append to existing log
echo "New entry" | tee -a logfile.txt
```

**Bioinformatics Use**:
```bash
# Run alignment and save output while watching
mafft --auto sequences.fasta | tee alignment.fasta

# Count sequences and save result
grep -c "^>" sequences.fasta | tee sequence_count.txt
# Displays: 150
# Also saves "150" to sequence_count.txt
```

**Try It**:
```bash
# List files and save to file
ls -l | tee file_list.txt
# You see the output AND it's saved to file_list.txt

# View the saved file
cat file_list.txt
```

---

## Redirecting Errors: `2>` and `2>&1`

**Purpose**: Control where error messages go.

**Syntax**:
```bash
command 2> error_log.txt           # Redirect errors to file
command > output.txt 2>&1          # Redirect output AND errors to same file
command > output.txt 2> error.txt  # Separate output and errors
command &> all_output.txt          # Redirect both (shorthand)
```

**Examples**:

```bash
# Save only errors
find / -name "*.fasta" 2> errors.txt

# Save output, discard errors
find / -name "*.fasta" > results.txt 2> /dev/null

# Save both to same file
find / -name "*.fasta" > all_output.txt 2>&1

# Save both (shorthand)
find / -name "*.fasta" &> all_output.txt
```

**Bioinformatics Use**:
```bash
# Run IQ-TREE and capture all output
iqtree -s alignment.fasta -m MFP -bb 1000 &> iqtree_log.txt

# Run BLAST and separate output from errors
blastn -query sequences.fasta -db nt -out results.txt 2> blast_errors.log
```

---

## Practical Bioinformatics Pipelines

### Pipeline 1: Process FASTA File
```bash
# Extract headers, clean them up, count unique sequences
grep "^>" sequences.fasta |
sed 's/>//' |
cut -d' ' -f1 |
sort |
uniq |
wc -l
```

### Pipeline 2: Combine Multiple Samples
```bash
# Merge all FASTA files, renaming sequences
cat sample1.fasta sample2.fasta sample3.fasta |
sed 's/>/>COMBINED_/' |
tee combined_sequences.fasta
```

### Pipeline 3: Filter by Sequence Length
```bash
# Keep only sequences longer than 500bp
awk 'BEGIN {RS=">"}
     NR>1 {
       header=$1;
       seq="";
       for(i=2; i<=NF; i++) seq=seq$i;
       if(length(seq)>500) print ">"header"\n"seq
     }' sequences.fasta > filtered.fasta
```

### Pipeline 4: Quality Control Report
```bash
# Generate QC statistics
echo "FASTA Quality Control Report" > qc_report.txt
echo "=============================" >> qc_report.txt
echo "" >> qc_report.txt

echo "Number of sequences:" >> qc_report.txt
grep -c "^>" sequences.fasta >> qc_report.txt

echo "" >> qc_report.txt
echo "Average sequence length:" >> qc_report.txt
awk '/^>/ {if(seqlen){sum+=seqlen; count++} seqlen=0; next}
     {seqlen+=length($0)}
     END {if(count>0) print sum/count}' sequences.fasta >> qc_report.txt

cat qc_report.txt
```

### Pipeline 5: Extract and Process BLAST Results
```bash
# Get top hit for each query, extract species name
cut -f1,2,3 blast_results.txt |
sort -k1,1 -k3,3rn |
awk '!seen[$1]++ {print $0}' |
cut -f2 |
sort |
uniq -c |
sort -rn > species_counts.txt
```

---

## Common Patterns

### Pattern 1: Save and Display
```bash
command | tee output.txt
```

### Pattern 2: Filter and Count
```bash
command | grep pattern | wc -l
```

### Pattern 3: Sort and Remove Duplicates
```bash
command | sort | uniq
```

### Pattern 4: Process Multiple Files
```bash
cat file1 file2 file3 | command > output.txt
```

### Pattern 5: Chain Processing
```bash
cat input.txt |
  command1 |
  command2 |
  command3 > output.txt
```

---

## Tips and Best Practices

### 1. Test Pipes Step by Step
```bash
# Build pipeline incrementally
command1
command1 | command2
command1 | command2 | command3
# Once working, add redirection
command1 | command2 | command3 > final_output.txt
```

### 2. Use Intermediate Files for Debugging
```bash
# Instead of:
cat data.txt | process1 | process2 | process3 > final.txt

# Use:
cat data.txt | process1 > step1.txt
cat step1.txt | process2 > step2.txt
cat step2.txt | process3 > final.txt
# Now you can inspect step1.txt and step2.txt
```

### 3. Always Check Before Overwriting
```bash
# DANGEROUS
important_command > important_file.txt

# SAFER - use different name
important_command > important_file_v2.txt

# SAFEST - check first
ls -lh important_file.txt
# Then decide to overwrite or append
```

### 4. Use Comments in Complex Pipelines
```bash
# Extract sequence IDs
grep "^>" sequences.fasta |
  # Remove the > character
  sed 's/>//' |
  # Get just the ID (first field)
  cut -d' ' -f1 |
  # Sort and remove duplicates
  sort | uniq > unique_ids.txt
```

---

## Common Mistakes

### Mistake 1: Redirecting to Source File
```bash
# DANGER - This will destroy your file!
sort input.txt > input.txt

# CORRECT - Use different output file
sort input.txt > input_sorted.txt

# Or use in-place tools like sponge
sort input.txt | sponge input.txt
```

### Mistake 2: Forgetting to Redirect Errors
```bash
# Errors still go to screen
command > output.txt

# CORRECT - Redirect both
command > output.txt 2>&1
# or
command &> output.txt
```

### Mistake 3: Using > Instead of >>
```bash
# First command
echo "Line 1" > data.txt

# OOPS - This overwrites!
echo "Line 2" > data.txt
# data.txt only contains "Line 2"

# CORRECT - Use >>
echo "Line 1" > data.txt
echo "Line 2" >> data.txt
# data.txt contains both lines
```

---

## Practice Exercises

### Exercise 1: Basic Redirection

```bash
# 1. List all files and save to file_list.txt
ls -lh > file_list.txt

# 2. Append current date to the file
date >> file_list.txt

# 3. View the file
cat file_list.txt
```

### Exercise 2: Pipes

```bash
# Using data/sequences.fasta:

# 1. Count number of sequences
grep -c "^>" data/sequences.fasta

# 2. Extract sequence IDs using pipe
grep "^>" data/sequences.fasta | sed 's/>//' > seq_ids.txt

# 3. Find sequences containing "Aedes" in header
grep ">" data/sequences.fasta | grep "Aedes"

# 4. Count how many sequences have "Aedes"
grep ">" data/sequences.fasta | grep -c "Aedes"
```

### Exercise 3: Complex Pipelines

```bash
# 1. Get all sequence IDs, sort, save to file
grep "^>" data/sequences.fasta | sed 's/>//' | sort > sorted_ids.txt

# 2. Count unique species (assuming species in header)
grep "^>" data/sequences.fasta | cut -d' ' -f2 | sort | uniq | wc -l

# 3. Create summary report
echo "Sequence Summary Report" > summary.txt
echo "======================" >> summary.txt
echo "Total sequences:" >> summary.txt
grep -c "^>" data/sequences.fasta >> summary.txt
cat summary.txt
```

### Exercise 4: Real Bioinformatics Task

```bash
# Process BLAST results (using data/blast_results.txt)

# 1. Get unique species from column 2
cut -f2 data/blast_results.txt | sort | uniq > unique_species.txt

# 2. Count hits per species
cut -f2 data/blast_results.txt | sort | uniq -c | sort -rn > species_hit_counts.txt

# 3. Get top 5 species
cut -f2 data/blast_results.txt | sort | uniq -c | sort -rn | head -5
```

---

## Quick Reference

```bash
>           # Redirect output (overwrite)
>>          # Redirect output (append)
|           # Pipe to next command
<           # Redirect input
2>          # Redirect errors
2>&1        # Redirect errors to stdout
&>          # Redirect both output and errors
tee         # Write to file AND stdout
tee -a      # Append to file AND stdout
```

---

## Real-World Example

Complete bioinformatics workflow using pipes:

```bash
#!/bin/bash
# Complete analysis pipeline

# 1. Quality control
echo "Starting QC..." | tee analysis.log

# 2. Combine all samples
cat sample*.fasta > all_sequences.fasta 2>> analysis.log

# 3. Count sequences
echo "Total sequences:" | tee -a analysis.log
grep -c "^>" all_sequences.fasta | tee -a analysis.log

# 4. Extract IDs
grep "^>" all_sequences.fasta |
  sed 's/>//' |
  sort > sequence_ids.txt

# 5. Calculate average length
echo "Average sequence length:" | tee -a analysis.log
awk '/^>/ {if(seqlen){sum+=seqlen; count++} seqlen=0; next}
     {seqlen+=length($0)}
     END {print sum/count}' all_sequences.fasta |
     tee -a analysis.log

# 6. Find sequences with quality issues
grep -v "^>" all_sequences.fasta |
  grep -c "N" |
  tee -a analysis.log

echo "Analysis complete!" | tee -a analysis.log
```

---

## Check Your Understanding

1. What's the difference between `>` and `>>`?
2. What does `|` do?
3. How do you save output to a file AND display it on screen?
4. What's the difference between `1>` and `2>`?
5. How do you redirect both stdout and stderr to the same file?

<details>
<summary>Click for answers</summary>

1. `>` overwrites file, `>>` appends to file
2. Pipes output from one command as input to next command
3. Use `tee`: `command | tee output.txt`
4. `1>` redirects stdout (default), `2>` redirects stderr
5. Use `&>` or `2>&1`: `command &> output.txt` or `command > output.txt 2>&1`

</details>

---

## What's Next?

**Lesson 7: Wildcards and Pattern Matching** - Use `*`, `?`, and `[]` for batch operations.

**File**: `07_wildcards_patterns.md`

---

**Master of streams and pipes? Proceed to Lesson 7!** →
