# Linux Command Line Cheatsheet for Bioinformatics

Quick reference guide for essential Linux commands.

---

## Navigation

| Command | Description | Example |
|---------|-------------|---------|
| `pwd` | Print working directory | `pwd` |
| `cd directory` | Change directory | `cd 01_linux_basics` |
| `cd ..` | Go up one level | `cd ..` |
| `cd ~` | Go to home directory | `cd ~` |
| `cd -` | Go to previous directory | `cd -` |
| `ls` | List files | `ls` |
| `ls -l` | Long format listing | `ls -l` |
| `ls -a` | Show hidden files | `ls -a` |
| `ls -lh` | Human-readable sizes | `ls -lh` |
| `ls -lt` | Sort by time | `ls -lt` |
| `tree -L 2` | Show directory tree (2 levels) | `tree -L 2` |

---

## File Manipulation

| Command | Description | Example |
|---------|-------------|---------|
| `touch file.txt` | Create empty file | `touch new_file.txt` |
| `mkdir dir` | Make directory | `mkdir new_folder` |
| `mkdir -p a/b/c` | Make nested directories | `mkdir -p data/results/figures` |
| `cp source dest` | Copy file | `cp file1.txt file2.txt` |
| `cp -r dir1 dir2` | Copy directory recursively | `cp -r folder1 folder2` |
| `mv source dest` | Move or rename | `mv old.txt new.txt` |
| `rm file` | Remove file | `rm unwanted.txt` |
| `rm -r dir` | Remove directory recursively | `rm -r old_folder` |
| `rm -i file` | Remove with confirmation | `rm -i important.txt` |

---

## Viewing Files

| Command | Description | Example |
|---------|-------------|---------|
| `cat file` | Display entire file | `cat sequences.fasta` |
| `cat file1 file2` | Concatenate files | `cat *.fasta > all.fasta` |
| `less file` | Page through file | `less large_file.txt` |
| `head file` | First 10 lines | `head sequences.fasta` |
| `head -n 20 file` | First 20 lines | `head -n 20 sequences.fasta` |
| `tail file` | Last 10 lines | `tail sequences.fasta` |
| `tail -n 50 file` | Last 50 lines | `tail -n 50 logfile.txt` |
| `wc file` | Count lines/words/chars | `wc sequences.fasta` |
| `wc -l file` | Count lines only | `wc -l sequences.fasta` |

---

## Searching

| Command | Description | Example |
|---------|-------------|---------|
| `grep pattern file` | Search for pattern | `grep "Aedes" sequences.fasta` |
| `grep -i pattern file` | Case-insensitive | `grep -i "aedes" sequences.fasta` |
| `grep -v pattern file` | Invert match (exclude) | `grep -v "^>" sequences.fasta` |
| `grep -c pattern file` | Count matches | `grep -c "^>" sequences.fasta` |
| `grep -n pattern file` | Show line numbers | `grep -n "ATGATG" sequences.fasta` |
| `grep -r pattern dir` | Recursive search | `grep -r "TODO" .` |
| `find . -name "*.fasta"` | Find files by name | `find . -name "*.fasta"` |
| `find . -type f` | Find files only | `find . -type f` |
| `find . -type d` | Find directories only | `find . -type d` |

---

## Text Processing

### cut - Extract Columns
```bash
cut -d',' -f1 file.csv           # Extract column 1 (comma-delimited)
cut -d$'\t' -f2,5 file.tsv       # Extract columns 2 and 5 (tab-delimited)
cut -c1-10 file.txt              # Extract characters 1-10
```

### sort - Sort Lines
```bash
sort file.txt                    # Sort alphabetically
sort -r file.txt                 # Reverse sort
sort -n file.txt                 # Numeric sort
sort -k2 file.txt                # Sort by column 2
sort -t',' -k3 -n file.csv       # Sort CSV by column 3 numerically
```

### uniq - Remove Duplicates
```bash
sort file.txt | uniq             # Remove adjacent duplicates (must sort first!)
sort file.txt | uniq -c          # Count occurrences
sort file.txt | uniq -d          # Show only duplicates
```

### paste - Merge Files Horizontally
```bash
paste file1.txt file2.txt        # Merge side-by-side
paste -d',' file1.txt file2.txt  # Use comma delimiter
```

### tr - Translate Characters
```bash
tr 'a-z' 'A-Z' < file.txt        # Convert to uppercase
tr -d ' ' < file.txt             # Delete spaces
tr '\n' ',' < file.txt           # Replace newlines with commas
```

### sed - Stream Editor
```bash
sed 's/old/new/' file.txt        # Replace first occurrence per line
sed 's/old/new/g' file.txt       # Replace all occurrences
sed 's/^>//' file.fasta          # Remove > from start of lines
sed -n '5,10p' file.txt          # Print lines 5-10
sed '/^#/d' file.txt             # Delete comment lines
```

### awk - Pattern Processing
```bash
awk '{print $1}' file.txt        # Print first column
awk '{print $1,$3}' file.txt     # Print columns 1 and 3
awk -F',' '{print $2}' file.csv  # Use comma as delimiter
awk '$3 > 100' file.txt          # Print lines where column 3 > 100
awk '{sum+=$1} END {print sum}'  # Sum column 1
```

---

## Streams and Redirection

| Operator | Description | Example |
|----------|-------------|---------|
| `>` | Redirect output (overwrite) | `ls > files.txt` |
| `>>` | Redirect output (append) | `echo "data" >> log.txt` |
| `\|` | Pipe to next command | `grep ">" file.fasta \| wc -l` |
| `<` | Redirect input | `sort < input.txt` |
| `2>` | Redirect errors | `find / -name "*.txt" 2> errors.log` |
| `2>&1` | Redirect errors to stdout | `command > all.log 2>&1` |
| `&>` | Redirect both (shorthand) | `command &> all.log` |
| `tee file` | Write to file AND stdout | `ls \| tee files.txt` |

---

## Wildcards

| Wildcard | Matches | Example |
|----------|---------|---------|
| `*` | Any characters | `*.fasta` (all FASTA files) |
| `?` | Single character | `sample?.txt` (sample1.txt, sample2.txt) |
| `[abc]` | One of a, b, or c | `file[123].txt` |
| `[a-z]` | Range (a through z) | `[A-Z]*.txt` |
| `{a,b,c}` | Exactly a, b, or c | `file.{txt,csv,tsv}` |

**Examples**:
```bash
ls *.fasta                       # All FASTA files
ls sample_*.txt                  # Files starting with "sample_"
ls file[1-9].txt                 # file1.txt through file9.txt
cat *.txt > combined.txt         # Combine all text files
mv *.fasta sequences/            # Move all FASTA files to directory
```

---

## Bioinformatics-Specific Commands

### FASTA File Operations

```bash
# Count sequences
grep -c "^>" sequences.fasta

# Extract sequence IDs
grep "^>" sequences.fasta | sed 's/>//'

# Extract sequence IDs (first word only)
grep "^>" sequences.fasta | sed 's/>//' | cut -d' ' -f1

# Get sequence without headers
grep -v "^>" sequences.fasta

# Count total nucleotides (excluding headers)
grep -v "^>" sequences.fasta | tr -d '\n' | wc -c

# Convert multi-line FASTA to single-line
awk '/^>/ {if(NR>1) print ""; print; next} {printf "%s", $0} END {print ""}' sequences.fasta

# Calculate average sequence length
awk '/^>/ {if(seqlen){sum+=seqlen; count++} seqlen=0; next}
     {seqlen+=length($0)}
     END {if(count>0) print sum/count}' sequences.fasta

# Find sequences containing specific motif
grep -v "^>" sequences.fasta | grep "ATGATG"
```

### CSV/TSV File Operations

```bash
# View first few rows with column alignment
head metadata.csv | column -t -s','

# Count rows (excluding header)
tail -n +2 metadata.csv | wc -l

# Extract specific columns
cut -d',' -f1,3,5 metadata.csv > subset.csv

# Sort by column
sort -t',' -k2 metadata.csv > sorted.csv

# Remove duplicate rows
sort metadata.csv | uniq > unique.csv

# Get unique values from column
cut -d',' -f2 metadata.csv | sort | uniq
```

### BLAST Results Processing

```bash
# Count unique hits
cut -f2 blast_results.txt | sort | uniq | wc -l

# Get top hit for each query
sort -k1,1 -k3,3rn blast_results.txt | awk '!seen[$1]++'

# Count hits per species
cut -f2 blast_results.txt | sort | uniq -c | sort -rn

# Filter by e-value (column 11 < 1e-50)
awk '$11 < 1e-50' blast_results.txt > filtered_blast.txt
```

---

## Common Pipelines

### Pipeline 1: Process FASTA
```bash
# Extract IDs, clean, sort, save
grep "^>" sequences.fasta | sed 's/>//' | cut -d' ' -f1 | sort > ids.txt
```

### Pipeline 2: Quality Control
```bash
# Count sequences with N's (ambiguous bases)
grep -v "^>" sequences.fasta | grep -c "N"
```

### Pipeline 3: Combine and Deduplicate
```bash
# Merge FASTA files and remove duplicates (by ID)
cat *.fasta | awk '/^>/ {id=$1; if(!seen[id]++) print; next} {print}'
```

### Pipeline 4: Species Summary
```bash
# From BLAST results, get species counts
cut -f2 blast_results.txt | sort | uniq -c | sort -rn | head -10
```

### Pipeline 5: Filter by Length
```bash
# Keep sequences longer than 500bp (requires awk)
awk 'BEGIN{RS=">"; ORS=""} NR>1 {split($0,a,"\n"); seq=""; for(i=2;i<=length(a);i++) seq=seq a[i]; if(length(seq)>500) print ">"$0}' sequences.fasta > filtered.fasta
```

---

## File Permissions

| Command | Description | Example |
|---------|-------------|---------|
| `chmod +x script.sh` | Make executable | `chmod +x run_analysis.sh` |
| `chmod 755 script.sh` | rwxr-xr-x | `chmod 755 script.sh` |
| `chmod 644 file.txt` | rw-r--r-- | `chmod 644 data.txt` |

**Permission Numbers**:
- 4 = read (r)
- 2 = write (w)
- 1 = execute (x)
- 7 = 4+2+1 = rwx
- 6 = 4+2 = rw-
- 5 = 4+1 = r-x

---

## Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| `Ctrl+C` | Stop current command |
| `Ctrl+D` | Exit terminal / EOF |
| `Ctrl+L` | Clear screen |
| `Ctrl+A` | Move to beginning of line |
| `Ctrl+E` | Move to end of line |
| `Ctrl+U` | Delete to beginning of line |
| `Ctrl+K` | Delete to end of line |
| `Ctrl+R` | Search command history |
| `Tab` | Autocomplete |
| `↑/↓` | Navigate command history |

---

## Useful One-Liners

```bash
# Find large files (>100MB)
find . -type f -size +100M

# Disk usage of directories
du -sh */

# Count files in directory
ls -1 | wc -l

# Find files modified in last 7 days
find . -type f -mtime -7

# Replace spaces with underscores in filenames
for f in *\ *; do mv "$f" "${f// /_}"; done

# Create backup with timestamp
cp important_file.txt important_file_$(date +%Y%m%d).txt

# Monitor file in real-time
tail -f logfile.txt

# Get file size
ls -lh file.txt | awk '{print $5}'

# Count unique values in column
cut -f2 data.txt | sort | uniq | wc -l
```

---

## Useful Variables

```bash
$HOME           # Your home directory
$PWD            # Current directory
$OLDPWD         # Previous directory
$USER           # Your username
$(date +%Y%m%d) # Today's date (YYYYMMDD)
$?              # Exit status of last command
```

---

## Getting Help

```bash
man command                      # Manual page for command
command --help                   # Quick help
apropos keyword                  # Search for commands
which command                    # Show command location
type command                     # Show command type
```

---

## Quick Tips

1. **Use Tab completion** - Start typing and press Tab
2. **Use `history`** - See all previous commands
3. **Use `!!`** - Repeat last command
4. **Use `!$`** - Last argument of previous command
5. **Use `Ctrl+R`** - Search command history
6. **Check before deleting** - Use `ls` before `rm`
7. **Test pipes step-by-step** - Build complex commands incrementally
8. **Read error messages** - They tell you what's wrong
9. **Google is your friend** - "bash how to [task]"
10. **Practice daily** - The only way to get comfortable

---

## Example Workflow

```bash
#!/bin/bash
# Complete analysis pipeline

# Navigate to project
cd ~/dna-barcoding-analysis

# Combine all samples
cat 04_data/student_sequences/*.fasta > combined.fasta

# Count sequences
total_seqs=$(grep -c "^>" combined.fasta)
echo "Total sequences: $total_seqs"

# Extract IDs
grep "^>" combined.fasta | sed 's/>//' | cut -d' ' -f1 > ids.txt

# Calculate average length
avg_length=$(awk '/^>/ {if(seqlen){sum+=seqlen; count++} seqlen=0; next}
                   {seqlen+=length($0)}
                   END {if(count>0) print sum/count}' combined.fasta)
echo "Average length: $avg_length bp"

# Find sequences with Ns
grep -v "^>" combined.fasta | grep -c "N" | tee sequences_with_Ns.txt

# Align sequences
mafft --auto combined.fasta > aligned.fasta

# Build tree
iqtree -s aligned.fasta -m MFP -bb 1000

echo "Analysis complete!"
```

---

**Print this cheatsheet and keep it handy!**

For detailed explanations, see the individual lesson files.
