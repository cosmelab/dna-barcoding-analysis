# Module 01: Linux Command Line Basics

**Duration**: 2-3 hours
**Difficulty**: Beginner (no prior experience required)
**Prerequisites**: None

---

## Overview

This module teaches essential Linux command-line skills for bioinformatics. You will learn to navigate filesystems, manipulate files, process text data, and chain commands together using pipes.

**Why Learn Command Line?**
- Most bioinformatics tools run on Linux
- Much faster than clicking through GUIs for repetitive tasks
- Required for working on high-performance computing (HPC) clusters
- Industry standard for computational biology
- Reproducible workflows through scripts

**What You'll Learn**:
- Navigate Linux filesystem confidently
- Create, move, copy, and delete files
- View and search file contents
- Process biological data files (FASTA, CSV, etc.)
- Use pipes to chain commands together
- Redirect input/output to files
- Use wildcards for pattern matching

---

## Learning Objectives

By the end of this module, you will be able to:

1. Navigate directories using `cd`, `ls`, and `pwd`
2. Create and organize files and directories
3. View file contents using `cat`, `less`, `head`, and `tail`
4. Search files with `grep`
5. Process text with `sed`, `awk`, `cut`, and `paste`
6. Redirect output with `>`, `>>`, and `|`
7. Use wildcards (`*`, `?`) for batch operations
8. Write simple bash scripts to automate tasks

---

## Module Contents

```
01_linux_basics/
├── README.md                    # This file
├── 01_navigation.md             # cd, ls, pwd, tree
├── 02_file_manipulation.md      # cp, mv, rm, mkdir, touch
├── 03_viewing_files.md          # cat, less, head, tail, wc
├── 04_searching.md              # grep, find
├── 05_text_processing.md        # sed, awk, cut, paste, sort, uniq
├── 06_streams_redirection.md    # >, >>, |, <, tee
├── 07_wildcards_patterns.md     # *, ?, [], {}
├── 08_scripting_basics.md       # bash scripts, variables, loops
├── cheatsheet.md                # Quick reference guide
├── exercises/                   # Practice problems
│   ├── exercise_01_navigation.md
│   ├── exercise_02_file_manipulation.md
│   ├── exercise_03_viewing.md
│   ├── exercise_04_searching.md
│   ├── exercise_05_text_processing.md
│   ├── exercise_06_pipes.md
│   ├── exercise_07_wildcards.md
│   └── exercise_08_final_project.md
├── solutions/                   # Answer keys
│   └── [corresponding solution files]
└── data/                        # Practice datasets
    ├── sequences.fasta
    ├── sample_metadata.csv
    ├── blast_results.txt
    └── README.md
```

---

## Quick Start

### Option 1: Using Docker Container (Recommended)

```bash
# Start the DNA Barcoding container
cd /path/to/dna-barcoding-analysis
docker run -it -v $(pwd):/workspace -w /workspace dna-barcoding

# Navigate to this module
cd 01_linux_basics

# You're ready to go!
```

### Option 2: Using Your Own Linux/Mac Terminal

```bash
# Open Terminal application
# Navigate to this module
cd /path/to/dna-barcoding-analysis/01_linux_basics

# Start with lesson 1
cat 01_navigation.md
```

### Option 3: Windows (WSL)

```bash
# Install Windows Subsystem for Linux (WSL2)
# Open Ubuntu terminal
# Clone this repository
git clone https://github.com/cosmelab/dna-barcoding-analysis.git
cd dna-barcoding-analysis/01_linux_basics
```

---

## Learning Path

### For Complete Beginners (3-4 hours)

1. **Read**: `01_navigation.md` (30 min)
   - Practice navigating directories
   - Complete Exercise 1

2. **Read**: `02_file_manipulation.md` (30 min)
   - Create, move, copy files
   - Complete Exercise 2

3. **Read**: `03_viewing_files.md` (30 min)
   - Learn to view file contents
   - Complete Exercise 3

4. **Read**: `04_searching.md` (30 min)
   - Search for text in files
   - Complete Exercise 4

5. **Read**: `05_text_processing.md` (45 min)
   - Process FASTA and CSV files
   - Complete Exercise 5

6. **Read**: `06_streams_redirection.md` (30 min)
   - Master pipes and redirection
   - Complete Exercise 6

7. **Read**: `07_wildcards_patterns.md` (20 min)
   - Use wildcards for batch operations
   - Complete Exercise 7

8. **Read**: `08_scripting_basics.md` (30 min)
   - Write your first bash script
   - Complete Exercise 8 (Final Project)

### For Those with Some Experience (1-2 hours)

1. Skim lessons 01-03
2. Focus on lessons 05-06 (text processing and pipes)
3. Complete exercises 5-8
4. Use `cheatsheet.md` as reference

---

## Essential Commands Summary

### Navigation
- `pwd` - Print working directory
- `ls` - List directory contents
- `cd` - Change directory
- `tree` - Display directory structure

### File Manipulation
- `mkdir` - Make directory
- `touch` - Create empty file
- `cp` - Copy files/directories
- `mv` - Move or rename files
- `rm` - Remove files/directories

### Viewing Files
- `cat` - Concatenate and display files
- `less` - Page through file
- `head` - View first lines
- `tail` - View last lines
- `wc` - Word, line, character count

### Searching
- `grep` - Search for patterns
- `find` - Find files by name/properties

### Text Processing
- `sed` - Stream editor
- `awk` - Pattern scanning and processing
- `cut` - Extract columns
- `paste` - Merge files horizontally
- `sort` - Sort lines
- `uniq` - Remove duplicates
- `tr` - Translate characters

### Redirection
- `>` - Redirect output (overwrite)
- `>>` - Redirect output (append)
- `|` - Pipe output to next command
- `<` - Redirect input from file
- `tee` - Write to file AND stdout

---

## Bioinformatics Use Cases

### Example 1: Count sequences in FASTA file
```bash
grep -c "^>" sequences.fasta
```

### Example 2: Extract sequence IDs
```bash
grep "^>" sequences.fasta | sed 's/>//' > sequence_ids.txt
```

### Example 3: Get sequence lengths
```bash
awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen+=length($0)} END {print seqlen}' sequences.fasta
```

### Example 4: Filter sequences by length
```bash
# Using awk to keep only sequences > 500bp
awk '/^>/ {if (seqlen > 500) {print header; print seq} header=$0; seq=""; seqlen=0; next} {seq=seq$0; seqlen+=length($0)} END {if (seqlen > 500) {print header; print seq}}' sequences.fasta
```

### Example 5: Combine multiple FASTA files
```bash
cat sample1.fasta sample2.fasta sample3.fasta > all_samples.fasta
```

### Example 6: Process CSV metadata
```bash
# Extract specific columns
cut -d',' -f1,3,5 metadata.csv > subset.csv

# Sort by column
sort -t',' -k2 metadata.csv > sorted_metadata.csv
```

---

## Tips for Success

### 1. Use Tab Completion
- Start typing a filename and press Tab
- The shell will autocomplete or show options
- Saves time and prevents typos

### 2. Use Command History
- Press Up Arrow to see previous commands
- Press Ctrl+R to search command history
- Type `history` to see all recent commands

### 3. Read Error Messages
- Error messages tell you what went wrong
- Look for the actual error (often at the bottom)
- Google the error message if stuck

### 4. Practice Regularly
- The only way to learn is by doing
- Try commands multiple times
- Experiment with different options

### 5. Always Check Before Deleting
- Use `ls` to verify files before `rm`
- Use `rm -i` for interactive confirmation
- There is no "Recycle Bin" in Linux!

### 6. Work in a Safe Directory
- Don't practice in your home directory
- Create a `practice/` folder
- Delete and recreate it when needed

### 7. Use `man` Pages
- Type `man command` to see manual
- Example: `man grep`
- Press `q` to quit

---

## Common Mistakes to Avoid

### 1. No Spaces in Variable Assignment
```bash
# WRONG
name = "John"

# RIGHT
name="John"
```

### 2. Use Quotes for Filenames with Spaces
```bash
# WRONG
cat my file.txt

# RIGHT
cat "my file.txt"
# Or avoid spaces altogether
cat my_file.txt
```

### 3. Forgetting File Extensions
```bash
# FASTA files might be .fasta, .fa, .fna
# Always check with ls first
ls *.fa*
```

### 4. Overwriting Files with Redirection
```bash
# DANGER: This will erase data.txt
command > data.txt

# SAFER: Append instead
command >> data.txt

# Or use different filename
command > data_processed.txt
```

### 5. Case Sensitivity
```bash
# Linux is case-sensitive
File.txt ≠ file.txt ≠ FILE.TXT
```

---

## Keyboard Shortcuts

- `Ctrl+C` - Stop current command
- `Ctrl+D` - Exit terminal (or end of input)
- `Ctrl+L` - Clear screen (or type `clear`)
- `Ctrl+A` - Move cursor to beginning of line
- `Ctrl+E` - Move cursor to end of line
- `Ctrl+U` - Delete from cursor to beginning
- `Ctrl+K` - Delete from cursor to end
- `Ctrl+R` - Search command history
- `Tab` - Autocomplete
- `Up/Down Arrows` - Navigate command history

---

## Assessment

After completing this module, you should be able to:

- [ ] Navigate to any directory on your system
- [ ] List files with various options (`ls -lh`, `ls -a`)
- [ ] Create directories and files
- [ ] Copy, move, and delete files safely
- [ ] View file contents using appropriate tools
- [ ] Search for text patterns in files
- [ ] Extract specific columns from CSV files
- [ ] Count sequences in FASTA files
- [ ] Use pipes to chain commands together
- [ ] Redirect output to files
- [ ] Use wildcards to process multiple files
- [ ] Write a simple bash script

**Test Yourself**: Complete Exercise 8 (Final Project) without looking at solutions!

---

## Next Steps

After mastering Linux basics, proceed to:

- **Module 02**: Python Basics for Bioinformatics
- **Module 03**: R Basics for Phylogenetics
- **Module 04**: Working with Reference Data
- **Module 05**: Quality Control of Sequences

---

## Additional Resources

### Online Tutorials
- [Explain Shell](https://explainshell.com/) - Paste command, get explanation
- [Linux Journey](https://linuxjourney.com/) - Interactive tutorials
- [Command Line Challenge](https://cmdchallenge.com/) - Practice problems

### Cheat Sheets
- See `cheatsheet.md` in this directory
- [Linux Command Line Cheat Sheet (PDF)](https://www.cheatography.com/davechild/cheat-sheets/linux-command-line/)

### Books
- "The Linux Command Line" by William Shotts (free PDF)
- "Unix and Perl Primer for Biologists" by Bradnam & Korf

### Videos
- Linux Tutorial for Beginners (YouTube)
- Command Line Basics for Bioinformatics

---

## Troubleshooting

### Problem: Permission Denied
```bash
# Solution: Use sudo (if you have admin rights)
sudo command

# Or change file permissions
chmod +x script.sh
```

### Problem: Command Not Found
```bash
# Solution: Check if installed
which command_name

# Install using package manager (if in container, already installed)
# Ubuntu/Debian
apt-get install package-name

# macOS
brew install package-name
```

### Problem: File Not Found
```bash
# Solution: Check if you're in the right directory
pwd
ls

# Use absolute path instead of relative
/full/path/to/file.txt
```

---

## Practice Datasets

The `data/` directory contains:

1. **sequences.fasta** - Sample COI barcode sequences (50 sequences)
2. **sample_metadata.csv** - Specimen information
3. **blast_results.txt** - Example BLAST output
4. **quality_scores.txt** - Phred scores from sequencing
5. **alignment.fasta** - Aligned sequences

Use these files to practice all the commands in the exercises.

---

## Getting Help

- Re-read the lesson files
- Check the `cheatsheet.md`
- Look at the solutions (but try first!)
- Use `man command` for detailed help
- Ask on the course discussion board
- Attend office hours

---

## Final Notes

**Remember**: Everyone struggles with command line at first. It takes practice! Don't get discouraged if you make mistakes - that's how you learn.

**Goal**: By the end of this module, you should feel comfortable opening a terminal and not being intimidated by the command line. You won't know everything, but you'll know how to find answers.

**Philosophy**: Learn by doing. Type every command yourself. Don't copy-paste until you understand what it does.

---

**Happy Learning!**

Start with `01_navigation.md` →
