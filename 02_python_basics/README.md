# Module 02: Python Basics for Bioinformatics

**Duration**: 3-4 hours
**Difficulty**: Beginner
**Prerequisites**: Module 01 (Linux Basics)

---

## Overview

Learn Python programming specifically for bioinformatics applications. This module focuses on practical skills: parsing biological file formats, manipulating sequences, and automating analysis tasks.

**Why Python for Bioinformatics?**
- Industry standard for computational biology
- BioPython library handles most file formats
- Easy to learn, powerful for automation
- Excellent for data analysis and visualization
- Works seamlessly with command-line tools

**What You'll Learn**:
- Python syntax and data types
- Read/write FASTA, GenBank, and CSV files
- Sequence manipulation with BioPython
- Quality score analysis
- Data visualization with matplotlib
- Write reusable bioinformatics scripts

---

## Learning Objectives

By the end of this module, you will be able to:

1. Write basic Python scripts from scratch
2. Parse FASTA files and extract information
3. Calculate sequence statistics (GC content, length, etc.)
4. Process quality scores from sequencing data
5. Visualize sequence data with plots
6. Automate repetitive bioinformatics tasks
7. Integrate Python scripts into analysis pipelines

---

## Module Contents

```
02_python_basics/
├── README.md                           # This file
├── notebooks/                          # Jupyter notebooks (interactive)
│   ├── 01_python_syntax.ipynb         # Variables, loops, functions
│   ├── 02_biopython_basics.ipynb      # BioPython introduction
│   ├── 03_fasta_parsing.ipynb         # Read/write FASTA files
│   ├── 04_sequence_manipulation.ipynb # Reverse complement, translation
│   ├── 05_quality_scores.ipynb        # Phred scores, filtering
│   └── 06_visualization.ipynb         # Plots and figures
├── exercises/                          # Practice problems
│   ├── exercise_01_basics.md
│   ├── exercise_02_fasta.md
│   ├── exercise_03_sequences.md
│   ├── exercise_04_quality.md
│   └── exercise_05_final_project.md
├── solutions/                          # Answer scripts
│   └── [solution files]
├── data/                               # Practice datasets
│   ├── sequences.fasta
│   ├── quality_scores.txt
│   └── metadata.csv
└── scripts/                            # Example scripts
    ├── count_sequences.py
    ├── calculate_gc_content.py
    ├── filter_by_length.py
    └── batch_process.py
```

---

## Quick Start

### Using Jupyter Lab (Recommended)

```bash
# Start container with Jupyter
docker run -p 8888:8888 -v $(pwd):/workspace dna-barcoding jupyter-lab

# Open browser to: http://localhost:8888
# Navigate to 02_python_basics/notebooks/
# Start with 01_python_syntax.ipynb
```

### Using Python Scripts

```bash
# Navigate to module
cd 02_python_basics

# Run example script
python scripts/count_sequences.py data/sequences.fasta

# Create your own script
nano my_script.py
# Write code, save, run:
python my_script.py
```

---

## Topics Covered

### 1. Python Syntax (Notebook 01)
- Variables and data types
- Lists, dictionaries, sets
- Loops (for, while)
- Conditionals (if/else)
- Functions
- File I/O

### 2. BioPython Basics (Notebook 02)
- Installing BioPython
- SeqIO module
- Seq objects
- SeqRecord objects
- Working with sequences

### 3. FASTA Parsing (Notebook 03)
- Read FASTA files
- Write FASTA files
- Extract sequence IDs
- Filter sequences
- Combine multiple FASTA files

### 4. Sequence Manipulation (Notebook 04)
- Calculate GC content
- Reverse complement
- Translation (DNA to protein)
- Find motifs
- Sequence slicing

### 5. Quality Scores (Notebook 05)
- Understanding Phred scores
- Parse .ab1 files
- Filter by quality
- Trim low-quality ends
- Quality visualization

### 6. Visualization (Notebook 06)
- matplotlib basics
- Plot sequence lengths
- GC content distribution
- Quality score plots
- Publication-ready figures

---

## Essential Python for Bioinformatics

### Example 1: Count Sequences in FASTA

```python
from Bio import SeqIO

def count_sequences(fasta_file):
    """Count number of sequences in FASTA file."""
    count = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        count += 1
    return count

# Usage
num_seqs = count_sequences("sequences.fasta")
print(f"Total sequences: {num_seqs}")
```

### Example 2: Calculate GC Content

```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def calculate_gc(fasta_file):
    """Calculate GC content for each sequence."""
    for record in SeqIO.parse(fasta_file, "fasta"):
        gc = gc_fraction(record.seq) * 100
        print(f"{record.id}: {gc:.2f}% GC")

calculate_gc("sequences.fasta")
```

### Example 3: Filter by Length

```python
from Bio import SeqIO

def filter_by_length(input_file, output_file, min_length=500):
    """Keep only sequences longer than min_length."""
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            if len(record.seq) >= min_length:
                SeqIO.write(record, output_handle, "fasta")

filter_by_length("sequences.fasta", "filtered.fasta", min_length=500)
```

### Example 4: Reverse Complement

```python
from Bio import SeqIO

def reverse_complement_all(input_file, output_file):
    """Generate reverse complement of all sequences."""
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            record.seq = record.seq.reverse_complement()
            record.id = record.id + "_RC"
            record.description = "reverse complement"
            SeqIO.write(record, output_handle, "fasta")

reverse_complement_all("sequences.fasta", "sequences_RC.fasta")
```

### Example 5: Quality Control Report

```python
from Bio import SeqIO
import matplotlib.pyplot as plt

def qc_report(fasta_file):
    """Generate QC report with statistics and plots."""
    lengths = []
    gc_contents = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths.append(len(record.seq))
        gc_contents.append(gc_fraction(record.seq) * 100)

    # Statistics
    print(f"Number of sequences: {len(lengths)}")
    print(f"Average length: {sum(lengths)/len(lengths):.1f} bp")
    print(f"Average GC content: {sum(gc_contents)/len(gc_contents):.2f}%")

    # Plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    ax1.hist(lengths, bins=30, edgecolor='black')
    ax1.set_xlabel('Sequence Length (bp)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Sequence Length Distribution')

    ax2.hist(gc_contents, bins=30, edgecolor='black')
    ax2.set_xlabel('GC Content (%)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('GC Content Distribution')

    plt.tight_layout()
    plt.savefig('qc_report.png', dpi=300)
    plt.show()

qc_report("sequences.fasta")
```

---

## Common BioPython Patterns

### Reading FASTA Files

```python
from Bio import SeqIO

# Method 1: Iterate through records
for record in SeqIO.parse("sequences.fasta", "fasta"):
    print(record.id, len(record.seq))

# Method 2: Load all into memory
records = list(SeqIO.parse("sequences.fasta", "fasta"))
print(f"Loaded {len(records)} sequences")

# Method 3: Index large files (doesn't load all into memory)
record_dict = SeqIO.index("sequences.fasta", "fasta")
print(record_dict["seq_001"].seq)
```

### Writing FASTA Files

```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Create a sequence record
record = SeqRecord(
    Seq("ATGATGATGATG"),
    id="my_sequence",
    description="Example sequence"
)

# Write to file
with open("output.fasta", "w") as output_handle:
    SeqIO.write(record, output_handle, "fasta")

# Write multiple records
records = [record1, record2, record3]
SeqIO.write(records, "output.fasta", "fasta")
```

---

## Tips for Success

1. **Use Jupyter Notebooks for Learning**
   - Interactive, see results immediately
   - Can mix code, text, and plots
   - Great for exploration

2. **Write Reusable Functions**
   - Break code into small functions
   - Test each function separately
   - Easier to debug and maintain

3. **Comment Your Code**
   - Future you will thank you
   - Explain *why*, not just *what*
   - Use docstrings for functions

4. **Test with Small Files First**
   - Debug with 10 sequences
   - Then run on full dataset
   - Saves time and frustration

5. **Use Virtual Environments**
   - Keep project dependencies isolated
   - Container already has everything
   - For local: `python -m venv venv`

---

## Common Errors and Solutions

### ImportError: No module named 'Bio'
```bash
# Solution: Install BioPython
pip install biopython

# Or use the container (already installed)
docker run -it dna-barcoding python
```

### File Not Found
```python
# WRONG - relative path might not work
with open("data.fasta") as f:

# BETTER - check file exists first
import os
if os.path.exists("data.fasta"):
    with open("data.fasta") as f:

# BEST - use absolute path or verify pwd
import os
print(os.getcwd())  # Where am I?
```

### Memory Error with Large Files
```python
# WRONG - loads entire file into memory
records = list(SeqIO.parse("huge_file.fasta", "fasta"))

# CORRECT - iterate without loading all
for record in SeqIO.parse("huge_file.fasta", "fasta"):
    # Process one record at a time
    pass
```

---

## Practice Exercises

See the `exercises/` directory for detailed problems:

1. **Exercise 01**: Python syntax basics
2. **Exercise 02**: FASTA file manipulation
3. **Exercise 03**: Sequence statistics
4. **Exercise 04**: Quality control
5. **Exercise 05**: Final project (complete analysis script)

---

## Additional Resources

- [BioPython Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
- [Python for Biologists](https://pythonforbiologists.com/)
- [Rosalind Bioinformatics Problems](http://rosalind.info/)
- [Real Python Tutorials](https://realpython.com/)

---

## Next Steps

After mastering Python basics, proceed to:

- **Module 03**: R Basics for Phylogenetics
- **Module 05**: Quality Control (apply Python skills)
- **Module 08**: Species Identification (Python scripts)

---

**Ready to code? Start with the Jupyter notebooks!** →
