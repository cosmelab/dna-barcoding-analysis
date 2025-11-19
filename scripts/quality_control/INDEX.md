# Quality Control Scripts - File Index

## Quick Navigation

### Getting Started (Read These First)
1. **QUICKSTART.txt** - 5 minute introduction to all scripts
2. **README.md** - Comprehensive guide with learning objectives

### Scripts (3 Python + 1 Bash)

#### Python Scripts
- **parse_ab1.py** - Parse ABI Sanger sequencer files (227 lines)
  - Extract sequences and quality scores from .ab1 binary files
  - Export to FASTA and QUAL formats
  - View metadata from sequencer

- **trim_quality.py** - Trim low-quality sequence ends (298 lines)
  - Remove poor quality bases from start and end
  - Work with FASTQ or FASTA+QUAL formats
  - Detailed trimming statistics

- **filter_sequences.py** - Filter sequences by quality (357 lines)
  - Remove entire bad sequences
  - Filter by: length, average quality, N characters
  - Object-oriented design with SequenceFilter class

#### Bash Script
- **batch_qc.sh** - Automate QC workflow (292 lines)
  - Process multiple .ab1 files
  - Complete pipeline: parse → trim → filter
  - Logging and reporting
  - Error handling

### Documentation (2,825 lines total)

- **README.md** (350 lines) - Full documentation with concepts
- **EXAMPLES.md** (400 lines) - Practical workflow scenarios
- **QUICKSTART.txt** (260 lines) - Quick reference guide
- **SUMMARY.md** (430 lines) - Overview of all components
- **INDEX.md** (this file) - File navigation

### Testing
- **test_scripts.sh** - Test suite to verify all scripts work

---

## Learning Path

### Absolute Beginner
Start here: **QUICKSTART.txt**
- Basic command syntax
- How to run each script
- When to use each one

### Getting Hands-On
Next: **EXAMPLES.md**
- Complete workflow examples
- Real project scenarios
- Parameter selection guide

### Understanding the Code
Read: **README.md**
- Why each script exists
- Bioinformatics concepts
- Quality score explanation

### Advanced Learning
Study: Script source code
- parse_ab1.py - File I/O
- trim_quality.py - Algorithms
- filter_sequences.py - OOP design
- batch_qc.sh - Bash scripting

### Reference
Use: **SUMMARY.md**
- Architecture overview
- Learning objectives
- Integration with other tools

---

## File Summary

### Executable Scripts (Runnable)
```
parse_ab1.py      - Extract from .ab1 files
trim_quality.py   - Trim by quality
filter_sequences.py - Filter by metrics
batch_qc.sh       - Batch processing
test_scripts.sh   - Run tests
```

### Documentation (Read)
```
README.md        - Main documentation
EXAMPLES.md      - Workflow examples
QUICKSTART.txt   - Quick reference
SUMMARY.md       - Overview
INDEX.md         - This file
```

---

## Common Tasks

### Task: I want to run QC on my sequences
1. Read **QUICKSTART.txt** (5 min)
2. Choose your script based on input format
3. Run the script with appropriate parameters
4. Check output statistics

### Task: I want to learn bioinformatics
1. Read **README.md** "Concepts Demonstrated" section
2. Study **parse_ab1.py** for file I/O concepts
3. Study **trim_quality.py** for algorithms
4. Study **filter_sequences.py** for OOP

### Task: I want to build a production pipeline
1. Read **EXAMPLES.md** "Batch Processing" section
2. Customize **batch_qc.sh** for your needs
3. Add error handling and logging
4. Test with real data

### Task: I want to understand quality scores
1. Read **README.md** "Understanding Quality Scores" section
2. Check **QUICKSTART.txt** "Quality Score Reference"
3. Run **trim_quality.py** with different thresholds
4. Compare results to understand impact

---

## Quick Reference

### Most Used Commands

**Parse AB1 file:**
```bash
python3 parse_ab1.py sample.ab1 --save-fasta output.fasta
```

**Trim sequences:**
```bash
python3 trim_quality.py input.fastq --quality 20 --output trimmed.fastq
```

**Filter sequences:**
```bash
python3 filter_sequences.py input.fastq --min-quality 20 --min-length 300
```

**Batch processing:**
```bash
bash batch_qc.sh /path/to/ab1_files
```

---

## Statistics

### Code
- **Total scripts:** 4 (3 Python + 1 Bash)
- **Total code lines:** 1,174 lines
- **Lines with comments:** ~40% (highly educational)
- **Functions/Classes:** 25+ well-documented

### Documentation
- **Total documentation:** 1,420 lines
- **Number of examples:** 20+ practical scenarios
- **Concepts covered:** 15+ bioinformatics topics

### Coverage
- DNA sequencing basics
- Quality scores (PHRED)
- File formats (FASTA, QUAL, FASTQ, AB1)
- Algorithm design
- Object-oriented programming
- Bash scripting
- Workflow automation
- Error handling
- Logging and reporting

---

## Dependencies

### Required
- Python 3.6+
- BioPython (`pip install biopython`)

### Optional
- GNU parallel (for fast parallel processing)
- BLAST (for downstream sequence analysis)
- Clustal/MAFFT (for alignment)
- RAxML/IQ-TREE (for phylogeny)

---

## Next Steps

After completing these scripts:

1. **Alignment** - Learn multiple sequence alignment
2. **Phylogeny** - Build evolutionary trees
3. **Identification** - Species assignment methods
4. **Visualization** - Create publication-quality figures
5. **Production** - Deploy real pipelines
6. **Best Practices** - QA/QC for bioinformatics

---

## Support & References

### Built-in Help
```bash
python3 parse_ab1.py -h
python3 trim_quality.py -h
python3 filter_sequences.py -h
bash batch_qc.sh --help  # Check script for help
```

### Key References in Documentation
- BioPython: https://biopython.org/
- FASTQ format: https://en.wikipedia.org/wiki/FASTQ_format
- PHRED quality: https://en.wikipedia.org/wiki/Phred_quality_score
- DNA barcoding: https://www.barcodinglife.org/

### Troubleshooting
See **README.md** section "Troubleshooting"
See **QUICKSTART.txt** "Troubleshooting Quick Reference"

---

## File Locations

All files are in: `/Users/lucianocosme/Projects/dna-barcoding-analysis/scripts/quality_control/`

### Scripts
- parse_ab1.py
- trim_quality.py
- filter_sequences.py
- batch_qc.sh
- test_scripts.sh

### Documentation
- README.md
- EXAMPLES.md
- QUICKSTART.txt
- SUMMARY.md
- INDEX.md (this file)

---

## How to Use This Index

1. **Find a task** → Look in "Common Tasks"
2. **Run a command** → Look in "Quick Reference"
3. **Learn a concept** → Look in corresponding documentation file
4. **Understand code** → Read the script source
5. **Troubleshoot** → Check README.md or QUICKSTART.txt

---

## Version Information

- Created: November 2024
- Python Version: 3.6+
- BioPython Version: 1.78+
- Bash Version: 4.0+
- Status: Educational release (fully functional)

---

## Tips for Students

1. Start with QUICKSTART.txt
2. Run test_scripts.sh to verify everything works
3. Try the examples in EXAMPLES.md
4. Read the source code and comments
5. Modify scripts to learn concepts
6. Document your modifications
7. Build your own tools using these as templates

Good luck with your bioinformatics journey!
