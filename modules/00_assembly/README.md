# Module 00: Forward/Reverse Sequence Assembly

## What Does This Module Do?

When you sequence DNA using Sanger sequencing, you typically get TWO reads from each sample: one read from the **forward** direction and one from the **reverse** direction. This module takes those two reads and combines them into a single, longer DNA sequence (consensus). Think of it like reading a sentence from left to right and also right to left, then combining both reads to get the complete picture!

---

## Why Is This Important?

- **Longer sequences** = More accurate identification
- **Better coverage** = More reliable results
- **Quality check** = If F and R disagree, something might be wrong

Without this assembly step, your DNA sequences would be too short for accurate species identification.

---

## What You Need (Inputs)

### File Format: AB1 Chromatograms

You need `.ab1` files from your sequencer. These contain:
- The DNA sequence (bases: A, T, G, C)
- Quality scores for each base
- Trace data (the squiggly lines you see)

### Important: Naming Convention

Your files **MUST** be named with `_F` or `_R` suffixes:

```
Sample01_F.ab1    ← Forward read
Sample01_R.ab1    ← Reverse read

Sample02_F.ab1    ← Forward read
Sample02_R.ab1    ← Reverse read
```

**Good examples:**
```
mosquito_F.ab1
mosquito_R.ab1

AedesAegypti_F.ab1
AedesAegypti_R.ab1
```

**Bad examples (won't work):**
```
mosquito_forward.ab1     ← Wrong suffix
mosquito.ab1             ← No direction marker
mosquitoF.ab1            ← Missing underscore
```

---

## How to Use This Module

### Step 1: Prepare Your Files

Place your AB1 files in a folder. For example:
```
data/
└── my_samples/
    ├── Sample01_F.ab1
    ├── Sample01_R.ab1
    ├── Sample02_F.ab1
    └── Sample02_R.ab1
```

### Step 2: Run the Script

From your project directory, use the command line:

```bash
python modules/00_assembly/merge_forward_reverse.py data/my_samples/ results/
```

**What does this command do?**
- `python` = Run Python code
- `modules/00_assembly/merge_forward_reverse.py` = The script file
- `data/my_samples/` = WHERE your AB1 files are
- `results/` = WHERE to put the output files

---

## What You'll Get (Outputs)

The script creates **3 files** in your results folder:

### 1. `consensus_sequences.fasta`
The most important file! Contains your assembled sequences in FASTA format:
```
>Sample01
ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA...

>Sample02
ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

This FASTA file is what you'll use for the next module (Quality Control).

### 2. `assembly_stats.json`
Technical information in machine-readable format:
- Number of merged pairs
- Number of forward-only sequences
- Number of reverse-only sequences
- Length of each consensus sequence

### 3. `assembly_report.html`
A nice summary report you can open in your web browser. Shows:
- How many sequences were successfully assembled
- Status of each sample (merged, forward-only, or reverse-only)
- Length of each assembled sequence

---

## How Does It Work? (The Simple Version)

### What Is a "Reverse Complement"?

DNA has direction. When we read it backwards, we need to flip the bases:
- A ↔ T (complementary)
- G ↔ C (complementary)

```
Forward:  5'→ ATGC →3'
Reverse:  3'→ GCAT →5' (as read by the machine)

When reversed:
5'→ TACG →3' (this is called the "reverse complement")
```

### The Assembly Process

1. **Read both sequences** (Forward and Reverse)
2. **Reverse-complement the reverse sequence** (flip it around)
3. **Combine them** into one long sequence
4. **Report what happened** (merged? missing one? both?)

### Three Possible Outcomes

| Outcome | What Happened | What We Do |
|---------|---------------|-----------|
| **Merged (F+R)** | You have both forward and reverse | Combine them - BEST! |
| **Forward only** | You're missing the reverse file | Use just the forward read |
| **Reverse only** | You're missing the forward file | Use just the reverse (flipped) |

---

## Troubleshooting: Common Problems

### Problem 1: "Input file not found"

**Error message:**
```
Input file not found: data/my_samples/
```

**Solution:**
- Check the path is correct
- Make sure your AB1 files are actually in that folder
- Use `ls data/my_samples/` to see what's there

---

### Problem 2: "Unpaired sequences (skipped)"

**Error message:**
```
Unpaired sequences (skipped): Sample01_F
```

**This means:** You have a Sample01_F file, but NO Sample01_R file.

**Solutions:**
- Check your file names match the naming convention (underscore before F/R)
- Make sure you have BOTH forward and reverse for each sample
- Example: If you have `Sample01_F.ab1`, you need `Sample01_R.ab1`

---

### Problem 3: File naming is wrong

**Common mistakes:**
```
❌ Sample01_forward.ab1     (says "forward" instead of "F")
❌ Sample01f.ab1            (lowercase, no underscore)
❌ Sample01-F.ab1           (hyphen instead of underscore)
❌ Sample01F.ab1            (no underscore at all)
```

**Correct format:**
```
✓ Sample01_F.ab1
✓ Sample01_R.ab1
```

---

### Problem 4: Empty output files

**This usually means:** None of your sequences matched the naming convention.

**Solutions:**
1. Double-check your file names
2. Make sure you used underscore: `_F` and `_R`
3. Check uppercase/lowercase (usually doesn't matter, but F and R must be capital)

---

## Next Step: Module 01 (Quality Control)

After assembly, your `consensus_sequences.fasta` file goes to **Module 01** where we:
- Check the quality of each base
- Remove low-quality regions
- Trim sequences to keep only the good parts

**Command for next step:**
```bash
python modules/01_quality_control/qc_chromatograms.py results/consensus_sequences.fasta results/
```

---

## Quick Reference

| Step | File | Purpose |
|------|------|---------|
| **Input** | `*.ab1` files | Raw chromatograms with `_F` and `_R` suffixes |
| **This Module** | `merge_forward_reverse.py` | Assemble F+R into consensus |
| **Output** | `consensus_sequences.fasta` | Your merged DNA sequences |
| **Next** | Module 01 | Quality control of assembled sequences |

---

## Questions?

If something doesn't work:
1. Check the error message carefully
2. Verify your file names follow the `SampleName_F.ab1` format
3. Make sure BOTH forward and reverse files exist for each sample
4. Check the `assembly_report.html` file for a summary

Still stuck? Ask your instructor!
