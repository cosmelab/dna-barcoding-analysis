# Module 01: Quality Control (QC)

## What is Quality Control and Why Does It Matter?

Before we can use DNA sequences for identification or building evolutionary trees, we need to make sure they're good quality. **Quality Control (QC)** is the process of checking your DNA sequences to ensure they are accurate and reliable.

Think of it like this: A DNA sequencing machine reads millions of tiny DNA fragments and predicts which base (A, T, G, or C) is at each position. But sometimes the machine isn't 100% confident in its prediction. Quality Control helps us identify which parts of the sequence we can trust and which parts might have errors.

### Why is QC critical?

- **Garbage in, garbage out**: Bad quality data leads to wrong identifications
- **Wasting time**: If you analyze poor quality sequences, you might get misleading results
- **Saving samples**: QC helps you identify which samples need re-sequencing
- **Professional standard**: Every professional lab and research project uses QC as the first step

---

## Quick Start: Running Quality Control

### Basic Command

```bash
python qc_chromatograms.py <input_directory> [output_directory]
```

### Example Usage

```bash
# Analyze all .ab1 files in the 'sequences' folder and save results to 'results'
python qc_chromatograms.py data/sequences/ results/

# If you don't specify output directory, results go to a 'results' folder
python qc_chromatograms.py data/my_samples/
```

### What You Need

Your input directory should contain **`.ab1 chromatogram files`**. These are the raw output files from Sanger DNA sequencers (like ABI 3730). Each file represents one DNA sequence that was read from one sample.

```
data/sequences/
├── sample_001.ab1
├── sample_002.ab1
├── sample_003.ab1
└── sample_004.ab1
```

---

## Understanding the Outputs

### Output Files Created

After running QC, you'll get three output files:

#### 1. **qc_report.html** (Interactive Report)
This is the main report you'll look at! Open it in your web browser.

**What you see:**
- A table showing results for each sequence file
- Color-coded rows:
  - **Green background** = PASS (good quality sequence)
  - **Red background** = FAIL (poor quality sequence)
  - **Yellow background** = ERROR (couldn't read the file)
- Expandable chromatogram visualizations
- Your full DNA sequence for each file

**How to read it:**
1. Look for the green "PASS" rows - these sequences are ready to use
2. Check the columns:
   - **File**: Name of your .ab1 file
   - **Length (bp)**: How many DNA bases the sequence contains (longer is usually better)
   - **Avg Quality**: Average quality score (higher is better, see "Understanding Quality Metrics" below)
   - **High Quality %**: What percentage of bases passed quality thresholds (we want >80%)
   - **Status**: PASS or FAIL

**Viewing the chromatogram:**
Click "Show Chromatogram & Sequence" under any sequence to see:
- **Chromatogram plot**: The raw data from the sequencer. Four colored lines show the signal for each DNA base (A=green, C=blue, G=black, T=red)
- **Base calls**: The letters above the colored lines show which base was called at that position, color-coded by quality
- **Full sequence**: Your complete DNA sequence in text form

#### 2. **qc_results.json** (Machine-Readable Data)
This is the same data as the HTML report, but in JSON format (a standard format for computers to read). Use this if you want to:
- Write a custom script to process the results
- Load the data into another analysis program
- Create your own visualizations

Example structure:
```json
{
  "file": "sample_001.ab1",
  "length": 850,
  "avg_quality": 28.5,
  "high_quality_bases": 750,
  "low_quality_bases": 100,
  "percent_high_quality": 88.24,
  "qc_status": "PASS",
  "sequence": "ATGCATGCAT..."
}
```

#### 3. **passed_sequences.fasta** (FASTA File)
This file contains ONLY the sequences that passed QC, in FASTA format. This is what you'll use as input for the next module (Alignment).

```
>sample_001
ATGCATGCATGCATGC...
>sample_003
GCTAGCTAGCTAGCTA...
```

Notice that `sample_002` and `sample_004` are missing (they failed QC and won't be used further).

---

## Understanding Quality Metrics

### What are Phred Scores?

Every DNA base called by the sequencer gets a **Phred quality score** (or Q score). This is a number that represents the confidence of the machine in that base call.

- **Q20**: 99% confidence the base is correct (1 error per 100 bases)
- **Q30**: 99.9% confidence the base is correct (1 error per 1000 bases)
- **Q40**: 99.99% confidence the base is correct (1 error per 10,000 bases)

Higher Q scores = higher confidence = better quality.

### QC Thresholds: What Does "Passed" Mean?

A sequence **PASSES QC** if BOTH of these conditions are met:

1. **Average Quality ≥ 20**: The average Q score across all bases is at least 20
2. **High Quality Percentage ≥ 80%**: At least 80% of the bases have a Q score ≥ 20

### How to Interpret Your Results

**Example 1: PASS**
```
- Avg Quality: 28.5 (✓ Above 20)
- High Quality: 88% (✓ Above 80%)
- Status: PASS ✓
This sequence is good quality and ready to use!
```

**Example 2: FAIL**
```
- Avg Quality: 18.2 (✗ Below 20)
- High Quality: 72% (✗ Below 80%)
- Status: FAIL ✗
This sequence has too many low-quality bases.
```

### Why These Thresholds?

- Q20 is the industry standard for DNA barcoding (required by most journals)
- 80% high-quality bases ensures most of your sequence is reliable
- These are conservative thresholds - they're designed to catch bad sequences before they cause problems downstream

---

## What the Chromatogram Shows

When you click "Show Chromatogram & Sequence", you see:

```
Colored lines represent signal strength:
- Black line = Guanine (G)
- Green line = Adenine (A)
- Red line = Thymine (T)
- Blue line = Cytosine (C)
```

**Base call colors (the letters on top):**
- Dark green text = Q ≥ 30 (high confidence)
- Orange text = Q ≥ 20 (medium confidence)
- Red text = Q < 20 (low confidence)

**What a good chromatogram looks like:**
- Strong, clear peaks for each base
- Minimal noise/background signal
- Consistent signal throughout the sequence

**What a bad chromatogram looks like:**
- Weak or double peaks
- Lots of background noise
- Degraded signal at the beginning or end

---

## Troubleshooting: What to Do If Sequences Fail

### Problem: Most sequences are failing QC

**Possible causes:**
1. **Poor DNA quality at source**: The DNA extracted from your sample might be degraded or contaminated
2. **PCR issues**: Problems during DNA amplification can reduce quality
3. **Sequencing machine issues**: The sequencer might need calibration

**What to do:**
- Re-extract DNA from your samples (use fresh material if possible)
- Check your PCR conditions (annealing temperature, primer design)
- Consult with your sequencing facility

### Problem: First 50 bases or last 50 bases fail, but middle region is good

**This is normal!** Sanger sequencing typically produces lower quality at the very start and very end of reads. The QC script focuses on bases 50-200 for this reason.

**What to do:**
- This is expected - don't worry if your first/last ~50 bases are low quality
- Your sequence can still PASS because the middle region is good

### Problem: One specific sequence failed unexpectedly

**Possible causes:**
1. **Sample contamination**: Another DNA sequence mixed with yours
2. **Poor PCR amplification**: Not enough DNA was made
3. **Sample degradation**: DNA broke down before sequencing

**What to do:**
- Check the chromatogram - look for strange peaks or weak signal
- Try re-amplifying the sample with PCR
- If it's contaminated, you may need to manually sequence the region of interest

### Problem: Getting an ERROR instead of PASS/FAIL

**This means the script couldn't read the .ab1 file.** Possible causes:
1. File is corrupted
2. File is not actually an .ab1 format
3. Insufficient permissions to read the file

**What to do:**
- Check that the file is a real .ab1 file from your sequencer
- Try re-exporting the file from your sequencing software
- Contact your sequencing facility if the file is corrupted

---

## Tips for Better Results

### Before You Sequence

1. **Start with good DNA quality**: Extract DNA carefully and check it with a spectrophotometer
2. **Use good primers**: Ensure primers are specific to your target gene
3. **Optimize PCR**: Test your PCR reaction before sending samples

### During Sequencing

1. **Send enough DNA**: Follow your sequencing facility's recommendations
2. **Clean up your PCR products**: Use a purification kit before sequencing
3. **Request bidirectional reads**: Sequence from both forward and reverse primers to get more high-quality bases in the middle

### After Sequencing

1. **Always run QC**: Every sequence, every time - no exceptions!
2. **Keep passing sequences**: Save the `passed_sequences.fasta` file for the next step
3. **Document failures**: Note which samples failed and why - helps you troubleshoot

---

## Next Steps: Moving to Module 02 (Alignment)

Once you have your high-quality sequences in `passed_sequences.fasta`, you're ready for the next module!

**Module 02: Alignment** will:
1. Take your good-quality sequences
2. Align them to a reference sequence
3. Find the optimal reading frame
4. Prepare data for identification

Go to `modules/03_alignment/` to get started.

**Remember**: Good QC now saves time and prevents errors later. Never skip this step!

---

## Quick Reference: Key Metrics

| Metric | Good | Acceptable | Poor |
|--------|------|-----------|------|
| Average Quality | ≥ 30 | 20-29 | < 20 |
| High Quality % | ≥ 90% | 80-89% | < 80% |
| Sequence Length | > 600 bp | 400-600 bp | < 400 bp |
| QC Status | PASS | - | FAIL |

---

## Getting Help

If you encounter issues:

1. **Check the HTML report** - Look at the chromatogram visualization for clues
2. **Review the error message** - The error description often tells you what went wrong
3. **Check your input files** - Ensure they are real .ab1 files
4. **Ask your instructor or lab manager** - They can help troubleshoot

---

## Summary

- Quality Control filters out bad sequences before they cause problems
- Run the QC script on all your .ab1 chromatogram files
- Review the HTML report to see which sequences passed
- Use the `passed_sequences.fasta` file as input for Module 02
- Always check the chromatogram visualization for context
- Better QC now = better results later!
