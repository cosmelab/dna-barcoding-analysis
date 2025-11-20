# Module 02: Consensus Sequences

## What are Consensus Sequences and Why Do We Need Them?

When you sequence DNA using Sanger sequencing, you typically sequence from **two directions**:

- **Forward (F) read**: Sequences from the 5' → 3' direction
- **Reverse (R) read**: Sequences from the 3' → 5' direction (the opposite strand)

**Why sequence both directions?**

1. **Accuracy**: If both reads agree on a base, you can be very confident it's correct
2. **Quality**: The middle region where F and R overlap typically has the highest quality
3. **Error correction**: If one read has a poor-quality region, the other might be good
4. **Professional standard**: Publications require bidirectional sequencing for DNA barcoding

A **consensus sequence** is created by combining the forward and reverse reads into a single, high-quality sequence that represents the true DNA sequence of your sample.

---

## Quick Start: Creating Consensus Sequences

### Basic Command

```bash
python create_consensus.py <input_fasta> <output_directory> [options]
```

### Example Usage

```bash
# Basic usage (keeps all sequences, even if only F or only R)
python create_consensus.py results/01_qc/passed_sequences.fasta results/02_consensus/

# Recommended: Only keep samples with BOTH F and R (highest quality)
python create_consensus.py results/01_qc/passed_sequences.fasta results/02_consensus/ --pairs-only

# Automatically open the HTML report when done
python create_consensus.py results/01_qc/passed_sequences.fasta results/02_consensus/ --pairs-only --open
```

### Important Flags

**`--pairs-only`** (RECOMMENDED)

Only create consensus sequences for samples that have BOTH forward AND reverse reads that passed QC.

- ✓ Higher quality consensus sequences
- ✓ More confident base calls
- ✓ Better for species identification and phylogenetics
- ✗ Fewer total sequences (only complete pairs)

**When to use `--pairs-only`:**
- For publication-quality data
- For species identification (DNA barcoding)
- For phylogenetic trees
- When you need maximum accuracy

**When NOT to use `--pairs-only`:**
- When you have very few samples and need every sequence
- For preliminary/exploratory analysis
- When you know one direction consistently fails for technical reasons

---

## How the Consensus Algorithm Works

### Step 1: Pair Matching

The script identifies forward/reverse pairs by sample name:

```
Input sequences:
- AT-HV1F  (forward read from sample AT-HV1)
- AT-HV1R  (reverse read from sample AT-HV1)
- AT-HV2F  (forward read from sample AT-HV2)
- AT-MA1R  (reverse read from sample AT-MA1, missing F)

Pairs identified:
✓ AT-HV1: F + R (complete pair)
✗ AT-HV2: F only (incomplete - missing R)
✗ AT-MA1: R only (incomplete - missing F)
```

**Naming requirements:**
- Forward reads must end with 'F' (e.g., Sample1F, AT-HV1F)
- Reverse reads must end with 'R' (e.g., Sample1R, AT-HV1R)
- Everything before the F/R must match exactly

### Step 2: Reverse Complement

The reverse read is converted to its reverse complement so both reads are in the same orientation (5' → 3').

```
Forward:  5'-ATGCATGCAT-3'
Reverse:  5'-ATGCATGCAT-3' (already reverse complemented)
```

### Step 3: Alignment and Consensus

The forward and reverse reads are aligned, and a consensus sequence is created:

**Rules for consensus calling:**
1. If both F and R agree on a base → use that base (high confidence)
2. If F and R disagree → use the base with higher quality score
3. If one position is a gap in F or R → use the non-gap base
4. If both are gaps → use gap

**Example:**
```
Position:  1 2 3 4 5 6 7 8 9 10
Forward:   A T G C A - G T A C
Reverse:   A T G C T C G T A C
Consensus: A T G C A C G T A C
           ^ ^ ^ ^ ^   ^ ^ ^ ^
           agree      disagree (use higher quality)
```

---

## Understanding the Outputs

### Output Files Created

After running consensus creation, you'll get:

#### 1. **consensus_report.html** (Interactive Report)

Open this in your web browser to see:

**Summary Table:**
- Shows each sample and its pairing status
- **Green rows**: Complete pairs (both F+R available)
- **Yellow rows**: Incomplete pairs (only F or only R)
- Columns show:
  - Sample name
  - Forward status (✓ or ✗)
  - Reverse status (✓ or ✗)
  - Consensus created (Yes/No)

**Pairing Statistics:**
```
Total samples: 15
Complete F+R pairs: 4 (26.7%)
Only Forward: 6 (40.0%)
Only Reverse: 5 (33.3%)

Consensus sequences created: 4
```

#### 2. **consensus_sequences.fasta** (Consensus FASTA File)

Contains the final consensus sequences ready for downstream analysis:

```
>AT-HV1_consensus
ATGCATGCATGCATGCAT...

>AT-HV3_consensus
GCTAGCTAGCTAGCTAGT...
```

**If using `--pairs-only`:** This file contains ONLY samples with complete F+R pairs.

**If NOT using `--pairs-only`:** This file contains all samples (some may be F-only or R-only, not true consensus).

#### 3. **pairing_report.json** (Machine-Readable Data)

JSON format data about the pairing results. Useful for:
- Custom analysis scripts
- Tracking which samples had complete pairs
- Quality control auditing

Example structure:
```json
{
  "sample": "AT-HV1",
  "has_forward": true,
  "has_reverse": true,
  "pair_status": "complete",
  "consensus_created": true
}
```

---

## Interpreting Your Results

### Good Results

```
Total samples: 15
Complete F+R pairs: 12 (80%)
Consensus sequences created: 12
```

**What this means:**
- Most of your samples had both F and R reads pass QC
- You'll have good quality consensus sequences
- Ready for downstream analysis

### Acceptable Results

```
Total samples: 15
Complete F+R pairs: 6 (40%)
Consensus sequences created: 6
```

**What this means:**
- About half your samples have complete pairs
- You'll have fewer sequences but good quality
- Some samples may need re-sequencing

### Poor Results

```
Total samples: 15
Complete F+R pairs: 2 (13%)
Consensus sequences created: 2
```

**What this means:**
- Very few complete pairs
- Likely a systematic problem with either F or R direction
- Check QC report to see which direction is failing
- May need to troubleshoot sequencing or primer issues

---

## Troubleshooting

### Problem: No F+R pairs found

**Possible causes:**
1. **Naming mismatch**: Forward and reverse reads don't have matching names
2. **One direction completely failed QC**: All F or all R reads failed quality control
3. **File naming issue**: Reads not ending in F/R

**What to do:**
- Check `results/01_qc/qc_report.html` - are both F and R passing?
- Verify your file naming: must end with F or R
- Example: `Sample1F.ab1` and `Sample1R.ab1` → names match

### Problem: Very few pairs (< 30%)

**Possible causes:**
1. **One primer failing**: Forward OR reverse primer not working well
2. **QC too strict**: Some good sequences failing QC thresholds
3. **Sample quality issues**: Some samples degraded before sequencing

**What to do:**
- Open QC report and check: are mostly F failing, or mostly R?
- If one direction consistently fails → primer or sequencing issue
- Consider re-sequencing samples with only one direction

### Problem: Many samples have only F or only R

**Possible causes:**
1. **Random QC failures**: Normal variation in sequence quality
2. **Sample-specific issues**: Some samples have poor DNA quality

**What to do:**
- This is somewhat normal (aim for >50% complete pairs)
- Prioritize re-sequencing samples where one direction was very close to passing
- Use `--pairs-only` flag to ensure only high-quality consensus sequences

### Problem: Consensus sequences seem wrong

**Possible causes:**
1. **Wrong reference orientation**: R read not properly reverse-complemented
2. **Contamination**: F and R from different samples
3. **Software bug**: Report this!

**What to do:**
- Check that F and R sequences are from the same sample
- Verify sample labeling at sequencing facility
- Compare consensus to expected sequence if you have a reference

---

## The `--pairs-only` Flag Explained

### Without `--pairs-only` (default):

```
Input:
- Sample1: F + R → Creates consensus ✓
- Sample2: F only → Uses F as "consensus" ⚠
- Sample3: R only → Uses R as "consensus" ⚠

Output: 3 sequences
Quality: Mixed (some true consensus, some single-direction)
```

### With `--pairs-only`:

```
Input:
- Sample1: F + R → Creates consensus ✓
- Sample2: F only → SKIPPED ✗
- Sample3: R only → SKIPPED ✗

Output: 1 sequence
Quality: High (all true bidirectional consensus)
```

**Recommendation:** Always use `--pairs-only` for final analysis unless you have very few samples.

---

## Tips for Better Consensus Sequences

### Before Sequencing

1. **Use both primers**: Always sequence with both forward and reverse primers
2. **Balance DNA amounts**: Send equal amounts of template for F and R reactions
3. **Same quality for both**: Don't prioritize one direction over the other

### During QC (Module 01)

1. **Check both directions**: Make sure F and R both pass QC
2. **Consistent quality**: Both F and R should have similar quality scores
3. **Re-sequence failures**: If only one direction fails, re-sequence just that direction

### During Consensus Creation

1. **Use `--pairs-only`**: Ensures highest quality data
2. **Check pairing stats**: Aim for >50% complete pairs
3. **Review the HTML report**: Look for systematic issues (all F failing, etc.)

---

## Next Steps: Moving to Module 03 (Combine with References)

Once you have consensus sequences, you're ready to combine them with reference sequences!

**Next step in the workflow:**
```bash
# Combine your consensus sequences with reference database
cat results/02_consensus/consensus_sequences.fasta \
    data/reference_sequences/socal_mosquitoes.fasta \
    > results/02_consensus/combined_with_references.fasta
```

This creates a single file with:
- YOUR consensus sequences (4 samples)
- Known reference sequences (52 SoCal mosquitoes)
- Total: 56 sequences ready for alignment and tree building

**Then proceed to Module 03: Alignment**

---

## Quick Reference

### Command Options

| Option | Description | When to Use |
|--------|-------------|-------------|
| `--pairs-only` | Only keep complete F+R pairs | Recommended for final analysis |
| `--open` | Auto-open HTML report | Convenient for immediate viewing |
| None | Keep all sequences | Exploratory analysis, low sample count |

### Quality Metrics

| Metric | Good | Acceptable | Poor |
|--------|------|-----------|------|
| Complete pairs % | > 70% | 40-70% | < 40% |
| Consensus created | All pairs | Most pairs | Few pairs |
| F+R agreement | > 95% | 90-95% | < 90% |

---

## Summary

- Consensus sequences combine F and R reads for maximum accuracy
- Use `--pairs-only` flag to ensure only high-quality consensus sequences
- Check the HTML report to see which samples have complete pairs
- Aim for >50% complete pairs (higher is better)
- Use `consensus_sequences.fasta` as input for alignment
- Always check pairing statistics before proceeding to next module

**Remember:** High-quality consensus sequences are the foundation for accurate species identification and phylogenetic analysis!
