# Module 02: Multiple Sequence Alignment

## What is Sequence Alignment?

Sequence alignment is the process of arranging DNA sequences so that similar regions line up with each other. Think of it like comparing multiple texts side-by-side to find matching words and phrases.

### Before and After Alignment

**Before alignment** (unaligned sequences):
```
Sequence 1: ATGCGATCG
Sequence 2: ATGCGTAATCG
Sequence 3: ATGATCG
```

**After alignment** (aligned sequences):
```
Sequence 1: ATG-CGATC--G
Sequence 2: ATGCGTAATCG
Sequence 3: ATG--ATCG---
```

Notice how the sequences are now lined up so that similar positions (like `ATG` at the start) match across all sequences. Gaps (the `-` symbols) are added where sequences differ to keep everything aligned.

## Why Do We Need Alignment?

Alignments are crucial for several reasons:

- **Finding Similarities**: Identifies which parts of sequences are conserved (similar) across different organisms
- **Detecting Variations**: Shows where sequences differ, which helps identify species or genetic differences
- **Evolutionary Studies**: Aligned sequences can be used to build evolutionary trees (phylogeny)
- **Gene Identification**: Helps find specific genes or regulatory regions by comparing to known sequences
- **Quality Control**: Verifies that sequences belong to the same gene or region

## How to Use This Module

### Basic Command

```bash
python align_sequences.py <input_fasta> [output_directory]
```

### Example

Assuming you have filtered sequences from Module 01:

```bash
python modules/02_alignment/align_sequences.py results/passed_sequences.fasta results/
```

This will:
1. Take your FASTA file with multiple sequences
2. Align them using MAFFT (see section below)
3. Create output files in the `results/` directory

## Input Requirements

- **File Format**: FASTA format (`.fasta` or `.fa` extension)
- **Content**: Multiple DNA sequences (minimum 2 sequences required)
- **Structure**: Each sequence must have:
  - A header line starting with `>`
  - Sequence data on the following lines

Example input file:
```
>Species_A
ATGCGATCGATCGATCG
>Species_B
ATGCGATCGATCGTAGC
>Species_C
ATGCGATCGTATCGATCG
```

## Output Files

After running the alignment, you'll get three output files:

### 1. `aligned_sequences.fasta`
The main alignment file containing all sequences aligned with gaps:
```
>Species_A
ATG-CGATCGATCGATC--G
>Species_B
ATGCGATCGATCGTA--GC
>Species_C
ATGCGATCGTATCGATCG-
```

This is the file you'll use for phylogenetic analysis in Module 03.

### 2. `alignment_report.html`
A nice HTML report you can open in your web browser showing:
- Summary statistics (number of sequences, alignment length)
- Details for each sequence (original length, number of gaps, gap percentage)

Open it like any web page - just double-click the file!

### 3. `alignment_stats.json`
A JSON file with raw statistics for programmatic use. This is useful if you want to analyze the alignment data with code.

## Understanding Your Results

### What Are Gaps?

Gaps (represented as `-`) indicate positions where:
- A sequence is shorter than others
- There's an insertion or deletion (called "indels") between species
- The algorithm predicts nucleotides were added or removed during evolution

**Gap percentages** tell you how much of a sequence is gaps. Low percentages are good; very high percentages might indicate poor sequence quality.

### Conservation vs. Variation

As you look at your alignment:

- **Conservation (Conserved regions)**: Positions where most or all sequences have the same nucleotide
  - Example: All sequences have `A` at position 5
  - These regions are usually important for function
  - They're highly conserved through evolution

- **Variation (Variable regions)**: Positions where sequences differ
  - Example: Position 10 has `G` in Species A, `T` in Species B, `C` in Species C
  - These regions accumulate mutations over time
  - They're useful for distinguishing between species

Example alignment visualization:
```
Position:  1 2 3 4 5 6 7 8 9 10 11 12
Species_A: A T G C G A T C G A  T  C  (conserved)
Species_B: A T G C G T T C G C  T  C  (variation at 6,10)
Species_C: A T G - G A T C G A  T  C  (gap at 4)
           ^ ^ ^ ^ ^ ^ ^ ^ ^ ^  ^  ^
         conservation, high, high, mostly conserved region
```

## The MAFFT Tool

This module uses **MAFFT** (Multiple Alignment using Fast Fourier Transform) to align sequences.

### What Makes MAFFT Special?

- **Fast**: Efficiently aligns many sequences
- **Accurate**: Produces high-quality alignments
- **Smart**: The `--auto` option automatically chooses the best alignment strategy based on your data
- **Widely Used**: Standard tool in bioinformatics for alignment

You don't need to understand MAFFT's inner workings - just know it's the industry standard for multiple sequence alignment!

## Troubleshooting

### Problem: "MAFFT not found"

**Solution**: Install MAFFT first. Instructions depend on your system:

**On macOS** (with Homebrew):
```bash
brew install mafft
```

**On Linux** (Ubuntu/Debian):
```bash
sudo apt-get install mafft
```

**On any system** (Conda):
```bash
conda install -c bioconda mafft
```

### Problem: "Need at least 2 sequences for alignment"

**Solution**: Your input FASTA file has fewer than 2 sequences. Make sure you have multiple sequences in your file. Check by opening the file and counting the lines starting with `>`.

### Problem: The alignment looks strange or has very high gap percentages

**Possible causes**:
1. Sequences might be from different genes (not aligned to same region)
2. Sequences might have poor quality
3. Sequences might have different lengths by a lot

**Solutions**:
- Go back to Module 01 and check your quality control
- Make sure all sequences are from the same target gene
- Review the alignment in a text editor to see where the gaps are

### Problem: The script runs but produces no output

**Solution**: Check that the output directory has write permissions and that the input file is a valid FASTA format.

## Next Steps

Once you have your aligned sequences, you're ready for:

### Module 03: Phylogeny
Use your aligned sequences to build a phylogenetic tree! This will show the evolutionary relationships between your sequences.

**What you'll do**:
- Build a tree showing which sequences are most similar
- Understand evolutionary distances
- Create publication-quality tree figures

**Input**: Use the `aligned_sequences.fasta` file from this module!

## Tips for Success

1. **Always review your input**: Open the input FASTA file to understand your sequences before alignment
2. **Check the HTML report**: The alignment report helps you understand if something went wrong
3. **Keep backups**: Save your aligned sequences - you'll need them for phylogenetic analysis
4. **Document your work**: Note what parameters you used (even though `--auto` handles it for you)

## Quick Reference

| File | What it contains | Use for |
|------|-----------------|---------|
| Input FASTA | Multiple unaligned sequences | Run through alignment |
| `aligned_sequences.fasta` | Aligned sequences with gaps | Phylogenetic tree building |
| `alignment_report.html` | Summary and statistics | Quality checking |
| `alignment_stats.json` | Raw statistics | Advanced analysis |

## Questions?

If something doesn't work:
1. Check the error message carefully
2. Review the troubleshooting section above
3. Verify your input file has valid FASTA format
4. Check that you have at least 2 sequences

Good luck with your sequence alignment!
