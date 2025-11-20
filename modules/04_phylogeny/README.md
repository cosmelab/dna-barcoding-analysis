# Module 03: Phylogenetic Tree Construction

## What is a Phylogenetic Tree?

A **phylogenetic tree** is a diagram that shows how different species (or DNA sequences) are related to each other. Think of it like a family tree, but for living organisms!

### Why is it useful for species identification?

When you sequence the DNA of an unknown organism, you can compare it to sequences from known species. A phylogenetic tree shows you:

- **Which known species is your unknown organism most closely related to?**
- **How similar are the DNA sequences?**
- **Could it be a new species, or a variation of an existing one?**

In the tree, the **closer two species are together (shorter branches between them), the more similar their DNA is**. This means they share a common ancestor more recently in evolutionary time.

---

## Quick Start: Building Your First Tree

### Command

To build a phylogenetic tree from aligned DNA sequences:

```bash
python modules/04_phylogeny/build_tree.py <aligned_sequences> [output_directory]
```

### Example

If you completed Module 02 (Sequence Alignment), you'll have an aligned FASTA file. Here's how to build the tree:

```bash
python modules/04_phylogeny/build_tree.py results/aligned_sequences.fasta results/
```

Or with a custom output directory:

```bash
python modules/04_phylogeny/build_tree.py data/my_alignment.fasta phylogeny_output/
```

---

## Inputs

### What you need to provide:

- **Aligned FASTA file** - This comes from Module 02 (Sequence Alignment)
  - The file must contain DNA sequences that are already aligned
  - All sequences must be the same length (aligned)
  - Example filename: `aligned_sequences.fasta`

### Example format:
```
>Species_A
ATGCTAGCTAG...
>Species_B
ATGCTAGCTAG...
>Unknown_Sample
ATGCTAGCTAG...
```

---

## Outputs

After running the script, you'll find these files in your output directory:

| File | Purpose |
|------|---------|
| **tree.treefile** | The phylogenetic tree in Newick format (for FigTree) |
| **tree.png** | A simple visualization of your tree |
| **phylogeny_report.html** | An HTML report you can open in your web browser |
| **tree.iqtree** | Detailed analysis results (what model was chosen, etc.) |
| **tree.log** | A log of the analysis run |

---

## Understanding Your Phylogenetic Tree

### How to read the tree:

```
          ┌─ Species_A
    ┌─────┤
    │     └─ Species_B
────┤
    │     ┌─ Unknown_Sample
    └─────┤
          └─ Species_C
```

In this example:
- **Species_A and Species_B are closely related** (they share a common branch point)
- **Unknown_Sample is less related to Species_A and B**, but closer to Species_C
- **All four species eventually connect** at a common ancestor (the base of the tree)

### What the branches mean:

- **Short branches** = Species are very similar (recently shared a common ancestor)
- **Long branches** = Species are quite different (ancestor was a long time ago)
- **Branch points (nodes)** = Where species split off from a common ancestor

### Bootstrap Support Values (numbers on branches)

The numbers you see on branches (in FigTree) represent **how confident we are** about that connection:
- **95-100** = Very confident this grouping is correct
- **70-94** = Moderately confident
- **Below 70** = Less confident (sequences might be ambiguous)

---

## How This Module Works

### IQ-TREE: The analysis tool

This module uses **IQ-TREE**, a powerful tool for building phylogenetic trees. Here's what it does:

1. **Auto-selects the best evolutionary model** (Model Finder Plus)
   - Different DNA sequences evolve at different rates
   - IQ-TREE finds the best mathematical model for YOUR data

2. **Builds the best-fitting tree** using Maximum Likelihood
   - Tests many possible tree arrangements
   - Keeps the one that best explains the sequence data

3. **Calculates bootstrap support** (1000 replicates)
   - A way to measure confidence in the tree
   - Tests how robust the tree is by re-sampling sequences

### Why this matters for species ID:

- **Accurate model** = Reliable tree that reflects true evolutionary relationships
- **Bootstrap values** = You know which connections to trust
- **Maximum Likelihood** = Statistically sound method, widely accepted in research

---

## Better Visualization: Using FigTree

The `tree.png` file is a quick preview, but **FigTree** gives you much better control over how the tree looks and shows bootstrap values clearly.

### Steps to visualize in FigTree:

1. **Download FigTree** (it's free!)
   - Visit: http://tree.bio.ed.ac.uk/software/figtree/
   - Download the version for your computer (Windows, Mac, or Linux)

2. **Open your tree file**
   - Launch FigTree
   - Go to File → Open
   - Select your `tree.treefile` from the output directory

3. **Enable bootstrap values**
   - Click the **"Node Labels"** checkbox on the left panel
   - You'll now see the confidence numbers on the branches!

4. **Customize your view**
   - Adjust line width, font size, colors
   - Rotate branches to see the tree from different angles
   - Zoom in/out to see details

---

## Troubleshooting

### Problem: "IQ-TREE not found"
**Solution:** IQ-TREE isn't installed on your computer.
```bash
# On macOS with Homebrew:
brew install iqtree

# On Linux (Ubuntu/Debian):
sudo apt-get install iqtree

# Or download from: http://www.iqtree.org/
```

### Problem: "Input file not found"
**Solution:** Check that your aligned FASTA file path is correct.
- Verify the file exists: `ls -la your_file.fasta`
- Use the correct path (relative or absolute)

### Problem: "All sequences are the same"
**Solution:** This can happen if sequences didn't align properly.
- Go back to Module 02 and check the alignment
- Make sure sequences are from the same gene region
- Verify that sequences are long enough (at least 100 bp recommended)

### Problem: Tree visualization looks strange
**Solutions:**
- Try opening in FigTree (better rendering)
- If tree is unrooted, you can root it in FigTree
- Very long sequences might take longer to visualize

### Problem: Extremely long running time
**Note:** On very large datasets (>100 sequences), this can take several minutes.
- This is normal! Tree construction is computationally intensive
- The script will display progress as it runs
- You can see `tree.log` to monitor progress

---

## What the Results Mean for Species Identification

### Scenario 1: Your unknown sequence is very close to a known species
```
Expected: >95% bootstrap support, short branch to known species
Conclusion: Likely the same species or a very close relative
```

### Scenario 2: Your unknown sequence is distant from known species
```
Expected: Separate branch position, >70 bootstrap support
Conclusion: Might be a new species, or species not in your reference database
```

### Scenario 3: Your unknown clusters within a species but is somewhat different
```
Expected: Groups with a species but on a longer sub-branch
Conclusion: Could be a subspecies, population variant, or individual variation
```

---

## Next Steps

Once you've built and examined your phylogenetic tree, you're ready for:

### **Module 04: Species Identification**
- Use your phylogenetic tree to identify your sequences
- Compare against reference databases
- Calculate similarity scores (BLAST, identity %)
- Generate a final species identification report

**Link to Module 04:** [See Module 04 documentation](../04_identification/README.md)

---

## Key Concepts Summary

| Concept | Meaning |
|---------|---------|
| **Phylogenetic tree** | Diagram showing evolutionary relationships |
| **Branch length** | How different two sequences are |
| **Bootstrap value** | Confidence in that part of the tree (0-100%) |
| **Newick format** | A text format for storing trees (like `.treefile`) |
| **Maximum Likelihood** | A statistical method for building trees |
| **Model selection** | Choosing the best evolutionary model for your data |

---

## Questions? Troubleshooting Help

- Check the `phylogeny_report.html` file for a visual summary
- Look at `tree.iqtree` for detailed technical information
- Read the `tree.log` for messages about what happened during analysis

---

**Ready to identify your species?** → Continue to Module 04: Species Identification
