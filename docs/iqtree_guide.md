# IQ-TREE Guide for DNA Barcoding

## What is IQ-TREE?

**IQ-TREE** (Efficient and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies) is a software tool that builds evolutionary trees (phylogenetic trees) from DNA sequences.

Think of it as a tool that:
- Compares DNA sequences from different organisms
- Figures out which organisms are most closely related
- Creates a family tree showing their evolutionary relationships

## Why Do We Use IQ-TREE?

IQ-TREE is the **gold standard** for phylogenetic analysis in 2025 because:

✓ **Fast** - Analyzes sequences much faster than older methods
✓ **Accurate** - Uses advanced statistical models to find the best tree
✓ **Automatic** - Picks the best evolutionary model for your data
✓ **Reliable** - Provides confidence scores (bootstrap values) for branches
✓ **Widely used** - Standard tool in scientific research

### Alternatives We Considered
- **RAxML** - Older, slower, but still good
- **FastTree** - Very fast but less accurate
- **MEGA** - GUI-based, good for teaching but less powerful
- **MrBayes** - Uses different method (Bayesian), very slow

**We chose IQ-TREE** because it gives the best balance of speed, accuracy, and ease of use for students.

## What is Maximum Likelihood (ML)?

**Maximum Likelihood** is a statistical method for building phylogenetic trees. It answers the question:

> "Given my DNA sequences, which tree arrangement is most likely to have produced these sequences through evolution?"

### Simple Explanation

Imagine you have DNA from 4 mosquito species. There are many ways to arrange them into a tree:

```
Tree 1:          Tree 2:          Tree 3:
  ┌─A              ┌─A              ┌─B
──┤               ──┤               ──┤
  │  ┌─B            │  ┌─C            │  ┌─A
  └──┤              └──┤              └──┤
     │  ┌─C            │  ┌─B            │  ┌─C
     └──┤              └──┤              └──┤
        └─D              └─D              └─D
```

**Maximum Likelihood asks:** "Which tree is most likely to be correct, given:
- The DNA differences we see
- A model of how DNA evolves (mutations, deletions, insertions)
- The laws of probability"

For each tree, ML calculates a **likelihood score** - how probable that tree is given your data. The tree with the highest score "wins."

### Why Maximum Likelihood Instead of Other Methods?

| Method | How It Works | Pros | Cons | When to Use |
|--------|--------------|------|------|-------------|
| **Maximum Likelihood (ML)** | Finds tree with highest probability given data | Most accurate; uses statistical model; handles complex evolution | Slower than simpler methods | Best for research & publication |
| **Neighbor-Joining (NJ)** | Joins closest sequences first | Very fast; simple | Less accurate; no statistical model | Quick previews, large datasets |
| **Maximum Parsimony (MP)** | Finds tree with fewest mutations | Fast; intuitive | Assumes all changes equally likely (wrong!) | Historical; rarely used now |
| **Bayesian (MrBayes)** | Samples many trees using probability | Most rigorous statistically | VERY slow (hours to days) | Research with critical questions |

**For DNA barcoding, Maximum Likelihood is the gold standard** because:
1. ✅ Accurate enough for species identification
2. ✅ Fast enough for student projects (1-5 minutes)
3. ✅ Uses a realistic model of DNA evolution
4. ✅ Provides statistical support (bootstrap values)
5. ✅ Widely accepted in scientific literature

### Why Not Bayesian?

Bayesian methods (like MrBayes) are theoretically superior BUT:
- **Time:** Takes hours to days for 50 sequences
- **Complexity:** Students must understand priors, MCMC, convergence
- **Overkill:** For species ID, ML is plenty accurate

**Our choice:** ML gives 95% of the accuracy in 1% of the time.

### In Practice

When IQ-TREE runs, you'll see:
```
Choosing best model... GTR+F+I+G4
Computing ML tree... done
Running 1000 bootstrap replicates... done
```

This means:
1. **Best model chosen:** GTR+F+I+G4 (complex model of DNA evolution)
2. **ML tree computed:** Found the most likely tree arrangement
3. **Bootstrap:** Tested how confident we are in each branch

## How IQ-TREE Works in Our Pipeline

### Input
IQ-TREE receives an **aligned FASTA file** from the alignment step:
```
>AT99_F
ATCGATCGATCG---ATCG
>AT99_R
ATCGATCGATCG---ATCG
>AT_ROCK_F
ATCGAT---ATCGATGATCG
>AT_ROCK_R
ATCGAT---ATCGATGATCG
```

All sequences must be the same length (gaps filled with `-`).

### What IQ-TREE Does

**Step 1: Choose the Best Model**
- DNA evolves in different ways (some positions change faster than others)
- IQ-TREE tests many models and picks the best one
- Common models: GTR+G, HKY+I, K2P

**Step 2: Build the Tree**
- Tries millions of possible tree arrangements
- Uses maximum likelihood to find the most probable tree
- This is like finding the family tree that best explains the DNA differences

**Step 3: Test the Tree (Bootstrap)**
- Creates 1000 "fake" datasets by resampling your data
- Builds a tree for each fake dataset
- Counts how often each branch appears
- This gives confidence scores (0-100%)

### Output
IQ-TREE creates several files:

```
results/phylogeny/
├── tree.treefile          # The best tree (Newick format)
├── tree.iqtree           # Detailed analysis report
├── tree.log              # What IQ-TREE did step-by-step
├── tree.png              # Simple tree visualization
└── phylogeny_report.html # Interactive HTML report
```

## Understanding the Tree

### Newick Format (tree.treefile)
```
((AT99_F:0.0123,AT99_R:0.0098)100:0.0456,(AT_ROCK_F:0.0234,AT_ROCK_R:0.0221)95:0.0567);
```

**Reading this:**
- `()` = groups organisms together
- `:0.0123` = branch length (amount of evolution)
- `100` or `95` = bootstrap support (% confidence)

### Tree Visualization

```
              ┌─ AT99_F
    ┌─100%───┤
    │         └─ AT99_R
────┤
    │         ┌─ AT_ROCK_F
    └─95%────┤
              └─ AT_ROCK_R
```

**Interpreting:**
- **Closer together** = more related (e.g., AT99_F and AT99_R are very similar)
- **Bootstrap 100%** = we're very confident in this grouping
- **Bootstrap <70%** = less confident, might not be real relationship
- **Longer branches** = more mutations = more evolution

## What the Numbers Mean

### Branch Length
- **0.001** = very few mutations (0.1% of positions different)
- **0.01** = moderate mutations (1% of positions different)
- **0.1** = many mutations (10% of positions different)

In DNA barcoding (COI gene):
- **Same species**: usually 0-2% difference (0.000-0.02)
- **Different species, same genus**: 2-10% difference (0.02-0.10)
- **Different genera**: >10% difference (>0.10)

### Bootstrap Support
- **90-100%** = Strong support - very confident
- **70-89%** = Moderate support - probably correct
- **50-69%** = Weak support - uncertain
- **<50%** = Very weak - might be wrong

**Rule of thumb:** Only trust branches with >70% bootstrap support.

## Reading the IQ-TREE Report (tree.iqtree)

When you open `tree.iqtree`, you'll see:

### 1. Model Selection
```
Best-fit model: GTR+F+I+G4
```
- **GTR** = General Time Reversible (allows different mutation rates)
- **+F** = Empirical base frequencies (uses your data's A/T/G/C ratios)
- **+I** = Invariable sites (some positions never change)
- **+G4** = Gamma rate heterogeneity (some positions change faster)

**What this means:** IQ-TREE chose a complex model that accounts for how DNA really evolves.

### 2. Likelihood Score
```
Log-likelihood: -2345.678
```
Higher (less negative) is better. Only useful for comparing trees from the same data.

### 3. Tree Statistics
```
Number of taxa: 52
Number of free parameters: 125
```
More sequences = more complex tree = needs more data to be confident.

### 4. Bootstrap Results
```
Branches with support >70%: 45/50 (90%)
```
If most branches have high support, your tree is reliable!

## Common Questions

### Q: Why does my tree look weird?
**A:** Check if:
- You have enough sequences (need at least 4, ideally 10+)
- Sequences are actually related (all from the same gene region)
- Alignment is good (run quality control first)
- You included reference sequences for comparison

### Q: What if bootstrap values are low (<70%)?
**A:** This means:
- Not enough data to be confident
- Sequences might be too similar (or too different)
- Need more sequences to resolve relationships

**Solution:** Add more reference sequences from GenBank.

### Q: Can I use the tree without reference sequences?
**A:** Yes, but it's limited:
- Shows relationships among YOUR sequences
- Can't identify species (need references for that)
- Can spot contamination (if one sequence is very different)

### Q: How long does IQ-TREE take?
**A:**
- 4-10 sequences: ~10 seconds
- 50 sequences: ~1 minute
- 100 sequences: ~5 minutes
- 1000 sequences: ~1 hour

Our pipeline uses 1000 ultrafast bootstrap replicates, which is much faster than traditional bootstrap (1000x replicates would take days!).

## Viewing Your Tree

### Option 1: FigTree (Recommended)
Download free from: http://tree.bio.ed.ac.uk/software/figtree/

**Steps:**
1. Open FigTree
2. File → Open → Select `tree.treefile`
3. Check "Node Labels" → "label" to see bootstrap values
4. Adjust appearance, colors, labels
5. Export as PDF or PNG

### Option 2: iTOL (Web-based)
Upload to: https://itol.embl.de/

**Steps:**
1. Upload `tree.treefile`
2. Customize colors and labels
3. Download publication-quality figure

### Option 3: Our HTML Report
Open `phylogeny_report.html` in your browser - includes basic visualization.

## Example: Interpreting Your Results

Let's say you sequenced mosquito samples and got this tree:

```
              ┌─ Your_Sample_1
    ┌─98%────┤
    │         └─ Aedes_aegypti_reference
────┤
    │         ┌─ Your_Sample_2
    └─95%────┤
              └─ Culex_quinquefasciatus_reference
```

**What this tells you:**
1. **Your_Sample_1** groups with *Aedes aegypti* (98% confidence)
   → Probably *Aedes aegypti*

2. **Your_Sample_2** groups with *Culex quinquefasciatus* (95% confidence)
   → Probably *Culex quinquefasciatus*

3. These two mosquitoes are from different genera
   → Confirms they're very different species

## Troubleshooting

### Error: "Alignment contains no parsimony-informative sites"
**Problem:** All sequences are identical (or only differ by gaps)
**Solution:** Check your alignment - you need sequences that actually differ

### Error: "Too few sequences"
**Problem:** IQ-TREE needs at least 4 sequences
**Solution:** Add more sequences or reference sequences

### Warning: "Some branch lengths are very long"
**Problem:** One sequence is very different from others
**Solution:**
- Might be contamination (check your sample)
- Might be from a different gene region
- Might be a different organism entirely

### Low bootstrap support
**Problem:** Can't confidently determine relationships
**Solution:**
- Add more sequences
- Use longer DNA sequences
- Check alignment quality
- Make sure sequences are actually related

## Further Reading

### Beginner-Friendly
- **IQ-TREE Tutorial**: http://www.iqtree.org/doc/Tutorial
- **Phylogenetics for Dummies**: Start with YouTube videos on "phylogenetic trees"

### Advanced
- **IQ-TREE Documentation**: http://www.iqtree.org/doc/
- **Model Selection Guide**: http://www.iqtree.org/doc/Substitution-Models
- **IQ-TREE Paper**: Nguyen et al. (2015) Molecular Biology and Evolution

### Our Course Materials
- See `tutorials/03_r_basics/` for tree visualization in R
- See `data/reference_sequences/` for example reference data

## Summary

**What you need to know:**
1. IQ-TREE builds evolutionary trees from DNA sequences
2. It automatically picks the best method for your data
3. Bootstrap values show confidence (>70% is good)
4. Branch lengths show how different sequences are
5. Use FigTree or iTOL to make pretty trees
6. Always include reference sequences for species ID

**What you DON'T need to worry about:**
- Complex mathematical models (IQ-TREE handles this)
- Manual parameter tuning (automatic mode works great)
- Programming (our pipeline runs it for you)

Just understand:
- What the tree shows (relationships)
- What the numbers mean (confidence and distance)
- How to interpret your results (which sequences are related)

Happy tree building! 🌳
