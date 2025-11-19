# Phylogeny Scripts - Quick Reference Guide

## One-Line Summaries

| Script | Purpose | Main Concept |
|--------|---------|--------------|
| `run_iqtree.sh` | Build phylogenetic tree | Maximum likelihood inference |
| `visualize_tree.py` | Draw tree as image | Tree topology & bootstrap values |
| `root_tree.py` | Add directionality to tree | Evolutionary ancestor-descendant |
| `calculate_distances.py` | Compute pairwise distances | Species identification threshold |

## Quick Commands

### Build a Tree
```bash
./run_iqtree.sh aligned_sequences.fasta
```

### View the Tree
```bash
python visualize_tree.py aligned_sequences_iqtree.treefile -o tree.png
```

### Root the Tree
```bash
python root_tree.py aligned_sequences_iqtree.treefile midpoint -o rooted.treefile
```

### Calculate Distances
```bash
python calculate_distances.py -a aligned_sequences.fasta --threshold-analysis
```

## Common Questions

### "My sequences aren't aligned yet"
1. Use MAFFT first: `mafft sequences.fasta > aligned.fasta`
2. Then run: `./run_iqtree.sh aligned.fasta`

### "How do I know if my tree is good?"
- Check bootstrap values (>70% = good)
- Read the .iqtree file for model fit statistics
- If most values <70%, add more sequences or use different gene

### "What does bootstrap value mean?"
- Percentage of times this group appears in resampled data
- >95% = certain this group exists
- <70% = uncertain (may not be real grouping)

### "Which outgroup should I use?"
- Species known to be distantly related to your study species
- Often a different genus or family
- List all species first: `python root_tree.py tree.treefile --list-species`

### "How different are same species supposed to be?"
- Same species: usually <0.003 (0.3%)
- Different species: usually >0.03 (3%)
- This is the "barcoding gap"

### "My tree looks too complex"
- Too many sequences? Subset to representatives
- Poor alignment? Recheck with MAFFT
- Some genes have low variation? Choose different gene

## Parameter Guide

### run_iqtree.sh
```bash
./run_iqtree.sh input.fasta [options]

-m <model>    Substitution model (default: MFP auto-select)
-bb <number>  Bootstrap replicates (default: 1000)
              500-1000 = fast (learning)
              2000-5000 = publication quality
-nt <number>  CPU threads to use (default: AUTO uses all)
```

### visualize_tree.py
```bash
python visualize_tree.py tree.treefile [options]

-o <file>     Save to file (png, pdf, jpg)
-w <number>   Figure width in inches (default: 10)
-h <number>   Figure height in inches (default: 8)
--dpi <number> Resolution (default: 100)
              Use 300 for publication-quality
```

### root_tree.py
```bash
python root_tree.py tree.treefile [method] [options]

method:
  midpoint      Root at midpoint (use if no outgroup)
  outgroup      Root using outgroup species

-g <name>     Outgroup species (use multiple for multiple species)
-o <file>     Save rooted tree
--list-species Show all sequences in tree
```

### calculate_distances.py
```bash
python calculate_distances.py [input] [options]

Input (choose one):
  -a <file>     Alignment FASTA
  -t <file>     Tree file

Options:
  -m <method>   Distance method
                p-distance (simple)
                jukes-cantor (corrected)
                kimura-2p (best for DNA, default)
  -o <file>     Save to CSV
  --threshold-analysis  Show species identification guidance
```

## Key Phylogenetic Terms

| Term | Meaning |
|------|---------|
| **Bootstrap** | Confidence in tree grouping (0-100%) |
| **Branch length** | Evolutionary distance (substitutions/site) |
| **Clade** | A group of related sequences |
| **Root** | Common ancestor of all sequences |
| **Topology** | The branching structure of the tree |
| **Model** | Mathematical description of DNA evolution |
| **Intraspecific** | Within a species (genetic variation) |
| **Interspecific** | Between different species |

## Visual Tree Guide

```
                    ┌─ Species_A ──┐
              ┌─────┤              │ 95% bootstrap
        ┌─────┤ 95  └─ Species_B ──┘  (strong support)
    ────┤     └─────────────────────── Species_C
        │
        │            ┌─ Species_D ──┐
        └────────────┤              │ 65% bootstrap
                 65  └─ Species_E ──┘  (weak support)

        ← Branch length →
        = evolutionary distance
```

## Data Flow

```
Your aligned FASTA file
        ↓
    ./run_iqtree.sh
        ↓
    .treefile (best ML tree)
        ├→ python visualize_tree.py (make figure)
        ├→ python root_tree.py (add direction)
        └→ python calculate_distances.py (identify species)
```

## Tips

- Start with small datasets (3-5 species) to learn
- Use default parameters for learning, adjust for publication
- Always check the .iqtree file for details
- Bootstrap >70% is good, >95% is very good
- Different genes have different variation rates

## Help Commands

```bash
# See all options for each script
./run_iqtree.sh -h
python visualize_tree.py -h
python root_tree.py -h
python calculate_distances.py -h

# Learn concepts
python visualize_tree.py --concepts
python root_tree.py --concepts
python calculate_distances.py --concepts
```

## Typical Workflow

```bash
# 1. Align sequences (use MAFFT separately)
mafft sequences.fasta > aligned.fasta

# 2. Build tree
./run_iqtree.sh aligned.fasta

# 3. Visualize
python visualize_tree.py aligned_iqtree.treefile -o my_tree.png

# 4. Root tree
python root_tree.py aligned_iqtree.treefile outgroup \
    --outgroup "Known_outgroup" -o rooted.treefile

# 5. Identify species
python calculate_distances.py -a aligned.fasta \
    --threshold-analysis -o distances.csv

# 6. Use distances to identify unknown specimens
# Distance < 0.03 to reference = same species!
```

## Interpreting Bootstrap Values

```
99-100%  ████████ Perfect! This group definitely exists
95-99%   ███████░ Excellent! Very confident
85-94%   ██████░░ Good! Fairly confident
70-84%   █████░░░ Moderate, probably okay
50-69%   ████░░░░ Weak, be cautious
<50%     ███░░░░░ Very weak, may not be real
```

## Interpreting Distance Values

```
0.000-0.003  ▌  Same individual or very recent divergence
0.003-0.01   ▌▌ Same species (normal intraspecific variation)
0.01-0.03    ▌▌▌ Between closely related species
0.03-0.1     ▌▌▌▌ Different species
0.1-0.5      ▌▌▌▌▌ Distant species
>0.5         ▌▌▌▌▌▌ Very distantly related
```

## File Extensions Explained

| Extension | Meaning |
|-----------|---------|
| `.fasta` | Sequence alignment (text, all same length) |
| `.treefile` | Phylogenetic tree (Newick format) |
| `.iqtree` | IQ-TREE analysis report (human-readable) |
| `.log` | Run log (mostly debug info) |
| `.ckp.gz` | Checkpoint (for resuming) |
| `.contree` | Consensus tree (from bootstrap) |

## Common Problems & Solutions

| Problem | Likely Cause | Solution |
|---------|--------------|----------|
| "IQ-TREE not found" | Not installed | `conda install iqtree` |
| Tree is unresolved (looks like star) | Too few sequences or similar sequences | Add more diverse sequences |
| Bootstrap values all <50% | Bad alignment or wrong model | Check alignment, more sequences |
| Cannot visualize tree | Too many sequences | Subset to key representatives |
| Root appears to be wrong | Chose wrong outgroup | Check outgroup selection |

## Learning Path

1. **Start here**: Run basic tree with example data
2. **Visualize**: Make a PNG figure and understand the structure
3. **Understand bootstrap**: Read tutorial on why >70% is good
4. **Rooting**: Learn about midpoint vs outgroup rooting
5. **Distances**: Use for species identification
6. **Advanced**: Adjust parameters, try different models

---

**For questions**: See README.md for detailed explanations of each script
