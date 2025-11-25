# Next Session - Start Here

## ✅ What's Ready:

1. **Container rebuilt and verified**
   - ✅ toytree 3.0.10 installed and working
   - ✅ phyTreeViz installed and working
   - ✅ tqdm available for progress bars
   - ✅ All previous tools still working

2. **Bugs fixed**
   - ✅ NC_ accession detection (NC_054318, NC_036006.1)
   - ✅ README paths (all use 01_qc, 02_consensus format)
   - ✅ Bootstrap values no longer colored
   - ✅ Legend simplified and positioned outside

3. **New features**
   - ✅ Visual ASCII progress bars with animations
   - ✅ Genus-based coloring (each genus has unique color)
   - ✅ Step-by-step progress indicators

## 🎯 Next Steps (in order):

### 1. Test toytree layouts (FIRST PRIORITY)

Test if toytree is good enough before adding ETE3:

```bash
# Generate tree with toytree to see quality
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 -c "
import toytree
from pathlib import Path

tree = toytree.tree('results/tutorial/04_phylogeny/tree.treefile')

# Try different layouts
for layout in ['r', 'c', 'd']:
    canvas = tree.draw(
        layout=layout,
        width=800,
        height=600,
        tip_labels_style={'font-size': 10}
    )
    canvas.save(f'results/tutorial/04_phylogeny/tree_toytree_{layout}.svg')
    print(f'✓ Generated {layout} layout')
"
```

Then open the SVG files and evaluate quality:
- Are the layouts clear?
- Are labels readable?
- Is circular tree useful?
- Do we need ETE3 or is toytree good enough?

### 2. If toytree is good → Implement in build_tree.py

Add toytree generation alongside Bio.Phylo:
- Keep Bio.Phylo for PNG (familiar format)
- Add toytree for SVG with circular layout
- Generate both in one run
- Embed both in HTML report

### 3. If toytree is not good → Add ETE3 for amd64

Modify Dockerfile:
```dockerfile
# Install ETE3 only for amd64
RUN if [ "$(uname -m)" = "x86_64" ]; then \
        pip3 install --no-cache-dir ete3; \
    fi
```

Add platform detection in build_tree.py:
```python
import platform

def generate_trees():
    if platform.machine() == 'x86_64':
        # Use ETE3 for multiple layouts
        use_ete3()
    else:
        # Use toytree on ARM64
        use_toytree()
```

### 4. Generate abbreviated trees

Run the abbreviation script and generate second set of trees:
```bash
python scripts/abbreviate_tree_names.py results/tutorial/04_phylogeny/tree.treefile
# Then regenerate trees with abbreviated names
```

### 5. Embed all trees in HTML

Update `generate_html_report()` to include:
- Full-name rectangular tree (PNG)
- Full-name circular tree (SVG if toytree, PNG if ETE3)
- Abbreviated rectangular tree (optional)
- Abbreviated circular tree (optional)
- Add tabs or toggle to switch between views

### 6. Test complete README workflow

Run all 5 steps from README:
```bash
# Step 1: QC
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/01_quality_control/qc_chromatograms.py \
  data/student_sequences/ \
  results/test_readme/01_qc/

# Step 2: Consensus
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/02_consensus/create_consensus.py \
  results/test_readme/01_qc/passed_sequences.fasta \
  results/test_readme/02_consensus/ \
  --pairs-only

# Step 3: Combine with references
cat results/test_readme/02_consensus/consensus_sequences.fasta \
    data/reference_sequences/socal_mosquitoes.fasta \
    > results/test_readme/02_consensus/combined_with_references.fasta

# Step 4: Alignment
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/03_alignment/align_sequences.py \
  results/test_readme/02_consensus/combined_with_references.fasta \
  results/test_readme/03_alignment/

# Step 4B: Phylogeny
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/04_phylogeny/build_tree.py \
  results/test_readme/03_alignment/aligned_sequences.fasta \
  results/test_readme/04_phylogeny/

# Step 5: BLAST
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/05_identification/identify_species.py \
  results/test_readme/02_consensus/consensus_sequences.fasta \
  results/test_readme/05_blast/
```

Verify all steps complete without errors and output files are created.

## 📋 Key Decisions to Make:

1. **Is toytree good enough?**
   - YES → Skip ETE3, use toytree for all layouts
   - NO → Add ETE3 for amd64 with architecture detection

2. **Which layouts to include?**
   - Minimum: Rectangular + Circular
   - Optional: Radial, Fan, Unrooted (if using ETE3)

3. **Do we need abbreviated trees?**
   - Test if labels are too crowded
   - If yes, generate both full and abbreviated
   - If no, skip abbreviation

## 🎨 Expected Output:

After implementing, students will see:
- Beautiful progress bars during analysis
- Multiple tree layouts in HTML report
- Genus-colored trees with clear legend
- Option to view abbreviated names if needed
- All working from simple README commands

## 📝 Files to Modify:

- `modules/04_phylogeny/build_tree.py` - Add toytree/ETE3 generation
- `container/Dockerfile` - Conditionally add ETE3 for amd64 (if needed)
- `README.md` - Update with new tree layout info (if needed)

## 🔗 Useful Links:

- Container build: https://github.com/cosmelab/dna-barcoding-analysis/actions
- Latest commit: 7ca7ee8
- toytree docs: https://toytree.readthedocs.io/
- ETE3 docs: http://etetoolkit.org/

---

**Safe to start fresh!** Everything is committed and pushed. Container is ready with toytree.
