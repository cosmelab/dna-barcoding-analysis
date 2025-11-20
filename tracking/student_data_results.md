# DNA Barcoding Analysis - Complete Workflow Summary

## Updated Workflow (5 Steps)

### ✅ STEP 1: Quality Control
**Command:**
```bash
python3 modules/01_quality_control/qc_chromatograms.py data/student_sequences/ results/student_qc/
```

**Results:**
- Input: 30 .ab1 files (15 samples, each with F and R)
- Passed QC: **12 sequences**
- Failed QC: 18 sequences
- Output: `results/student_qc/passed_sequences.fasta`

---

### ✅ STEP 2: Consensus Sequences
**Command:**
```bash
python3 modules/02_consensus/create_consensus.py results/student_qc/passed_sequences.fasta results/student_consensus/ --pairs-only
```

**Results:**
- Input: 12 passed sequences (mixed F and R)
- Consensus created (F+R pairs): **4 samples**
  - ✓ AT-HV1
  - ✓ AT-HV3
  - ✓ AT-JM2
  - ✓ AT-WL2
- Skipped (missing F or R): 4 samples
- Output: `results/student_consensus/consensus_sequences.fasta`

---

### ✅ STEP 3: Combine with References
**Command:**
```bash
cat results/student_consensus/consensus_sequences.fasta data/reference_sequences/socal_mosquitoes.fasta > results/student_consensus/combined_with_references.fasta
```

**Results:**
- Student consensus: 4 sequences
- Reference sequences: 52 sequences (known SoCal mosquitoes)
- **Total: 56 sequences**
- Output: `results/student_consensus/combined_with_references.fasta`

---

### ✅ STEP 4: Alignment & Phylogenetic Tree
**Commands:**
```bash
# Alignment
python3 modules/03_alignment/align_sequences.py results/student_consensus/combined_with_references.fasta results/student_alignment/

# Tree
python3 modules/04_phylogeny/build_tree.py results/student_alignment/aligned_sequences.fasta results/student_phylogeny/
```

**Results:**
- Aligned: 56 sequences (4 students + 52 references)
- Alignment length: 2435 bp
- Tree: Shows where student samples cluster with known species
- Outputs:
  - `results/student_alignment/alignment_report.html`
  - `results/student_phylogeny/tree.png`
  - `results/student_phylogeny/tree.treefile` (for FigTree)

---

### ✅ STEP 5: Species Identification (BLAST)
**Command:**
```bash
python3 modules/05_identification/identify_species.py results/student_consensus/consensus_sequences.fasta results/student_blast/
```

**Results:**

| Sample | Species Identified | % Identity | Common Name |
|--------|-------------------|------------|-------------|
| AT-HV1 | *Aedes albopictus* | 99.55% | Asian tiger mosquito |
| AT-HV3 | *Culex pipiens* | 98.12% | Northern house mosquito |
| AT-JM2 | *Culex pipiens* | 99.25% | Northern house mosquito |
| AT-WL2 | *Culex pipiens* | 98.67% | Northern house mosquito |

- Output: `results/student_blast/identification_report.html`

---

## Summary Statistics

### Sample Success Rate
- Total samples submitted: 15 (30 .ab1 files)
- Samples with both F & R passing QC: **4 (26.7%)**
- Species identified: **4**
  - 1 *Aedes albopictus* (Asian tiger mosquito)
  - 3 *Culex pipiens* (Northern house mosquito)

### Why Some Samples Failed
1. **Low quality scores** - Poor chromatogram quality
2. **Short sequences** - Less than minimum length threshold
3. **Missing F or R** - Only one direction passed QC

### Key Insight
Having both forward AND reverse reads pass QC is critical for creating high-quality consensus sequences. Only 4/15 samples (26.7%) had both F and R pass QC in this dataset.

---

## Files Created

```
results/
├── student_qc/
│   ├── qc_report.html                    # QC results with chromatograms
│   ├── qc_results.json
│   └── passed_sequences.fasta            # 12 sequences that passed
│
├── student_consensus/
│   ├── consensus_report.html             # Consensus pairing results
│   ├── consensus_sequences.fasta         # 4 consensus sequences
│   └── combined_with_references.fasta    # 56 sequences (4+52)
│
├── student_alignment/
│   ├── alignment_report.html             # Visual alignment viewer
│   ├── aligned_sequences.fasta           # 56 aligned sequences
│   └── alignment_stats.json
│
├── student_phylogeny/
│   ├── phylogeny_report.html             # Tree visualization
│   ├── tree.png                          # Tree image
│   ├── tree.treefile                     # Newick format (for FigTree)
│   └── tree.iqtree                       # Full IQ-TREE analysis
│
└── student_blast/
    ├── identification_report.html         # Species ID results
    └── identification_results.json        # Raw BLAST data
```

---

## Workflow Improvements Made

1. **Consensus Sequences** - Now combines F and R reads for better accuracy
2. **--pairs-only Flag** - Only uses samples with complete F+R pairs
3. **Reference Integration** - Tree includes 52 known SoCal mosquito sequences
4. **Updated Documentation** - Visual ASCII workflow shows all 5 steps
5. **Multi-Architecture Container** - Works on Intel and Apple Silicon Macs

---

## Next Steps for Students

1. Run QC on their .ab1 files
2. Create consensus sequences (--pairs-only)
3. Combine with references
4. Build alignment and tree
5. Run BLAST for species ID
6. Fill in assignment table with results
7. Interpret tree to see how samples cluster with known species

