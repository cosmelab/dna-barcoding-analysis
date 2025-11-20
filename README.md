# 🧬 DNA Barcoding Analysis

**ENTM201L Week 8 - Mosquito Species Identification**

Complete workflow for analyzing Sanger sequencing chromatograms (.ab1 files) to identify mosquito species using COI barcoding.

---

## 🎯 Quick Start for Students

### Prerequisites

1. **Docker Desktop** must be running on your computer
2. **Docker login** (required to pull the container):
   ```bash
   docker login
   ```
3. **Pull the analysis container**:
   ```bash
   docker pull cosmelab/dna-barcoding-analysis:latest
   ```

**Linux users:** You can use Podman instead of Docker Desktop:
```bash
sudo apt-get install podman
podman pull ghcr.io/cosmelab/dna-barcoding-analysis:latest
```

### Get the Repository

```bash
git clone https://github.com/cosmelab/dna-barcoding-analysis.git
cd dna-barcoding-analysis
```

---

## ⚠️ STEP 0: Complete the Tutorial FIRST (REQUIRED!)

**Before analyzing your own data, you MUST run the tutorial:**

```bash
./tutorial.sh
```

**This tutorial:**
- Uses test data (you can't break anything)
- Takes 15-20 minutes
- Teaches you all 5 steps of the workflow
- Shows you what results should look like
- Makes the actual assignment much easier

**DO NOT SKIP THIS!** Students who skip the tutorial get confused.

---

## 📋 The 5-Step Workflow

Once you've completed the tutorial, here's how to analyze YOUR mosquito sequences:

### Your Data

Put your `.ab1` chromatogram files in: `data/student_sequences/`

You should have pairs:
- Forward reads (ending in F): `AT-HV1F_A01.ab1`
- Reverse reads (ending in R): `AT-HV1R_B01.ab1`

### Step 1: Quality Control

Check which sequences are good enough to use:

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/01_quality_control/qc_chromatograms.py \
  data/student_sequences/ \
  results/my_analysis/qc/ \
  --open
```

**Look at the HTML report:**
- How many sequences passed QC?
- How many samples have BOTH F and R that passed?

### Step 2: Create Consensus Sequences

Combine forward and reverse reads:

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/02_consensus/create_consensus.py \
  results/my_analysis/qc/passed_sequences.fasta \
  results/my_analysis/consensus/ \
  --pairs-only \
  --open
```

The `--pairs-only` flag keeps only samples where BOTH F and R passed QC.

### Step 3: Combine with Reference Sequences

Add your sequences to the database of known Southern California mosquitoes:

```bash
cat results/my_analysis/consensus/consensus_sequences.fasta \
    data/reference_sequences/socal_mosquitoes.fasta \
    > results/my_analysis/consensus/combined_with_references.fasta
```

This creates a file with your sequences + 52 reference sequences.

### Step 4: Alignment and Phylogenetic Tree

**Part A: Align sequences**

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/03_alignment/align_sequences.py \
  results/my_analysis/consensus/combined_with_references.fasta \
  results/my_analysis/alignment/
```

**Part B: Build tree** (takes ~2-3 minutes)

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/04_phylogeny/build_tree.py \
  results/my_analysis/alignment/aligned_sequences.fasta \
  results/my_analysis/phylogeny/
```

Look at `results/my_analysis/phylogeny/tree.png` to see where your samples cluster!

### Step 5: Species Identification (BLAST)

Compare your sequences to GenBank database:

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/05_identification/identify_species.py \
  results/my_analysis/consensus/consensus_sequences.fasta \
  results/my_analysis/blast/
```

Look at the HTML report to see species matches and % identity scores.

---

## 📊 Understanding Your Results

### Quality Control
- **PASSED**: Sequence is good quality (>500bp, good Phred scores)
- **FAILED**: Too short, low quality, or unreadable

### Consensus Sequences
- Combines F+R reads for better accuracy
- Only keeps samples with complete pairs (both F and R passed)

### Phylogenetic Tree
- Shows evolutionary relationships
- Your samples cluster near related species
- Bootstrap values show confidence (higher = better)

### BLAST Results
- **>98% identity**: Same species
- **95-98% identity**: Same genus, possibly different species
- **<95% identity**: Different genus or poor quality

---

## 📁 Repository Structure

```
dna-barcoding-analysis/
├── data/
│   ├── student_sequences/       # PUT YOUR .ab1 FILES HERE
│   ├── test_data/               # Tutorial test data (8 .ab1 files)
│   └── reference_sequences/     # 52 SoCal mosquito COI sequences
├── results/
│   ├── tutorial/                # Tutorial results (from tutorial.sh)
│   └── my_analysis/             # Your real analysis results
├── modules/
│   ├── 01_quality_control/
│   ├── 02_consensus/
│   ├── 03_alignment/
│   ├── 04_phylogeny/
│   └── 05_identification/
├── docs/
│   ├── pipeline_workflow.md     # Visual guide to the workflow
│   └── iqtree_guide.md          # Understanding phylogenetic trees
├── tutorial.sh                  # INTERACTIVE TUTORIAL (RUN THIS FIRST!)
└── assignment.md                # Assignment worksheet
```

---

## 🛠️ Tools Included in Container

- **BioPython**: Chromatogram parsing and sequence handling
- **MAFFT**: Multiple sequence alignment
- **IQ-TREE2**: Maximum likelihood phylogenetic inference
- **BLAST+**: Species identification via NCBI
- **R packages**: ape, ggtree for tree visualization

---

## 📚 Reference Dataset

85 mosquito COI sequences from **Hoque et al. 2022**:
- 19 species from 6 genera
- Genera: *Aedes*, *Anopheles*, *Culex*, *Deinocerites*, *Psorophora*, *Uranotaenia*
- All sequences validated and published

**Citation**: Hoque MM, Valentine MJ, Kelly PJ, et al. Modification of the Folmer primers for the cytochrome c oxidase gene facilitates identification of mosquitoes. *Parasites Vectors*. 2022;15:437. doi:[10.1186/s13071-022-05494-2](https://doi.org/10.1186/s13071-022-05494-2)

---

## 🐛 Troubleshooting

### Docker Issues

**Problem**: "Cannot connect to the Docker daemon"
- **Solution**: Make sure Docker Desktop is running

**Problem**: "permission denied"
- **Solution**: Run `docker login` and enter your credentials

**Problem**: Container is slow
- **Solution**: Allocate more RAM to Docker (4GB+) in Docker Desktop settings

### Analysis Issues

**Problem**: No sequences pass QC
- **Solution**: Check chromatogram quality - may need re-sequencing

**Problem**: BLAST returns no hits
- **Solution**: Sequence may be contamination or poor quality

**Problem**: Tree has low bootstrap values
- **Solution**: Need more sequences or higher quality data

---

## 🆘 Need Help?

1. **Re-run the tutorial**: `./tutorial.sh`
2. **Read the visual workflow**: `docs/pipeline_workflow.md`
3. **Understand IQ-TREE**: `docs/iqtree_guide.md`
4. **Check assignment**: `assignment.md`
5. **Ask your TA or instructor**

---

## 🎓 Learning Goals

This workflow teaches you:
- How DNA barcoding identifies species
- Why quality control matters in sequencing
- How consensus sequences improve accuracy
- How to interpret phylogenetic trees
- How to use bioinformatics tools for real research

---

## 📜 License

**Code**: MIT License
**Educational Materials**: CC BY 4.0
**Reference Data**: See individual citations

---

## 🙏 Acknowledgments

- **UC Riverside Department of Entomology**
- **ENTM201L Students** (Fall 2025)
- **Hoque et al. 2022** for mosquito COI reference sequences
- **Tool developers**: BioPython, MAFFT, IQ-TREE, BLAST+

---

**Last Updated**: November 19, 2025
**Status**: Ready for Deployment
**Container**: `cosmelab/dna-barcoding-analysis:latest` (multi-architecture: amd64 + arm64)
