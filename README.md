<div align="center">

# 🧬 DNA Barcoding Analysis Pipeline

### From Chromatogram to Species Discovery

![Status](https://img.shields.io/badge/Status-Production_Ready-brightgreen?style=for-the-badge)
![Course](https://img.shields.io/badge/Course-ENTM201L-bd93f9?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-8be9fd?style=for-the-badge)
![Docker](https://img.shields.io/badge/Docker-Multi--Arch-50fa7b?style=for-the-badge)

**[🚀 Quick Start](start_here.md)** | **[📖 Assignment](assignment.md)** | **[📚 Docs](docs/)** | **[🐳 Container](https://hub.docker.com/r/cosmelab/dna-barcoding-analysis)**

</div>

---

## 🎯 Overview

A complete automated workflow for analyzing Sanger sequencing chromatograms (.ab1 files) to identify mosquito species using **COI DNA barcoding**. Built for ENTM201L students at **UC Riverside** with zero coding experience required.

<table>
<tr>
<td><strong>🏛️ Institution</strong></td>
<td>University of California, Riverside</td>
</tr>
<tr>
<td><strong>📖 Course</strong></td>
<td>ENTM201L - Molecular Biology Laboratory</td>
</tr>
<tr>
<td><strong>👨‍🎓 Target Users</strong></td>
<td>Undergraduate students (zero coding experience)</td>
</tr>
<tr>
<td><strong>⏱️ Analysis Time</strong></td>
<td>~5 minutes (all 5 steps automated)</td>
</tr>
<tr>
<td><strong>🧑‍🏫 Instructor</strong></td>
<td>Luciano Cosme, Department of Entomology</td>
</tr>
</table>

---

## 🚀 Quick Start

### 👨‍🎓 For Students

**→ See [start_here.md](start_here.md) for complete beginner's guide**

```bash
# STEP 1: Learn with test data (15 min)
./tutorial.sh

# STEP 2: Analyze YOUR mosquito sequences (5 min)
./run-analysis.sh

# STEP 3: Fill out assignment.md
```

### 👨‍🏫 For Instructors

```bash
# Clone as GitHub Classroom template
git clone https://github.com/cosmelab/dna-barcoding-analysis.git

# Students get their own repos:
# github.com/cosmelab/dna-barcoding-analysis-STUDENT-USERNAME
```

**GitHub Classroom Ready** ✅ Use as template repository

---

## 🔬 The 5-Step Pipeline

<details>
<summary><strong>📊 Step 1: Quality Control</strong></summary>
<br>

**What it does:**
- Analyzes .ab1 chromatogram files
- Checks quality scores (Phred Q30+ required)
- Validates sequence length (>500bp required)
- Filters out low-quality reads

**Output:**
- `qc_report.html` - Interactive quality control report
- `passed_sequences.fasta` - High-quality sequences only

**Why it matters:** Garbage in = garbage out. Bad sequences produce unreliable species IDs.

</details>

<details>
<summary><strong>🧬 Step 2: Consensus Sequences</strong></summary>
<br>

**What it does:**
- Pairs forward (F) and reverse (R) reads
- Reverse-complements the R read
- Creates consensus sequence from F+R alignment
- Filters for complete pairs only

**Output:**
- `consensus_sequences.fasta` - Final consensus sequences
- `consensus_report.html` - Alignment visualization

**Why it matters:** Combining F+R reads doubles coverage and accuracy.

</details>

<details>
<summary><strong>📚 Step 3: Combine with References</strong></summary>
<br>

**What it does:**
- Adds your sequences to 52 reference mosquito COI sequences
- References include 19 species from 6 genera (*Aedes*, *Anopheles*, *Culex*, *Deinocerites*, *Psorophora*, *Uranotaenia*)
- All references trimmed to ~700bp barcode region

**Output:**
- `combined_with_references.fasta` - Your sequences + references

**Why it matters:** Can't build a tree without known species for context.

</details>

<details>
<summary><strong>🌳 Step 4: Phylogenetic Tree</strong></summary>
<br>

**What it does:**
- Aligns all sequences with MAFFT
- Builds maximum likelihood tree with IQ-TREE2
- Calculates 1000 ultrafast bootstrap replicates
- Generates 4 tree layouts (rectangular, circular, unrooted, radial)

**Output:**
- `tree.png`, `tree_circular.pdf`, etc. - Tree visualizations
- `phylogeny_report.html` - Interactive tree explorer

**Why it matters:** Shows evolutionary relationships. Your samples cluster with related species.

</details>

<details>
<summary><strong>🔍 Step 5: Species Identification</strong></summary>
<br>

**What it does:**
- BLASTs your sequences against NCBI GenBank
- Returns top 10 matches with % identity
- Interprets results (>98% = same species)

**Output:**
- `identification_report.html` - BLAST results table
- Top hits with accession numbers and % identity

**Why it matters:** Confirms species ID from tree with global database.

</details>

---

## 📦 What's Included

### 🛠️ Bioinformatics Tools
- **BioPython** - Chromatogram parsing and sequence handling
- **MAFFT** - Multiple sequence alignment (industry standard)
- **IQ-TREE2** - Maximum likelihood phylogenetic inference
- **BLAST+** - Species identification via NCBI GenBank
- **toytree** - Beautiful tree visualizations with genus coloring

### 🎨 Interactive HTML Reports
- Quality control dashboard with chromatogram viewer
- Consensus sequence comparisons (F vs R alignment)
- Alignment heatmaps (conservation visualization)
- Phylogenetic trees (4 layouts, genus-colored)
- BLAST results tables (sortable, interactive)

### 🎨 Beautiful Terminal
- **Zsh** with oh-my-zsh framework
- **Dracula theme** - professional dark colors
- **Git integration** - see status in prompt
- **Aliases**: `ll` (detailed view), `lt` (tree view)

### 📚 Reference Dataset
- 52 Southern California mosquito COI sequences
- 19 species from 6 genera
- All trimmed to ~700bp barcode region
- Published sequences from Hoque et al. 2022

---

## 📂 Repository Structure

```
dna-barcoding-analysis/
├── 📄 start_here.md              # 👈 START HERE! Complete beginner's guide
├── 📄 assignment.md              # Student assignment questions
├── 🚀 tutorial.sh                # Step 1: Learn with test data
├── 🚀 run-analysis.sh            # Step 2: Analyze YOUR data
│
├── 📁 data/
│   ├── student_sequences/        # PUT YOUR .ab1 FILES HERE
│   ├── test_data/                # 8 test chromatograms (for tutorial)
│   └── reference_sequences/      # 52 known mosquito sequences
│
├── 📁 results/
│   ├── tutorial/                 # Tutorial output (test data)
│   └── my_analysis/              # YOUR analysis output
│       ├── 01_qc/                # Quality control results
│       ├── 02_consensus/         # Consensus sequences
│       ├── 03_alignment/         # MAFFT alignment
│       ├── 04_phylogeny/         # Trees (4 layouts)
│       └── 05_blast/             # Species identification
│
├── 📁 modules/                   # Python analysis scripts
│   ├── 01_quality_control/
│   ├── 02_consensus/
│   ├── 03_alignment/
│   ├── 04_phylogeny/
│   └── 05_identification/
│
├── 📁 docs/                      # Documentation
│   ├── pipeline_workflow.md
│   ├── iqtree_guide.md
│   └── reference_trimming.md
│
└── 📁 intro_to_cli/              # Optional CLI tutorials
```

---

## 💻 System Requirements

### ✅ Required Software

<table>
<tr>
<td><strong>🐋 Docker Desktop</strong></td>
<td>Windows 10+, macOS 10.15+, or Linux</td>
</tr>
<tr>
<td><strong>🐙 Git</strong></td>
<td>For cloning the repository</td>
</tr>
<tr>
<td><strong>🐳 Docker Hub Account</strong></td>
<td>Free account (no payment needed)</td>
</tr>
</table>

### 💡 Recommended (Optional)

- **VS Code** - Best experience with integrated terminal
- **4GB+ RAM** allocated to Docker for faster tree building
- **Internet connection** - For BLAST searches

### 🖥️ Platform Support

| Platform | Status | Notes |
|----------|--------|-------|
| 🍎 **macOS** (Intel) | ✅ Fully Supported | Native amd64 |
| 🍎 **macOS** (Apple Silicon) | ✅ Fully Supported | Native arm64 |
| 🪟 **Windows 10/11** | ✅ Fully Supported | Requires WSL2 |
| 🐧 **Linux** | ✅ Fully Supported | Native support |

**Multi-architecture container:** Automatically uses correct version for your system!

---

## 🎓 Learning Outcomes

Upon completing this workflow, students will be able to:

✅ **Assess DNA sequence quality** from chromatogram data
✅ **Interpret quality metrics** (Phred scores, base calling)
✅ **Understand consensus sequences** and why F+R reads matter
✅ **Read phylogenetic trees** and identify evolutionary relationships
✅ **Perform species identification** using BLAST and % identity
✅ **Use Docker containers** for reproducible bioinformatics
✅ **Navigate the command line** with confidence

---

## 🔄 Adaptable for Other Projects

**This pipeline is generic!** Use it for any Sanger sequencing project:

### 🦋 Different Organisms
- Insects, plants, fungi, bacteria, fish, mammals
- Any organism with reference sequences in GenBank

### 🧬 Different Barcode Regions
- **COI** (animals) - current default
- **ITS** (fungi)
- **rbcL, matK** (plants)
- **16S rRNA** (bacteria)

### 🔧 How to Customize

```bash
# 1. Replace reference sequences
python3 data/reference_sequences/trim_references_to_barcode.py \
  your_genbank_refs.fasta your_trimmed_refs.fasta --start 50 --end 750

# 2. Replace student sequences (put your .ab1 files here)
cp ~/my_chromatograms/*.ab1 data/student_sequences/

# 3. Run the pipeline (same workflow!)
./run-analysis.sh
```

**BLAST automatically searches NCBI** for any organism - species ID works for everything!

---

## 🐛 Troubleshooting

<details>
<summary><strong>🐳 Docker Issues</strong></summary>
<br>

**"Cannot connect to Docker daemon"**
- ✅ Make sure Docker Desktop is running
- ✅ Check system tray/menu bar for Docker icon

**"Permission denied"**
- ✅ Run `docker login` with your Docker Hub credentials
- ✅ Windows: Make sure WSL2 integration is enabled in Docker settings

**Container is slow**
- ✅ Allocate more RAM to Docker (Settings → Resources)
- ✅ Recommended: 4GB+ for tree building

</details>

<details>
<summary><strong>📊 Analysis Issues</strong></summary>
<br>

**No sequences pass QC**
- ✅ Check chromatogram quality - may need re-sequencing
- ✅ Look at the QC report HTML to see failure reasons

**BLAST returns no hits**
- ✅ Sequence may be contamination or very poor quality
- ✅ Check alignment - might be wrong reading frame

**Tree has low bootstrap values**
- ✅ Normal for closely related species
- ✅ Add more reference sequences for better resolution

</details>

<details>
<summary><strong>🪟 Windows/WSL Issues</strong></summary>
<br>

**Enable Virtualization in BIOS**
- ✅ Restart → Enter BIOS (F2, F10, Del, or Esc)
- ✅ Enable "Intel Virtualization Technology" or "AMD-V"

**Docker Commands Hang**
- ✅ Restart WSL: `wsl --shutdown` in PowerShell (as Admin)
- ✅ Reopen WSL terminal and try again

**Verify Docker Works**
```bash
docker run hello-world
```

</details>

---

## 🆘 Getting Help

1. 📖 **Read [start_here.md](start_here.md)** - Complete beginner's guide
2. 🔬 **Check [docs/pipeline_workflow.md](docs/pipeline_workflow.md)** - Visual workflow
3. 🌳 **Read [docs/iqtree_guide.md](docs/iqtree_guide.md)** - Understanding trees
4. 🎓 **Ask your TA or instructor** - Office hours available

---

## 📊 Pipeline Statistics

| Metric | Value |
|--------|-------|
| ⏱️ **Tutorial Time** | ~3 minutes (all 5 steps) |
| ⚡ **Analysis Time** | ~5 minutes (automated) |
| 🧬 **Reference Sequences** | 52 mosquito COI sequences |
| 🦟 **Species Covered** | 19 species from 6 genera |
| 🌳 **Tree Layouts** | 4 visualizations (rectangular, circular, unrooted, radial) |
| 📊 **HTML Reports** | 5 interactive dashboards |
| 🐳 **Container Size** | ~2.5GB (includes all tools) |

---

## 🔬 Citation

If you use this pipeline in your research or teaching, please cite:

**Reference sequences:**
> Hoque MM, Valentine MJ, Kelly PJ, et al. Modification of the Folmer primers for the cytochrome c oxidase gene facilitates identification of mosquitoes. *Parasites Vectors*. 2022;15:437. [doi:10.1186/s13071-022-05494-2](https://doi.org/10.1186/s13071-022-05494-2)

---

## 📜 License

<div align="center">

![MIT License](https://img.shields.io/badge/Code-MIT-brightgreen?style=for-the-badge)
![CC BY 4.0](https://img.shields.io/badge/Docs-CC_BY_4.0-8be9fd?style=for-the-badge)

**Code:** MIT License
**Educational Materials:** Creative Commons Attribution 4.0 (CC BY 4.0)
**Reference Data:** See individual citations

**You are free to:**
- ✅ Use for teaching and research
- ✅ Modify and adapt for your needs
- ✅ Share with attribution

</div>

---

## 🙏 Acknowledgments

- **UC Riverside Department of Entomology**
- **ENTM201L Students** (Fall 2025)
- **Hoque et al. 2022** for Southern California mosquito COI sequences
- **Open-source developers**: BioPython, MAFFT, IQ-TREE, BLAST+, toytree teams

---

<div align="center">

## 🔗 Quick Links

**[🚀 Get Started](start_here.md)** | **[📖 Assignment](assignment.md)** | **[🐳 Docker Hub](https://hub.docker.com/r/cosmelab/dna-barcoding-analysis)** | **[📚 Course Website](https://cosmelab.github.io/entm201l-fall2025/)**

---

**Last Updated**: November 25, 2025
**Status**: ✅ Production Ready - Student Tested
**Container**: `cosmelab/dna-barcoding-analysis:latest` (multi-arch: amd64 + arm64)
**GitHub Classroom**: ✅ Template Ready

**Instructor**: Luciano Cosme | Department of Entomology | UC Riverside

![UCR Entomology](https://img.shields.io/badge/UCR-Entomology-FFB81C?style=for-the-badge)

</div>
