# 🧬 DNA Barcoding Analysis

**ENTM201L Week 8 - Mosquito Species Identification**

Complete workflow for analyzing Sanger sequencing chromatograms (.ab1 files) to identify mosquito species using COI barcoding.

---

## ⚠️ PREREQUISITES (Complete BEFORE Starting)

**You must have these installed and configured:**

1. **Docker Desktop** - Must be running before you start
   - See: [Docker Setup Guide](https://cosmelab.github.io/entm201l-fall2025/setup/cli-tools.html)

2. **Docker Hub Account** - Required to pull the analysis container
   - **Create account**: Go to [https://hub.docker.com/signup](https://hub.docker.com/signup)
   - **Free account is fine** - no payment needed
   - Remember your username and password!

3. **VS Code** - **HIGHLY RECOMMENDED** for this project
   - See: [CLI Tools Setup](https://cosmelab.github.io/entm201l-fall2025/setup/cli-tools.html)
   - Includes integrated terminal, Docker support, Dev Containers extension
   - Makes running commands much easier than switching between apps

4. **GitHub Account** - For accessing your assignment
   - See: [GitHub Setup Guide](https://cosmelab.github.io/entm201l-fall2025/setup/github-setup.html)

5. **Git** - For cloning the repository
   - See: [Software Setup Guide](https://cosmelab.github.io/entm201l-fall2025/setup/index.html)

**Platform-specific notes:**
- **Windows**: Use PowerShell, Git Bash, or **VS Code integrated terminal** (WSL2 recommended for Docker)
- **macOS**: Use Terminal or **VS Code integrated terminal** (recommended)
- **Linux**: Use your default terminal or **VS Code integrated terminal**

**🚨 If you haven't completed these setup steps, STOP and do them first!**

---

## 💻 Recommended Workflow: Use VS Code

**If you completed the ENTM201L setup, you already have VS Code configured!**

### Why Use VS Code?
- ✅ **Integrated terminal** - Run all commands without switching apps
- ✅ **Docker integration** - See containers, images, and logs
- ✅ **Dev Containers** - Open project inside Docker (advanced, optional)
- ✅ **Git integration** - Commit, push, pull with visual interface
- ✅ **File browser** - Navigate files easily
- ✅ **Syntax highlighting** - Read code and reports better

### Quick Start with VS Code

**1. Open this project in VS Code:**
```bash
# Navigate to the project directory first
cd ~/Desktop/dna-barcoding-analysis  # Adjust path to where you cloned it

# Open in VS Code
code .
```

**2. Open the integrated terminal:**
- **macOS/Linux**: Press `` Ctrl+` `` or View → Terminal
- **Windows**: Press `` Ctrl+` `` or View → Terminal

**3. Run all commands in the VS Code terminal:**
```bash
# Everything works the same as regular terminal
./tutorial.sh
./run-analysis.sh
```

**4. View HTML reports:**
- They open automatically in your default browser
- Or right-click HTML files → "Open with Live Server" (if you have the extension)

### Advanced: Dev Containers (Optional)

**For advanced users:** This project includes Dev Container configuration.

**What this means:**
- VS Code can run **inside** the Docker container
- No need to type `docker run` commands
- Terminal automatically runs in the container
- One-click setup

**To use:**
1. Install the "Dev Containers" extension (you should already have it from ENTM201L)
2. Open this project in VS Code
3. Click the green button in bottom-left corner
4. Select "Reopen in Container"

**For beginners:** Stick with the regular workflow for now. You can explore Dev Containers later!

### 🎨 Terminal Experience

**The container includes a beautiful terminal setup!**

When you run commands inside the Docker container, you get:
- **Zsh** - Modern shell with better autocompletion than bash
- **Oh-My-Zsh** - Popular framework with helpful plugins
- **Dracula theme** - Professional dark theme
- **Colorful output** - File listings with icons and colors
- **Git integration** - See git status in your prompt
- **Smart autosuggestions** - Suggests commands as you type

**To use the nice terminal:**

```bash
# Run the container interactively with zsh (instead of running scripts directly)
docker run --rm -it -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest zsh

# Now you're inside the container with the fancy terminal!
# Try these commands to see the colorful output:
ls         # Colorful file listings
ll         # Detailed view with git status
lt         # Tree view

# Run your analysis:
./tutorial.sh
./run-analysis.sh
```

**Why does this matter?**
- Makes working in the terminal more enjoyable
- Helps you see file types and git status at a glance
- Teaches you professional development tools
- Same setup used by many developers

**Note:** The fancy terminal only shows when you're inside the container interactively. When running scripts directly with `docker run ... ./tutorial.sh`, you won't see it (but everything still works!).

---

## 🎯 For Students: START HERE

**👉 See [start_here.md](start_here.md) for the complete beginner's guide**

### Quick Version (3 commands):

```bash
./tutorial.sh              # STEP 1: Learn with test data (15 min)
./run-analysis.sh          # STEP 2: Analyze YOUR data (5 min)
# STEP 3: Fill out assignment.md
```

**That's it!** Everything else below is reference documentation.

---

## 🔐 STEP 0: Log in to Docker (REQUIRED!)

**Before you can pull the container, you must log in to Docker Hub:**

```bash
docker login
```

You'll be prompted for:
- **Username**: Your Docker Hub username
- **Password**: Your Docker Hub password (typing is hidden for security)

**Success looks like:**
```
Login Succeeded
```

**Troubleshooting:**
- **"permission denied"**: Make sure Docker Desktop is running
- **"unauthorized"**: Double-check your username and password
- **Don't have an account?**: See Prerequisites above - create free account at [https://hub.docker.com/signup](https://hub.docker.com/signup)

**✅ Only need to do this ONCE per computer!** Docker will remember your login.

---

## 📚 STEP 1: Complete the Tutorial FIRST (REQUIRED!)

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
  results/my_analysis/01_qc/ \
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
  results/my_analysis/01_qc/passed_sequences.fasta \
  results/my_analysis/02_consensus/ \
  --pairs-only \
  --open
```

The `--pairs-only` flag keeps only samples where BOTH F and R passed QC.

### Step 3: Combine with Reference Sequences

Add your sequences to the database of known Southern California mosquitoes:

```bash
cat results/my_analysis/02_consensus/consensus_sequences.fasta \
    data/reference_sequences/socal_mosquitoes.fasta \
    > results/my_analysis/02_consensus/combined_with_references.fasta
```

This creates a file with your sequences + 52 reference sequences.

### Step 4: Alignment and Phylogenetic Tree

**Part A: Align sequences**

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/03_alignment/align_sequences.py \
  results/my_analysis/02_consensus/combined_with_references.fasta \
  results/my_analysis/03_alignment/
```

**Part B: Build tree** (takes ~2-3 minutes)

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/04_phylogeny/build_tree.py \
  results/my_analysis/03_alignment/aligned_sequences.fasta \
  results/my_analysis/04_phylogeny/
```

Look at `results/my_analysis/04_phylogeny/tree.png` to see where your samples cluster!

### Step 5: Species Identification (BLAST)

Compare your sequences to GenBank database:

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/05_identification/identify_species.py \
  results/my_analysis/02_consensus/consensus_sequences.fasta \
  results/my_analysis/05_blast/
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
├── start_here.md                # 👈 START HERE! Complete beginner's guide
├── tutorial.sh                  # Step 1: Interactive tutorial (test data)
├── run-analysis.sh              # Step 2: Analyze YOUR data
├── assignment.md                # Step 3: Fill this out
│
├── data/
│   ├── student_sequences/       # PUT YOUR .ab1 FILES HERE
│   ├── test_data/               # Tutorial uses this (8 .ab1 files)
│   └── reference_sequences/     # 52 known SoCal mosquito sequences
│
├── results/
│   ├── tutorial/                # Tutorial output (test data)
│   │   ├── 01_qc/
│   │   ├── 02_consensus/
│   │   ├── 03_alignment/
│   │   ├── 04_phylogeny/
│   │   └── 05_blast/
│   └── my_analysis/             # YOUR analysis output
│       ├── 01_qc/
│       ├── 02_consensus/
│       ├── 03_alignment/
│       ├── 04_phylogeny/
│       └── 05_blast/
│
├── modules/                     # Analysis scripts (used by container)
│   ├── 01_quality_control/
│   ├── 02_consensus/
│   ├── 03_alignment/
│   ├── 04_phylogeny/
│   └── 05_identification/
│
├── scripts/                     # Utility scripts
│   └── trim_references_to_barcode.py  # Trim GenBank references to ~700bp
│
├── docs/                        # Reference documentation
│   ├── pipeline_workflow.md
│   ├── iqtree_guide.md
│   └── reference_trimming.md    # Why we trim references to barcode region
│
└── intro_to_cli/                # Optional CLI tutorials (separate course)
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

50 mosquito COI sequences trimmed to the **barcode region** (~700bp):
- 19 species from 6 genera
- Genera: *Aedes*, *Anopheles*, *Culex*, *Deinocerites*, *Psorophora*, *Uranotaenia*
- All sequences validated and published
- **Trimmed from GenBank to match the 712bp AUCOS amplicon** (see `docs/reference_trimming.md`)

**Why trimmed?** GenBank contains COI sequences of varying lengths (640bp - 2,300bp). Using sequences of different lengths creates poor alignments with excessive gaps (~1,000 gaps!). We trimmed all references to the ~700bp barcode region that matches your student sequences.

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
4. **Reference sequence trimming**: `docs/reference_trimming.md`
5. **Check assignment**: `assignment.md`
6. **Ask your TA or instructor**

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

## 🌐 Documentation Website

A comprehensive documentation website is planned for future release, which will include:
- Interactive setup guides
- Step-by-step tutorials with screenshots
- Detailed module documentation
- Video walkthroughs
- FAQ and troubleshooting guides

For now, this README and the in-repo documentation provide everything you need to complete the analysis.

---

## 🔧 Development

### For Instructors and Developers

This repository includes:
- **Dev Container** configuration (`.devcontainer/`) for VS Code development
- **Modular CSS** design system (`tracking/styles/`) for HTML reports
- **Comprehensive tracking** (`tracking/`) for project management
- **Multi-architecture** Docker container (amd64 + arm64)

See `tracking/` directory for development documentation.

---

**Last Updated**: November 21, 2025
**Status**: ✅ Public Release - Student Ready
**Container**: `cosmelab/dna-barcoding-analysis:latest` (multi-architecture: amd64 + arm64)
**GitHub Classroom**: ✅ Compatible - use as template repository
**VS Code**: ✅ Dev Container included - optional advanced feature
