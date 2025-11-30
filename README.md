<div align="center">

# DNA Barcoding Analysis Pipeline

### Automated COI Sequencing Analysis for Species Identification

![Status](https://img.shields.io/badge/status-production_ready-50fa7b?style=flat-square&labelColor=282a36)
![Course](https://img.shields.io/badge/course-ENTM201L-bd93f9?style=flat-square&labelColor=282a36)
![License](https://img.shields.io/badge/license-MIT-8be9fd?style=flat-square&labelColor=282a36)
![Platform](https://img.shields.io/badge/platform-macOS_|_Linux_|_Windows-ff79c6?style=flat-square&labelColor=282a36)
![Docker](https://img.shields.io/badge/docker-multi--arch-ffb86c?style=flat-square&labelColor=282a36)
![Codespaces](https://img.shields.io/badge/Codespaces-ready-50fa7b?style=flat-square&labelColor=282a36)

[Quick Start](start_here.md) · [Assignment](assignment.md) · [Documentation](docs/) · [Docker Hub](https://hub.docker.com/r/cosmelab/dna-barcoding-analysis)

</div>

---

## Overview

A complete automated workflow for analyzing Sanger sequencing chromatograms (.ab1 files) to identify mosquito species using COI DNA barcoding. Built for ENTM201L students at UC Riverside with zero coding experience required.

```
Institution       University of California, Riverside
Course            ENTM201L - Molecular Biology Laboratory
Target Users      Graduate students (no coding required)
Analysis Time     ~5 minutes (fully automated)
Instructor        Luciano Cosme, Department of Entomology
```

---

## Quick Start

### For Students

See **[start_here.md](start_here.md)** for the complete beginner's guide.

**Two ways to run:**

#### Option A: GitHub Codespaces (Recommended - No Installation!)

**Step-by-step to open Codespaces:**

1. Go to your repository page on GitHub
2. Look for the green **"<> Code"** button near the top-right of the page
3. Click it to open a dropdown menu
4. You'll see two tabs: **"Local"** and **"Codespaces"** — click **"Codespaces"**
5. Click the green **"Create codespace on main"** button
6. Wait 2-3 minutes while GitHub builds your environment (you'll see a loading screen)
7. When ready, you'll see VS Code in your browser with a terminal at the bottom

```
┌─────────────────────────────────────────────────────────────┐
│  [Your Repository Name]                                     │
│                                                             │
│  ┌──────────────────────────────────────┐                   │
│  │  <> Code ▼  │  (green button)        │   ← CLICK HERE    │
│  └──────────────────────────────────────┘                   │
│                                                             │
│  ┌──────────────────────────────────────┐                   │
│  │  Local    │  Codespaces             │   ← SELECT TAB    │
│  ├──────────────────────────────────────┤                   │
│  │                                      │                   │
│  │  ┌────────────────────────────────┐  │                   │
│  │  │ + Create codespace on main     │  │   ← CLICK THIS   │
│  │  └────────────────────────────────┘  │                   │
│  │                                      │                   │
│  └──────────────────────────────────────┘                   │
└─────────────────────────────────────────────────────────────┘
```

**Once Codespaces opens**, run these commands in the terminal:

> **Tip: Nice Terminal** — Type `zsh` for a colorful terminal with the Dracula theme!

> **Refresh Files** — After running scripts, click the 🔄 refresh icon in the Explorer panel to see new results.

> **Viewing HTML Reports** — Right-click any `.html` file → **"Download"** → Open in your browser. All reports work offline.

```bash
# STEP 1: Learn with test data (5 min)
./tutorial-cs.sh

# STEP 2: Analyze the class mosquito sequences (5 min)
./run-analysis-cs.sh

# STEP 3: Answer questions interactively (10 min)
python3 answer_assignment.py

# STEP 4: Submit to GitHub (auto-graded!)
git add submission/answers.json
git commit -m "Complete assignment"
git push origin main

# OPTIONAL: Generate your personal lab report
./student_report-cs.sh YOUR_CODE    # e.g., ./student_report-cs.sh HV

# OPTIONAL: Check your progress
python3 check_progress.py
```

#### Option B: Local Docker (Mac/Windows/Linux)

**IMPORTANT:** Run these commands **on your computer** (Mac/Windows/Linux), NOT inside Docker.

**Requirements:** Docker Desktop must be running!

> **Tip: Nice Terminal** — Want to explore inside the container? Run:
> ```bash
> docker run --rm -it -v $(pwd):/workspace -w /workspace cosmelab/dna-barcoding-analysis:latest zsh
> ```
> This gives you a colorful Dracula-themed terminal with the analysis tools.

```bash
# STEP 1: Learn with test data (5 min)
./tutorial.sh

# STEP 2: Analyze the class mosquito sequences (5 min)
./run-analysis.sh

# STEP 3: Answer questions interactively (10 min)
python3 answer_assignment.py

# STEP 4: Submit to GitHub (auto-graded!)
git add submission/answers.json
git commit -m "Complete assignment"
git push origin main

# OPTIONAL: Generate your personal lab report
./student_report.sh YOUR_CODE    # e.g., ./student_report.sh HV

# OPTIONAL: Check your progress
python3 check_progress.py
```

### Script Reference Table

| Script | Codespaces | Local Docker | Purpose |
|--------|------------|--------------|---------|
| `tutorial-cs.sh` | ✓ | | Learn pipeline with test data |
| `tutorial.sh` | | ✓ | Learn pipeline with test data |
| `run-analysis-cs.sh` | ✓ | | Analyze class sequences |
| `run-analysis.sh` | | ✓ | Analyze class sequences |
| `student_report-cs.sh` | ✓ | | Generate personal lab report |
| `student_report.sh` | | ✓ | Generate personal lab report |
| `answer_assignment.py` | ✓ | ✓ | Answer assignment questions |
| `check_progress.py` | ✓ | ✓ | Check your progress |

### For Instructors

```bash
# Clone as GitHub Classroom template
git clone https://github.com/cosmelab/dna-barcoding-analysis.git

# Students get their own repos:
# github.com/cosmelab/dna-barcoding-analysis-STUDENT-USERNAME
```

**GitHub Classroom compatible** — use as template repository

---

## The 6-Step Pipeline

<details>
<summary><strong>Step 1: Quality Control</strong></summary>

<br>

**What it does:**
- Analyzes .ab1 chromatogram files
- Checks quality scores (Phred Q30+ required)
- Validates sequence length (>500bp required)
- Filters out low-quality reads

**Output:**
- `qc_report.html` — Interactive quality control report
- `passed_sequences.fasta` — High-quality sequences only

**Why it matters:** Garbage in = garbage out. Bad sequences produce unreliable species IDs.

</details>

<details>
<summary><strong>Step 2: Consensus Sequences</strong></summary>

<br>

**What it does:**
- Pairs forward (F) and reverse (R) reads
- Reverse-complements the R read
- Creates consensus sequence from F+R alignment
- Filters for complete pairs only

**Output:**
- `consensus_sequences.fasta` — Final consensus sequences
- `consensus_report.html` — Alignment visualization

**Why it matters:** Combining F+R reads doubles coverage and accuracy.

</details>

<details>
<summary><strong>Step 3: Combine with References</strong></summary>

<br>

**What it does:**
- Adds your sequences to 52 reference mosquito COI sequences
- References include 19 species from 6 genera (*Aedes*, *Anopheles*, *Culex*, *Deinocerites*, *Psorophora*, *Uranotaenia*)
- All references trimmed to ~700bp barcode region

**Output:**
- `combined_with_references.fasta` — Your sequences + references

**Why it matters:** Can't build a tree without known species for context.

</details>

<details>
<summary><strong>Step 4: Phylogenetic Tree</strong></summary>

<br>

**What it does:**
- Aligns all sequences with MAFFT
- Builds maximum likelihood tree with IQ-TREE2
- Calculates 1000 ultrafast bootstrap replicates
- Generates 4 tree layouts (rectangular, circular, unrooted, radial)

**Output:**
- `tree.png`, `tree_circular.pdf`, etc. — Tree visualizations
- `phylogeny_report.html` — Interactive tree explorer

**Why it matters:** Shows evolutionary relationships. Your samples cluster with related species.

</details>

<details>
<summary><strong>Step 5: Species Identification</strong></summary>

<br>

**What it does:**
- BLASTs your sequences against NCBI GenBank
- Returns top 10 matches with % identity
- Interprets results (>98% = same species)

**Output:**
- `identification_report.html` — BLAST results table
- Top hits with accession numbers and % identity

**Why it matters:** Confirms species ID from tree with global database.

</details>

<details>
<summary><strong>Step 6: Lab Data Analysis (Interactive Visualizations)</strong></summary>

<br>

**What it does:**
- Analyzes class lab data (DNA extraction, PCR, sequencing)
- Creates interactive Plotly visualizations
- Generates personalized student reports
- Compares Team Spin vs Team Magnet performance

**Output:**
- `lab_report.html` — Interactive dashboard with all class results
- Individual student reports in `student_reports/`
- DNA yield comparisons, PCR success rates, sequencing QC metrics

**Why it matters:** Visualize and understand the entire lab workflow, from extraction to species ID.

</details>

---

## What's Included

### Bioinformatics Tools

- **BioPython** — Chromatogram parsing and sequence handling
- **MAFFT** — Multiple sequence alignment (industry standard)
- **IQ-TREE2** — Maximum likelihood phylogenetic inference
- **BLAST+** — Species identification via NCBI GenBank
- **toytree** — Beautiful tree visualizations with genus coloring

### Interactive HTML Reports

- Quality control dashboard with chromatogram viewer
- Consensus sequence comparisons (F vs R alignment)
- Alignment heatmaps (conservation visualization)
- Phylogenetic trees (4 layouts, genus-colored)
- BLAST results tables (sortable, interactive)

### Terminal Environment

- **Zsh** with oh-my-zsh framework
- **Dracula theme** — professional dark colors
- **Git integration** — see status in prompt
- **Aliases**: `ll` (detailed view), `lt` (tree view)

### Reference Dataset

- 52 Southern California mosquito COI sequences
- 19 species from 6 genera
- All trimmed to ~700bp barcode region
- Published sequences from Hoque et al. 2022

---

## Repository Structure

```
dna-barcoding-analysis/
├── start_here.md                 # Complete beginner's guide (START HERE!)
├── assignment.md                 # Student assignment questions
├── tutorial.sh                   # Learn with test data (Docker)
├── tutorial-cs.sh                # Learn with test data (Codespaces)
├── run-analysis.sh               # Analyze YOUR data (Docker)
├── run-analysis-cs.sh            # Analyze YOUR data (Codespaces)
│
├── data/
│   ├── student_sequences/        # PUT YOUR .ab1 FILES HERE
│   ├── test_data/                # 8 test chromatograms (for tutorial)
│   └── reference_sequences/      # 52 known mosquito sequences
│
├── results/
│   ├── tutorial/                 # Tutorial output (test data)
│   ├── my_analysis/              # YOUR analysis output
│   │   ├── 01_qc/                # Quality control results
│   │   ├── 02_consensus/         # Consensus sequences
│   │   ├── 03_alignment/         # MAFFT alignment
│   │   ├── 04_phylogeny/         # Trees (4 layouts)
│   │   └── 05_blast/             # Species identification
│   ├── lab_analysis/             # Lab data visualizations
│   └── student_reports/          # Individual student reports
│
├── modules/                      # Python analysis scripts
│   ├── 01_quality_control/
│   ├── 02_consensus/
│   ├── 03_alignment/
│   ├── 04_phylogeny/
│   ├── 05_identification/
│   └── 06_lab_data_analysis/     # Lab data visualizations
│
├── docs/                         # Documentation
│   ├── pipeline_workflow.md
│   ├── iqtree_guide.md
│   └── reference_trimming.md
│
├── .devcontainer/                # Codespaces configuration
│
└── intro_to_cli/                 # Optional CLI tutorials
```

---

## System Requirements

### Option A: GitHub Codespaces (Recommended - No Installation!)

Just a web browser and GitHub account. Everything runs in the cloud.

### Option B: Local Docker

```
Docker Desktop      Windows 10+, macOS 10.15+, or Linux
Git                 For cloning the repository
Docker Hub Account  Free account (no payment needed)
```

### Recommended (Optional)

- **VS Code** — Best experience with integrated terminal
- **4GB+ RAM** allocated to Docker for faster tree building
- **Internet connection** — For BLAST searches

### Platform Support

| Platform | Status | Architecture |
|----------|--------|--------------|
| **macOS** (Intel) | Fully Supported | Native amd64 |
| **macOS** (Apple Silicon) | Fully Supported | Native arm64 |
| **Windows 10/11** | Fully Supported | Requires WSL2 |
| **Linux** | Fully Supported | Native support |

Multi-architecture container automatically uses the correct version for your system.

---

## Learning Outcomes

Upon completing this workflow, students will be able to:

- **Assess DNA sequence quality** from chromatogram data
- **Interpret quality metrics** (Phred scores, base calling)
- **Understand consensus sequences** and why F+R reads matter
- **Read phylogenetic trees** and identify evolutionary relationships
- **Perform species identification** using BLAST and % identity
- **Use Docker containers** for reproducible bioinformatics
- **Navigate the command line** with confidence

---

## Adaptable for Other Projects

**Current implementation:** COI barcoding for mosquito identification

**Potential adaptations:** This pipeline could be modified for other Sanger sequencing projects by replacing reference sequences and .ab1 files. The workflow (QC → Consensus → Alignment → Tree → BLAST) works for any organism with GenBank data.

### Possible Barcode Regions

Examples of what this pipeline could be adapted for:
- **ITS** (fungi)
- **rbcL, matK** (plants)
- **16S rRNA** (bacteria)
- **COI** (other animals)

### How to Adapt

```bash
# 1. Replace reference sequences
python3 data/reference_sequences/trim_references_to_barcode.py \
  your_genbank_refs.fasta your_trimmed_refs.fasta --start 50 --end 750

# 2. Replace student sequences (put your .ab1 files here)
cp ~/my_chromatograms/*.ab1 data/student_sequences/

# 3. Run the pipeline (same workflow!)
./run-analysis.sh
```

BLAST automatically searches NCBI for any organism.

---

## Troubleshooting

<details>
<summary><strong>Docker Issues</strong></summary>

<br>

**"Cannot connect to Docker daemon"**
- Make sure Docker Desktop is running
- Check system tray/menu bar for Docker icon

**"Permission denied"**
- Run `docker login` with your Docker Hub credentials
- Windows: Make sure WSL2 integration is enabled in Docker settings

**Container is slow**
- Allocate more RAM to Docker (Settings → Resources)
- Recommended: 4GB+ for tree building

</details>

<details>
<summary><strong>Analysis Issues</strong></summary>

<br>

**No sequences pass QC**
- Check chromatogram quality — may need re-sequencing
- Look at the QC report HTML to see failure reasons

**BLAST returns no hits**
- Sequence may be contamination or very poor quality
- Check alignment — might be wrong reading frame

**Tree has low bootstrap values**
- Normal for closely related species
- Add more reference sequences for better resolution

</details>

<details>
<summary><strong>Windows/WSL Issues</strong></summary>

<br>

**Enable Virtualization in BIOS**
- Restart → Enter BIOS (F2, F10, Del, or Esc)
- Enable "Intel Virtualization Technology" or "AMD-V"

**Docker Commands Hang**
- Restart WSL: `wsl --shutdown` in PowerShell (as Admin)
- Reopen WSL terminal and try again

**Verify Docker Works**
```bash
docker run hello-world
```

</details>

---

## Getting Help

1. **Read [start_here.md](start_here.md)** — Complete beginner's guide
2. **Check [docs/pipeline_workflow.md](docs/pipeline_workflow.md)** — Visual workflow
3. **Read [docs/iqtree_guide.md](docs/iqtree_guide.md)** — Understanding trees
4. **Ask your instructor** — Office hours available

---

## Pipeline Statistics

| Metric | Value |
|--------|-------|
| **Tutorial Time** | ~3 minutes (all 6 steps) |
| **Analysis Time** | ~5 minutes (automated) |
| **Reference Sequences** | 52 mosquito COI sequences |
| **Species Covered** | 19 species from 6 genera |
| **Tree Layouts** | 4 visualizations (rectangular, circular, unrooted, radial) |
| **HTML Reports** | 6 interactive dashboards + student reports |
| **Container Size** | ~2.5GB (includes all tools) |

---

## Citation

If you use this pipeline in your research or teaching, please cite:

**Reference sequences:**
> Hoque MM, Valentine MJ, Kelly PJ, et al. Modification of the Folmer primers for the cytochrome c oxidase gene facilitates identification of mosquitoes. *Parasites Vectors*. 2022;15:437. [doi:10.1186/s13071-022-05494-2](https://doi.org/10.1186/s13071-022-05494-2)

---

## For Developers & Advanced Users

### GitHub CLI Authentication

Need to interact with GitHub Packages or use `gh` CLI? Set up authentication.

#### 1. Create GitHub Personal Access Token (Classic)

**Important:** Use a **classic token** (not fine-grained) for full API access.

Go to: https://github.com/settings/tokens/new

**Configure your token:**

- **Note:** Give it a descriptive name (e.g., `gh-cli-packages`)
- **Expiration:** Choose 90 days or custom (tokens expire for security)

**Select scopes** (check these boxes):

- **`repo`** (Full control of private repositories)
  - This automatically checks all sub-scopes under `repo`
- **`read:packages`** (Download packages from GitHub Package Registry)
- **`write:packages`** (Upload packages to GitHub Package Registry)
- **`delete:packages`** (Delete packages — optional, for cleanup)

Click **"Generate token"** at the bottom.

**IMPORTANT:** Copy the token immediately — you won't see it again.

Token format: Starts with `ghp_` followed by 36 characters
- Example: `ghp_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx`

#### 2. Store Token Securely

```bash
# Create secure token file (only you can read it)
touch ~/.github_token
chmod 600 ~/.github_token

# Add your token (replace YOUR_TOKEN with actual token)
echo "export GITHUB_TOKEN=YOUR_TOKEN" > ~/.github_token

# Load in shell profiles
echo 'source ~/.github_token 2>/dev/null' >> ~/.zshrc
echo 'source ~/.github_token 2>/dev/null' >> ~/.bashrc

# Reload shell
source ~/.zshrc
```

#### 3. Verify Setup

```bash
# Check if token is loaded
echo $GITHUB_TOKEN

# Test GitHub CLI
gh api /user --jq '.login'
# Should print your username
```

#### 4. Use with Docker

```bash
# Login to GitHub Container Registry
echo $GITHUB_TOKEN | docker login ghcr.io -u YOUR_USERNAME --password-stdin

# Pull from GitHub Packages
docker pull ghcr.io/cosmelab/dna-barcoding-analysis:latest
```

Full documentation: [GitHub CLI Setup Guide](docs/github_cli_setup.md)

### Container Development

Modify the Docker container:

```bash
# Edit container/Dockerfile
vim container/Dockerfile

# Build locally (test before pushing)
cd container
./build.sh

# Push to main branch → GitHub Actions auto-builds and publishes
git add container/Dockerfile
git commit -m "Update container"
git push origin main
```

Auto-publish to:
- Docker Hub: `docker.io/cosmelab/dna-barcoding-analysis:latest`
- GitHub Packages: `ghcr.io/cosmelab/dna-barcoding-analysis:latest`

### Contributing

Want to improve the pipeline or add features?

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/amazing-feature`
3. Make your changes
4. Test with both `./tutorial.sh` and `./run-analysis.sh`
5. Commit: `git commit -m "Add amazing feature"`
6. Push: `git push origin feature/amazing-feature`
7. Open a Pull Request

---

## License

![MIT License](https://img.shields.io/badge/code-MIT-50fa7b?style=flat-square&labelColor=282a36)
![CC BY 4.0](https://img.shields.io/badge/docs-CC_BY_4.0-8be9fd?style=flat-square&labelColor=282a36)

- **Code:** MIT License
- **Educational Materials:** Creative Commons Attribution 4.0 (CC BY 4.0)
- **Reference Data:** See individual citations

You are free to:
- Use for teaching and research
- Modify and adapt for your needs
- Share with attribution

---

## Acknowledgments

- UC Riverside Department of Entomology
- ENTM201L Students (Fall 2025)
- Hoque et al. 2022 for Southern California mosquito COI sequences
- Open-source developers: BioPython, MAFFT, IQ-TREE, BLAST+, toytree teams
- GitHub for hosting and GitHub Classroom infrastructure
- Docker Hub for container distribution

---

<div align="center">

## Quick Links

[Get Started](start_here.md) · [Assignment](assignment.md) · [Docker Hub](https://hub.docker.com/r/cosmelab/dna-barcoding-analysis) · [Course Website](https://cosmelab.github.io/entm201l-fall2025/)

---

**Last Updated**: November 29, 2025
**Status**: Production Ready — Student Tested
**Container**: `cosmelab/dna-barcoding-analysis:latest` (multi-arch: amd64 + arm64)
**GitHub Classroom**: Template Ready
**GitHub Codespaces**: Fully Supported

**Instructor**: Luciano Cosme | Department of Entomology | UC Riverside

![UCR Entomology](https://img.shields.io/badge/UCR-Entomology-FFB81C?style=flat-square&labelColor=282a36)

</div>
