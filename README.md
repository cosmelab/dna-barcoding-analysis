# DNA Barcoding Analysis for COI Gene Sequencing

**Simple, automated analysis of Sanger sequences from the COI gene**

A containerized workflow for students with no coding experience to analyze their DNA barcoding data from chromatograms to phylogenetic trees and species identification.

![Status](https://img.shields.io/badge/Status-In_Development-orange?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-blue?style=flat-square)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED?style=flat-square&logo=docker&logoColor=white)

---

## What This Does

**One Command. Complete Analysis.**

Drop your `.ab1` chromatogram files into a folder → Run one command → Get a beautiful HTML report with:

- ✅ Quality control metrics
- ✅ Sequence alignment visualization
- ✅ Phylogenetic tree
- ✅ Species identification (BLAST)

**No coding required.** Everything runs in a Docker container with all tools pre-installed.

---

## Quick Start

### For Students (GitHub Classroom)

1. **Accept the assignment** (link provided by instructor)
2. **Clone your repo**:
   ```bash
   git clone https://github.com/entm201l-fall2025/dna-barcoding-YOUR-USERNAME.git
   cd dna-barcoding-YOUR-USERNAME
   ```

3. **Add your chromatograms**:
   ```bash
   # Copy your .ab1 files to the data folder
   cp ~/Desktop/*.ab1 data/my_sequences/
   ```

4. **Run the analysis**:
   ```bash
   docker compose up
   # OR
   docker run -v $(pwd):/workspace ghcr.io/cosmelab/dna-barcoding-analysis analyze-sequences
   ```

5. **View results**:
   ```bash
   open results/index.html
   # Your browser will open with interactive dashboard
   ```

---

## Features

### Automated Quality Control
- Parse `.ab1` chromatogram files
- Calculate quality metrics (Phred scores)
- Identify low-quality regions
- Generate interactive QC plots

### Alignment Visualization
- Multiple sequence alignment with MAFFT
- Visual representation of conserved vs variable positions
- Interactive HTML plots (zoom, pan, explore)
- Based on Hoque et al. 2022 methodology

### Phylogenetic Analysis
- IQ-TREE2 with automatic model selection
- Bootstrap support values (1000 replicates)
- Publication-quality tree figures
- R visualization with ggtree

### Species Identification
- Automatic BLAST against NCBI database
- Top 5 matches with % identity
- Confidence assessment
- Reference comparison with mosquito COI sequences

---

## What's Inside

### Container (~700MB)

Built on `mambaorg/micromamba:1.5.0` with:

**Analysis Tools**:
- BioPython (chromatogram parsing)
- MAFFT (sequence alignment)
- IQ-TREE2 (phylogenetic inference)
- BLAST+ (species identification)

**Visualization**:
- R with ape, ggtree, tidyverse
- Python with plotly, matplotlib
- Interactive HTML dashboard

**Beautiful Terminal**:
- zsh with oh-my-zsh
- Dracula theme
- colorls, lsd, starship prompt
- Auto-completion and syntax highlighting

### Reference Data

**85 mosquito COI sequences** from Hoque et al. 2022:
- 19 species from 6 genera
- *Aedes*, *Anopheles*, *Culex*, *Deinocerites*, *Psorophora*, *Uranotaenia*
- Improved AUCOS primers (67.5% success vs 16.7% Folmer)
- All sequences validated and published

**Citation**: Hoque MM, Valentine MJ, Kelly PJ, et al. Modification of the Folmer primers for the cytochrome c oxidase gene facilitates identification of mosquitoes. Parasites Vectors. 2022;15:437. doi:[10.1186/s13071-022-05494-2](https://doi.org/10.1186/s13071-022-05494-2)

---

## Repository Structure

```
dna-barcoding-analysis/
├── data/
│   ├── my_sequences/          # Put your .ab1 files here
│   └── reference/             # Mosquito COI sequences (Hoque et al)
├── results/                   # Auto-generated outputs
│   ├── qc_report.html        # Quality control metrics
│   ├── alignment.html        # Sequence alignment visualization
│   ├── tree.pdf              # Phylogenetic tree
│   ├── species_id.html       # BLAST results
│   └── index.html            # Main dashboard (open this!)
├── scripts/
│   ├── qc-chromatograms      # Quality control script
│   ├── align-sequences       # MAFFT wrapper
│   ├── build-tree            # IQTree2 wrapper
│   ├── identify-species      # BLAST automation
│   └── analyze-sequences     # Master pipeline (run this!)
├── container/
│   ├── Dockerfile            # Container definition
│   └── README.md             # Container documentation
├── docker-compose.yml        # Easy container startup
└── README.md                 # This file
```

---

## Installation

### Option 1: Docker Desktop (Recommended)

1. **Install Docker Desktop**: https://www.docker.com/products/docker-desktop/
2. **Pull the container**:
   ```bash
   docker pull ghcr.io/cosmelab/dna-barcoding-analysis:latest
   ```
3. **Ready to go!**

### Option 2: Build Locally

```bash
git clone https://github.com/cosmelab/dna-barcoding-analysis.git
cd dna-barcoding-analysis
docker build -t dna-barcoding container/
```

---

## Usage

### Full Pipeline (One Command)

```bash
# Inside the container or via docker run
analyze-sequences

# This runs:
# 1. QC your chromatograms
# 2. Trim low-quality regions
# 3. Align with reference sequences
# 4. Build phylogenetic tree
# 5. BLAST for species ID
# 6. Generate HTML dashboard
```

### Individual Steps

```bash
# Quality control only
qc-chromatograms data/my_sequences/

# Alignment only
align-sequences data/my_sequences/ data/reference/

# Tree only
build-tree results/alignment.fasta

# Species ID only
identify-species results/cleaned/
```

---

## Output Explained

### index.html Dashboard

Open `results/index.html` in your browser to see:

**Tab 1: Quality Control**
- Per-base quality scores
- Read length distribution
- Pass/fail summary
- Trimming recommendations

**Tab 2: Alignment View**
- Conserved positions (LARGE letters)
- Variable positions (small letters)
- Color-coded nucleotides
- Interactive zoom and pan

**Tab 3: Phylogenetic Tree**
- Maximum likelihood tree
- Bootstrap support values
- Your samples highlighted
- Reference species labeled

**Tab 4: Species Identification**
- BLAST top hits table
- % Identity scores
- GenBank accessions
- Confidence assessment

---

## GitHub Classroom Workflow

### For Instructors

1. **Create assignment** from this template repo
2. **Students accept link** → get personal copy
3. **Students analyze data** → commit results
4. **Review student repos** for grading

### For Students

```bash
# 1. Clone your assignment repo
git clone <your-classroom-repo-url>

# 2. Add your data
cp ~/Desktop/*.ab1 data/my_sequences/

# 3. Run analysis
docker compose up

# 4. Commit results
git add results/
git commit -m "Complete DNA barcoding analysis"
git push origin main
```

---

## Troubleshooting

### Container Issues

**Problem**: Docker won't start
- **Solution**: Restart Docker Desktop, check system resources

**Problem**: Permission denied
- **Solution**: On Linux, add user to docker group: `sudo usermod -aG docker $USER`

**Problem**: Container is slow
- **Solution**: Allocate more RAM (4GB+) in Docker Desktop settings

### Analysis Issues

**Problem**: No sequences pass QC
- **Solution**: Check chromatogram quality, may need re-sequencing

**Problem**: BLAST returns no hits
- **Solution**: Sequence may be contamination, check chromatogram

**Problem**: Tree has low bootstrap values
- **Solution**: Need more sequences or higher quality data

---

## Citation

If you use this workflow in your research, please cite:

```
Cosme, L. (2025). DNA Barcoding Analysis: Automated COI Gene Analysis Workflow.
GitHub: https://github.com/cosmelab/dna-barcoding-analysis
```

And the reference dataset:

```
Hoque MM, Valentine MJ, Kelly PJ, Barua S, Murillo DFB, Wang C. Modification
of the Folmer primers for the cytochrome c oxidase gene facilitates
identification of mosquitoes. Parasites Vectors. 2022;15(1):437.
doi:10.1186/s13071-022-05494-2
```

---

## Contributing

We welcome contributions! Areas where you can help:

- 🐛 Report bugs or issues
- 📚 Improve documentation
- 🧪 Add example datasets
- 🎨 Enhance visualizations
- 🔧 Optimize scripts

**See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines**

---

## License

**Code**: MIT License
**Educational Materials**: CC BY 4.0
**Reference Data**: See individual dataset citations

---

## Acknowledgments

- **UC Riverside Department of Entomology**
- **ENTM201L Students** (beta testers)
- **Hoque et al. 2022** for mosquito COI reference sequences
- **Tool developers**: MAFFT, IQ-TREE, BioPython, ggtree

---

## Support

**Questions?** Open an issue: https://github.com/cosmelab/dna-barcoding-analysis/issues

**Instructor**: Luciano Cosme
**Course**: ENTM201L - Molecular Biology Laboratory
**Institution**: UC Riverside, Department of Entomology

---

**Status**: 🚧 In Development | Target: Week 8, Fall 2025

**Last Updated**: November 18, 2025
