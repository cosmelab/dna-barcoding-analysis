# DNA Barcoding Analysis Container

A minimal, beautiful Docker container for DNA barcoding analysis with essential bioinformatics tools and a gorgeous terminal experience.

## Features

### Bioinformatics Tools
- **MAFFT**: Multiple sequence alignment
- **IQ-TREE2**: Maximum likelihood phylogenetic inference
- **BLAST**: Species identification via sequence similarity
- **BioPython**: Sequence manipulation and analysis
- **R packages**: `ape`, `ggtree` for phylogenetic visualization

### Python Packages
- `biopython`: Biological sequence analysis
- `pandas`: Data manipulation
- `matplotlib`: Static plotting
- `plotly`: Interactive HTML visualizations

### Beautiful Terminal
- **zsh** with **oh-my-zsh**
- **Dracula theme** for consistent colors
- **colorls**: Colorful file listings with icons
- **lsd**: Modern ls replacement with tree view
- **starship**: Fast, minimal prompt
- Auto-completion plugins

## Quick Start

### Build the container
```bash
docker-compose build
```

### Run the container
```bash
docker-compose up -d
docker exec -it dna-barcoding zsh
```

Or use Docker directly:
```bash
docker build -t dna-barcoding:latest container/
docker run -it -v $(pwd):/workspace dna-barcoding:latest
```

## Available Commands

Once inside the container, use these commands:

- `analyze-sequences` - Run the complete DNA barcoding pipeline
- `qc-sequences` - Quality control and filtering only
- `align-sequences` - Multiple sequence alignment only
- `build-tree` - Phylogenetic tree construction only
- `blast-identify` - BLAST-based species identification only

## Directory Structure

```
/workspace/
├── data/
│   ├── my_sequences/    # Your input FASTA files
│   └── reference/       # Reference sequences for comparison
├── results/             # Analysis outputs
└── scripts/             # Analysis scripts
```

## Terminal Features

### File Listing
- `ls` - Modern listing with lsd
- `ll` - Detailed listing with git status
- `la` - Show all files including hidden
- `lt` - Tree view

### Navigation
The Dracula theme provides:
- Color-coded file types
- Git status integration
- Fast, responsive prompt with starship
- Intelligent auto-completion

## Container Size

Target: <800MB (minimal installation)

## What's Included

### Core Tools (from conda-forge + bioconda)
```
mafft        - Multiple sequence alignment
iqtree       - Phylogenetic tree inference
blast        - Sequence similarity search
r-base       - R statistical computing
r-ape        - Phylogenetic analysis
ggtree       - Tree visualization
```

### Python Environment
```
Python 3.10
biopython
pandas
matplotlib
plotly
```

## What's NOT Included

This is a minimal container. We removed:
- Jupyter notebooks (not needed for command-line analysis)
- RNA-seq specific tools (STAR, salmon, DESeq2, etc.)
- Web development tools (Jekyll, Node.js, etc.)
- Heavy R packages (tidyverse, ggplot2, etc.)
- Unnecessary system libraries

## For Students

This container is designed for learning DNA barcoding:
1. Place your sequences in `/workspace/data/my_sequences/`
2. Run `analyze-sequences` to start
3. Find results in `/workspace/results/`
4. Enjoy the beautiful terminal while you work!

## Building from Scratch

```bash
cd container/
docker build -t dna-barcoding:latest .
```

Check the size:
```bash
docker images dna-barcoding:latest
```

## Troubleshooting

### Container won't start
```bash
docker-compose down
docker-compose up -d
```

### Permission issues
Make sure your user has access to Docker:
```bash
sudo usermod -aG docker $USER
```

### Need to rebuild
```bash
docker-compose build --no-cache
```

## License

This container is provided for educational purposes.

## Credits

- Based on the RNA-seq analysis container template
- Beautiful terminal setup inspired by modern development environments
- Bioinformatics tools from conda-forge and bioconda
