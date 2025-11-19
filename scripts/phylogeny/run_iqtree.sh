#!/bin/bash
################################################################################
# IQ-TREE Phylogenetic Analysis Wrapper
#
# PURPOSE:
#   User-friendly wrapper for IQ-TREE that handles common phylogenetic analysis
#   tasks. Provides default settings optimized for DNA barcoding (COI gene).
#
# PHYLOGENETIC CONCEPTS:
#   - Model Selection (MFP): Automatically selects the best substitution model
#   - Bootstrap Resampling: Tests confidence in tree topology (1000 replicates)
#   - Automatic Threading: Uses all available CPU cores
#   - Output: Best Maximum Likelihood tree with bootstrap support values
#
# PREREQUISITES:
#   - IQ-TREE 2.0 or newer (installed in PATH)
#   - Aligned FASTA file (from MAFFT or similar)
#
# USAGE:
#   ./run_iqtree.sh <alignment.fasta> [options]
#
# EXAMPLES:
#   # Quick analysis with default settings
#   ./run_iqtree.sh aligned_sequences.fasta
#
#   # Specify model explicitly
#   ./run_iqtree.sh aligned_sequences.fasta -m GTR+G+I
#
#   # More bootstrap replicates for publication
#   ./run_iqtree.sh aligned_sequences.fasta -bb 2000
#
#   # Use only 4 threads (for resource-limited systems)
#   ./run_iqtree.sh aligned_sequences.fasta -nt 4
#
# OUTPUT FILES:
#   *.treefile   - Best ML tree (Newick format) - USE THIS FOR DOWNSTREAM ANALYSIS
#   *.iqtree     - Full analysis report with statistics
#   *.log        - Run log (useful for debugging)
#   *.ckp.gz     - Checkpoint (for resuming interrupted runs)
#
# AUTHOR: Educational Materials for DNA Barcoding
# VERSION: 1.0
################################################################################

# Exit on any error
set -e

# Color codes for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

################################################################################
# FUNCTION: Print help message
################################################################################
print_help() {
    cat << EOF
${BLUE}IQ-TREE Wrapper for DNA Barcoding Analysis${NC}

${YELLOW}USAGE:${NC}
    ./run_iqtree.sh <alignment.fasta> [options]

${YELLOW}REQUIRED:${NC}
    <alignment.fasta>      Input multiple sequence alignment (FASTA format)

${YELLOW}OPTIONAL IQ-TREE PARAMETERS:${NC}
    -m <model>            Substitution model (default: MFP = auto-select)
                          Common models: GTR+G, GTR+G+I, TN93+G, HKY+G

    -bb <value>           Bootstrap replicates (default: 1000)
                          Higher = more robust but slower
                          Publication-quality: 2000-5000

    -nt <value>           Number of threads (default: AUTO)
                          Use 1 for single-core machines
                          Use number of available cores for speed

${YELLOW}COMMON EXAMPLES:${NC}
    # Basic analysis (recommended for learning)
    ./run_iqtree.sh coi_alignment.fasta

    # DNA barcoding with more bootstrap replicates
    ./run_iqtree.sh coi_alignment.fasta -bb 2000

    # Faster analysis (fewer bootstraps, fewer threads)
    ./run_iqtree.sh coi_alignment.fasta -bb 500 -nt 2

    # Specific model (if you already know which fits best)
    ./run_iqtree.sh coi_alignment.fasta -m GTR+G+I -bb 1000

${YELLOW}WHAT THE OUTPUT MEANS:${NC}
    The .treefile contains bootstrap values at each node.
    - Values > 95% = Strong support (species clearly different)
    - Values 70-95% = Moderate support (likely different)
    - Values < 70% = Weak support (uncertain relationship)

${YELLOW}OUTPUT FILES TO KEEP:${NC}
    → *.treefile - This is your phylogenetic tree! Use for visualization or downstream
    → *.iqtree - Read this for model selection details and tree statistics

${YELLOW}NEXT STEPS:${NC}
    1. Check the .iqtree file for model selection results
    2. Visualize with: python visualize_tree.py <treefile>
    3. Extract species relationships for identification

EOF
}

################################################################################
# FUNCTION: Check if IQ-TREE is installed
################################################################################
check_iqtree() {
    if ! command -v iqtree &> /dev/null; then
        echo -e "${RED}ERROR: IQ-TREE not found in PATH${NC}"
        echo "Install with: conda install -c bioconda iqtree"
        exit 1
    fi

    # Show version for transparency
    IQTREE_VERSION=$(iqtree --version 2>&1 | head -1)
    echo -e "${BLUE}Using: $IQTREE_VERSION${NC}"
}

################################################################################
# FUNCTION: Validate input file
################################################################################
validate_input() {
    local alignment_file=$1

    # Check file exists
    if [ ! -f "$alignment_file" ]; then
        echo -e "${RED}ERROR: File not found: $alignment_file${NC}"
        exit 1
    fi

    # Check file is not empty
    if [ ! -s "$alignment_file" ]; then
        echo -e "${RED}ERROR: File is empty: $alignment_file${NC}"
        exit 1
    fi

    # Basic FASTA validation (has > symbols)
    if ! grep -q "^>" "$alignment_file"; then
        echo -e "${RED}ERROR: File does not appear to be FASTA format${NC}"
        echo "FASTA files should have lines starting with '>' followed by sequences"
        exit 1
    fi

    # Count sequences
    SEQ_COUNT=$(grep -c "^>" "$alignment_file")
    echo -e "${BLUE}Input validation:${NC}"
    echo "  File: $alignment_file"
    echo "  Sequences: $SEQ_COUNT"

    # Warn if only 2-3 sequences
    if [ $SEQ_COUNT -lt 4 ]; then
        echo -e "${YELLOW}⚠ Warning: Only $SEQ_COUNT sequences detected${NC}"
        echo "  Phylogenetic trees need at least 4 sequences for meaningful topology"
    fi
}

################################################################################
# FUNCTION: Print IQ-TREE parameters being used
################################################################################
print_parameters() {
    local model=$1
    local bootstrap=$2
    local threads=$3
    local output_prefix=$4

    echo -e "${BLUE}IQ-TREE Parameters:${NC}"
    echo "  Model selection: $model"
    echo "  Bootstrap replicates: $bootstrap"
    echo "  Threads: $threads"
    echo "  Output prefix: $output_prefix"
    echo ""
    echo -e "${YELLOW}Understanding these parameters:${NC}"
    echo "  • Model: MFP tests multiple substitution models, picks the best"
    echo "  • Bootstrap: Resamples data 1000 times to estimate confidence"
    echo "  • Threads: Uses multiple CPU cores for faster computation"
    echo ""
}

################################################################################
# FUNCTION: Run IQ-TREE with educational output
################################################################################
run_analysis() {
    local alignment=$1
    local model=$2
    local bootstrap=$3
    local threads=$4
    local prefix=$5

    echo -e "${GREEN}Starting IQ-TREE analysis...${NC}"
    echo "This may take several minutes depending on sequence length and count."
    echo ""

    # Run IQ-TREE with specified parameters
    # Documentation of each flag:
    #   -s         Input alignment file
    #   -m         Substitution model (MFP = ModelFinder Plus)
    #   -bb        Ultrafast bootstrap replicates
    #   -nt        Number of threads (AUTO uses all available)
    #   -pre       Output prefix

    iqtree \
        -s "$alignment" \
        -m "$model" \
        -bb "$bootstrap" \
        -nt "$threads" \
        -pre "$prefix"

    echo ""
    echo -e "${GREEN}Analysis complete!${NC}"
}

################################################################################
# FUNCTION: Print summary of results
################################################################################
print_summary() {
    local prefix=$1
    local treefile="${prefix}.treefile"
    local iqtreefile="${prefix}.iqtree"

    if [ ! -f "$treefile" ]; then
        echo -e "${RED}ERROR: Tree file not created${NC}"
        exit 1
    fi

    echo -e "${GREEN}Results Summary:${NC}"
    echo ""
    echo "OUTPUT FILES CREATED:"
    echo "  ${BLUE}${prefix}.treefile${NC}  ← Your phylogenetic tree (most important!)"
    echo "  ${BLUE}${prefix}.iqtree${NC}     ← Detailed analysis report"
    echo "  ${BLUE}${prefix}.log${NC}        ← Run log file"
    echo ""

    # Extract some useful info from iqtree file
    if [ -f "$iqtreefile" ]; then
        echo "KEY RESULTS FROM ANALYSIS:"

        # Try to extract model selection result
        if grep -q "Best-fit model" "$iqtreefile"; then
            echo ""
            echo "MODEL SELECTION:"
            grep "Best-fit model" "$iqtreefile" | head -1 | sed 's/^/  /'
        fi

        # Try to extract log-likelihood
        if grep -q "Log-likelihood" "$iqtreefile"; then
            echo ""
            echo "TREE QUALITY:"
            grep "Log-likelihood" "$iqtreefile" | head -1 | sed 's/^/  /'
        fi
    fi

    echo ""
    echo -e "${BLUE}NEXT STEPS:${NC}"
    echo ""
    echo "1. Examine the tree:"
    echo "   ${YELLOW}cat ${prefix}.iqtree${NC}  (read full analysis)"
    echo ""
    echo "2. Visualize the tree:"
    echo "   ${YELLOW}python visualize_tree.py ${prefix}.treefile${NC}"
    echo ""
    echo "3. Extract relationships for species identification:"
    echo "   ${YELLOW}python root_tree.py ${prefix}.treefile${NC}"
    echo ""
    echo "4. Calculate genetic distances between sequences:"
    echo "   ${YELLOW}python calculate_distances.py ${prefix}.treefile${NC}"
    echo ""
}

################################################################################
# MAIN SCRIPT
################################################################################

# Show help if no arguments
if [ $# -eq 0 ]; then
    print_help
    exit 0
fi

# Parse help flag
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    print_help
    exit 0
fi

# Extract alignment file (first positional argument)
ALIGNMENT_FILE="$1"
shift  # Remove first argument, leaving only optional parameters

# Check IQ-TREE is available
check_iqtree

# Validate input alignment file
validate_input "$ALIGNMENT_FILE"

# Set defaults
MODEL="MFP"
BOOTSTRAP="1000"
THREADS="AUTO"
OUTPUT_PREFIX="${ALIGNMENT_FILE%.*}_iqtree"  # Remove extension, add suffix

# Parse optional arguments
# This allows IQ-TREE parameters to pass through
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m)
            MODEL="$2"
            shift 2
            ;;
        -bb)
            BOOTSTRAP="$2"
            shift 2
            ;;
        -nt)
            THREADS="$2"
            shift 2
            ;;
        *)
            # Unknown parameter - pass to IQ-TREE
            shift
            ;;
    esac
done

# Print analysis parameters
echo ""
print_parameters "$MODEL" "$BOOTSTRAP" "$THREADS" "$OUTPUT_PREFIX"

# Run the analysis
run_analysis "$ALIGNMENT_FILE" "$MODEL" "$BOOTSTRAP" "$THREADS" "$OUTPUT_PREFIX"

# Print results summary
print_summary "$OUTPUT_PREFIX"

echo -e "${BLUE}For more information:${NC}"
echo "  • Read about IQ-TREE: http://www.iqtree.org"
echo "  • Understand bootstrap values: http://www.iqtree.org/doc/Tutorial"
echo "  • Join discussions: https://github.com/iqtree/iqtree2/discussions"
echo ""
