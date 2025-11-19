#!/bin/bash

################################################################################
# MAFFT Alignment Wrapper Script
#
# Purpose: Wrapper script for running MAFFT sequence alignments with different
#          strategies and options. Demonstrates basic to advanced MAFFT usage.
#
# Requirements:
#   - MAFFT installed and available in PATH
#   - Input FASTA file(s)
#
# Usage:
#   ./run_mafft.sh -i input.fasta -o output.fasta [-m strategy] [-t threads]
#
# Examples:
#   # Basic FFT-NS-2 (default)
#   ./run_mafft.sh -i sequences.fasta -o aligned.fasta
#
#   # High accuracy for small alignments
#   ./run_mafft.sh -i sequences.fasta -o aligned.fasta -m accurate
#
#   # Fast mode for large alignments
#   ./run_mafft.sh -i sequences.fasta -o aligned.fasta -m fast
#
#   # With parallel processing (4 threads)
#   ./run_mafft.sh -i sequences.fasta -o aligned.fasta -m balanced -t 4
#
# MAFFT Strategies (from https://mafft.cbrc.jp/alignment/software/manual/manual.html):
#   - FFT-NS-1 (fastest, least accurate)
#   - FFT-NS-2 (default)
#   - FFT-NS-i (iterative refinement)
#   - L-INS-i (local alignment, high accuracy for <200 sequences)
#   - G-INS-i (global alignment, high accuracy for <200 sequences)
#   - E-INS-i (multiple domain alignment, for long sequences)
#
# References:
#   - MAFFT Manual: https://mafft.cbrc.jp/alignment/software/manual/manual.html
#   - Original Paper: Katoh et al. (2002) NAR, Katoh et al. (2018) MSB
#
################################################################################

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Default values
INPUT_FILE=""
OUTPUT_FILE=""
STRATEGY="default"
NUM_THREADS=1
QUIET=false
VERBOSE=false

################################################################################
# Function: Print usage information
################################################################################
usage() {
    cat << EOF
Usage: $0 -i INPUT_FILE -o OUTPUT_FILE [OPTIONS]

Required arguments:
  -i, --input FILE          Input FASTA file to align

Optional arguments:
  -o, --output FILE         Output alignment file (default: input.aligned.fasta)
  -m, --mode STRATEGY       Alignment strategy (default: default)
                            Options: fast, balanced, default, accurate
  -t, --threads NUM         Number of threads to use (default: 1)
  -q, --quiet               Suppress progress information
  -v, --verbose             Show detailed MAFFT output
  -h, --help                Show this help message

Examples:
  # Fast alignment with 4 threads
  $0 -i input.fasta -o output.fasta -m fast -t 4

  # Accurate alignment (slower, better for few sequences)
  $0 -i input.fasta -o output.fasta -m accurate

EOF
    exit 1
}

################################################################################
# Function: Print colored messages
################################################################################
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

################################################################################
# Function: Check if file exists and is readable
################################################################################
check_file() {
    if [[ ! -f "$1" ]]; then
        log_error "Input file not found: $1"
        exit 1
    fi
    if [[ ! -r "$1" ]]; then
        log_error "Input file is not readable: $1"
        exit 1
    fi
}

################################################################################
# Function: Check if MAFFT is installed
################################################################################
check_mafft_installed() {
    if ! command -v mafft &> /dev/null; then
        log_error "MAFFT is not installed or not in PATH"
        log_error "Install MAFFT from: https://mafft.cbrc.jp/alignment/software/"
        exit 1
    fi
    log_info "MAFFT found: $(command -v mafft)"
}

################################################################################
# Function: Build MAFFT command based on strategy
################################################################################
build_mafft_command() {
    local strategy="$1"
    local threads="$2"
    local input_file="$3"
    local output_file="$4"

    # Base command components
    local mafft_cmd="mafft"
    local mafft_args=""
    local thread_arg=""

    # Add threading if specified
    if [[ $threads -gt 1 ]]; then
        thread_arg="--thread $threads"
    fi

    # Select strategy and build command
    case "$strategy" in
        fast)
            # FFT-NS-1: Fastest, least accurate
            # Suitable for: Very large alignments (>5000 sequences)
            log_info "Using FAST mode (FFT-NS-1): Fastest, least accurate"
            mafft_args="--retree 1 --maxiterate 0"
            ;;
        balanced)
            # FFT-NS-2: Good balance of speed and accuracy
            # Suitable for: Medium-sized alignments (100-5000 sequences)
            log_info "Using BALANCED mode (FFT-NS-2): Default algorithm"
            mafft_args="--retree 2 --maxiterate 2"
            ;;
        default)
            # FFT-NS-2 with iterative refinement
            # Suitable for: Most general cases
            log_info "Using DEFAULT mode with iterative refinement"
            mafft_args="--retree 2 --maxiterate 100"
            ;;
        accurate)
            # L-INS-i or G-INS-i depending on sequence length
            # Suitable for: Small alignments (<200 sequences)
            log_info "Using ACCURATE mode (L-INS-i): High accuracy"
            mafft_args="--localpair --maxiterate 1000"
            ;;
        *)
            log_error "Unknown strategy: $strategy"
            exit 1
            ;;
    esac

    # Build final command
    echo "$mafft_cmd $mafft_args $thread_arg $input_file > $output_file"
}

################################################################################
# Function: Get sequence count from FASTA file
################################################################################
get_sequence_count() {
    # Count lines starting with '>' (FASTA headers)
    grep -c "^>" "$1" 2>/dev/null || echo "unknown"
}

################################################################################
# Parse command line arguments
################################################################################
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -m|--mode)
            STRATEGY="$2"
            shift 2
            ;;
        -t|--threads)
            NUM_THREADS="$2"
            shift 2
            ;;
        -q|--quiet)
            QUIET=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            ;;
    esac
done

################################################################################
# Validate input
################################################################################
if [[ -z "$INPUT_FILE" ]]; then
    log_error "Input file is required"
    usage
fi

check_file "$INPUT_FILE"
check_mafft_installed

# Set output file if not specified
if [[ -z "$OUTPUT_FILE" ]]; then
    OUTPUT_FILE="${INPUT_FILE%.fasta}.aligned.fasta"
    OUTPUT_FILE="${OUTPUT_FILE%.fa}.aligned.fa"
    if [[ "$OUTPUT_FILE" == "$INPUT_FILE" ]]; then
        OUTPUT_FILE="${INPUT_FILE}.aligned"
    fi
fi

################################################################################
# Display alignment information
################################################################################
SEQ_COUNT=$(get_sequence_count "$INPUT_FILE")
if [[ "$QUIET" != "true" ]]; then
    log_info "Input file: $INPUT_FILE"
    log_info "Output file: $OUTPUT_FILE"
    log_info "Strategy: $STRATEGY"
    log_info "Threads: $NUM_THREADS"
    log_info "Sequences to align: $SEQ_COUNT"
    echo ""
fi

################################################################################
# Build and execute MAFFT command
################################################################################
MAFFT_COMMAND=$(build_mafft_command "$STRATEGY" "$NUM_THREADS" "$INPUT_FILE" "$OUTPUT_FILE")

if [[ "$VERBOSE" == "true" ]]; then
    log_info "Executing: $MAFFT_COMMAND"
fi

# Check if output file already exists
if [[ -f "$OUTPUT_FILE" ]]; then
    log_warn "Output file already exists: $OUTPUT_FILE"
    read -p "Overwrite? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log_info "Aborted."
        exit 0
    fi
fi

# Execute alignment
if [[ "$QUIET" == "true" ]]; then
    eval "$MAFFT_COMMAND" 2>/dev/null
else
    eval "$MAFFT_COMMAND"
fi

################################################################################
# Verify output
################################################################################
if [[ -f "$OUTPUT_FILE" ]] && [[ -s "$OUTPUT_FILE" ]]; then
    ALIGNED_SEQ_COUNT=$(get_sequence_count "$OUTPUT_FILE")
    ALIGNMENT_LENGTH=$(grep -v "^>" "$OUTPUT_FILE" | head -1 | wc -c)
    alignment_length=$((ALIGNMENT_LENGTH - 1)) # subtract newline

    log_info "Alignment complete!"
    log_info "Output sequences: $ALIGNED_SEQ_COUNT"
    log_info "Alignment length: $alignment_length bp"
    log_info "Output saved to: $OUTPUT_FILE"
    exit 0
else
    log_error "Alignment failed or output file is empty"
    exit 1
fi
