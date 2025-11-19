#!/bin/bash

################################################################################
# batch_qc.sh - Batch quality control processing for DNA sequences
#
# Educational example: This bash script demonstrates how to automate
# quality control workflows across multiple sequence files.
#
# What this script does:
# 1. Processes all .ab1 files in a directory
# 2. Converts them to FASTQ format
# 3. Trims low-quality bases
# 4. Filters sequences by quality metrics
# 5. Generates a summary report
#
# Bash scripting concepts demonstrated:
# - Directory traversal with find and for loops
# - Command piping and output redirection
# - Error handling with conditional statements
# - Variable expansion and string manipulation
# - Function definitions for code reuse
# - Logging and reporting
#
# Usage:
#     # Run QC on AB1 files in current directory
#     bash batch_qc.sh
#
#     # Run QC on specific directory
#     bash batch_qc.sh /path/to/ab1_files
#
#     # Run QC with custom quality threshold
#     QUALITY_THRESHOLD=25 bash batch_qc.sh
#
# Requirements:
#     - Python 3
#     - BioPython library
#     - parse_ab1.py, trim_quality.py, filter_sequences.py in same directory
#
# Environment variables (optional):
#     QUALITY_THRESHOLD - Minimum quality score (default: 20)
#     MIN_LENGTH - Minimum sequence length in bp (default: 300)
#     MAX_N - Maximum N percentage allowed (default: 5)
################################################################################

set -e  # Exit on any error

# ============================================================================
# CONFIGURATION
# ============================================================================

# Directory containing this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Input and output directories
INPUT_DIR="${1:-.}"  # Use first argument or current directory
OUTPUT_DIR="${INPUT_DIR}/qc_results"

# Quality control parameters
QUALITY_THRESHOLD="${QUALITY_THRESHOLD:-20}"  # Minimum quality score
MIN_LENGTH="${MIN_LENGTH:-300}"                # Minimum sequence length (bp)
MAX_N="${MAX_N:-5}"                            # Maximum N percentage

# Log file for recording what happened
LOG_FILE="${OUTPUT_DIR}/qc_log.txt"

# ============================================================================
# FUNCTIONS
# ============================================================================

# Print a message with timestamp to both console and log file
log_message() {
    local message="$1"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] ${message}" | tee -a "${LOG_FILE}"
}

# Print a section header
print_header() {
    local text="$1"
    echo ""
    echo "=========================================="
    echo "${text}"
    echo "=========================================="
    echo "" | tee -a "${LOG_FILE}"
}

# Check if a Python script exists and is executable
check_script() {
    local script="$1"
    local script_path="${SCRIPT_DIR}/${script}"

    if [ ! -f "${script_path}" ]; then
        echo "Error: ${script} not found in ${SCRIPT_DIR}" >&2
        return 1
    fi

    return 0
}

# Convert a single AB1 file to FASTQ
convert_ab1_to_fastq() {
    local ab1_file="$1"
    local fastq_file="$2"

    log_message "Converting ${ab1_file} to FASTQ..."

    # Use parse_ab1.py to extract sequence and save as FASTA
    # Note: This is a simplified version - a full implementation would
    # convert to FASTQ with quality scores
    python3 "${SCRIPT_DIR}/parse_ab1.py" "${ab1_file}" \
        --save-fasta "${fastq_file%.fastq}.fasta" \
        --save-qual "${fastq_file%.fastq}.qual" \
        >> "${LOG_FILE}" 2>&1

    if [ $? -eq 0 ]; then
        log_message "  Successfully converted to ${fastq_file}"
        return 0
    else
        log_message "  ERROR: Failed to convert ${ab1_file}"
        return 1
    fi
}

# Trim a single sequence file
trim_sequence_file() {
    local input_file="$1"
    local output_file="$2"

    log_message "Trimming ${input_file} (quality >= ${QUALITY_THRESHOLD})..."

    python3 "${SCRIPT_DIR}/trim_quality.py" \
        "${input_file}" \
        --quality "${QUALITY_THRESHOLD}" \
        --output "${output_file}" \
        >> "${LOG_FILE}" 2>&1

    if [ $? -eq 0 ]; then
        log_message "  Successfully trimmed to ${output_file}"
        return 0
    else
        log_message "  ERROR: Failed to trim ${input_file}"
        return 1
    fi
}

# Filter sequences by quality criteria
filter_sequence_file() {
    local input_file="$1"
    local output_file="$2"

    log_message "Filtering ${input_file} (min_length=${MIN_LENGTH}, max_n=${MAX_N}%)..."

    python3 "${SCRIPT_DIR}/filter_sequences.py" \
        "${input_file}" \
        --min-length "${MIN_LENGTH}" \
        --min-quality "${QUALITY_THRESHOLD}" \
        --max-n "${MAX_N}" \
        --output "${output_file}" \
        >> "${LOG_FILE}" 2>&1

    if [ $? -eq 0 ]; then
        log_message "  Successfully filtered to ${output_file}"
        return 0
    else
        log_message "  ERROR: Failed to filter ${input_file}"
        return 1
    fi
}

# Process a single sequence file through the QC pipeline
process_sequence() {
    local sequence_file="$1"
    local base_name=$(basename "${sequence_file}" | sed 's/\.[^.]*$//')

    log_message "Processing: ${base_name}"

    # Define output files
    local trimmed_file="${OUTPUT_DIR}/${base_name}_trimmed.fastq"
    local filtered_file="${OUTPUT_DIR}/${base_name}_qc.fastq"

    # Step 1: Trim low-quality bases
    trim_sequence_file "${sequence_file}" "${trimmed_file}" || return 1

    # Step 2: Filter sequences by quality criteria
    filter_sequence_file "${trimmed_file}" "${filtered_file}" || return 1

    log_message "  COMPLETE: ${base_name}"
    return 0
}

# Generate a summary report of all processed files
generate_report() {
    print_header "QC SUMMARY REPORT"

    log_message "Report generated: $(date)"
    log_message "Input directory: ${INPUT_DIR}"
    log_message "Output directory: ${OUTPUT_DIR}"
    log_message ""
    log_message "QC Parameters:"
    log_message "  Minimum quality: ${QUALITY_THRESHOLD}"
    log_message "  Minimum length: ${MIN_LENGTH} bp"
    log_message "  Maximum N percentage: ${MAX_N}%"
    log_message ""

    # Count files
    local ab1_count=$(find "${INPUT_DIR}" -maxdepth 1 -name "*.ab1" 2>/dev/null | wc -l)
    local qc_count=$(find "${OUTPUT_DIR}" -name "*_qc.fastq" 2>/dev/null | wc -l)

    log_message "Results:"
    log_message "  AB1 files found: ${ab1_count}"
    log_message "  QC files generated: ${qc_count}"
    log_message ""
    log_message "Output files:"
    find "${OUTPUT_DIR}" -name "*_qc.fastq" -exec ls -lh {} \; | tee -a "${LOG_FILE}"
}

# ============================================================================
# MAIN WORKFLOW
# ============================================================================

main() {
    print_header "DNA SEQUENCE QUALITY CONTROL PIPELINE"

    # Verify input directory exists
    if [ ! -d "${INPUT_DIR}" ]; then
        echo "Error: Input directory '${INPUT_DIR}' not found" >&2
        exit 1
    fi

    # Create output directory
    if [ ! -d "${OUTPUT_DIR}" ]; then
        mkdir -p "${OUTPUT_DIR}"
        log_message "Created output directory: ${OUTPUT_DIR}"
    fi

    # Initialize log file
    log_message "Starting QC pipeline"
    log_message "Input directory: ${INPUT_DIR}"
    log_message "Output directory: ${OUTPUT_DIR}"

    # Verify required Python scripts exist
    print_header "CHECKING DEPENDENCIES"
    check_script "parse_ab1.py" || exit 1
    check_script "trim_quality.py" || exit 1
    check_script "filter_sequences.py" || exit 1
    log_message "All required scripts found"

    # Find and process all AB1 files
    print_header "PROCESSING SEQUENCE FILES"

    local ab1_files=()
    mapfile -t ab1_files < <(find "${INPUT_DIR}" -maxdepth 1 -name "*.ab1" -type f)

    if [ ${#ab1_files[@]} -eq 0 ]; then
        log_message "Warning: No .ab1 files found in ${INPUT_DIR}"
        log_message "The script will exit without processing any files."
        exit 0
    fi

    local processed=0
    local failed=0

    # Process each AB1 file
    for ab1_file in "${ab1_files[@]}"; do
        if process_sequence "${ab1_file}"; then
            ((processed++))
        else
            ((failed++))
        fi
    done

    # Generate summary report
    print_header "PIPELINE COMPLETE"
    generate_report

    log_message "Processing complete: ${processed} succeeded, ${failed} failed"
    log_message "Log file: ${LOG_FILE}"

    # Exit with error status if any files failed
    if [ ${failed} -gt 0 ]; then
        exit 1
    fi

    exit 0
}

# ============================================================================
# SCRIPT ENTRY POINT
# ============================================================================

# Run the main function
main "$@"
