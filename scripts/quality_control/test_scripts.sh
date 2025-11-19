#!/bin/bash

################################################################################
# test_scripts.sh - Test QC scripts with sample data
#
# This script demonstrates how to use the QC scripts by creating test data
# and running through a complete workflow.
#
# Run this to verify all scripts are working correctly!
#
# Usage:
#     bash test_scripts.sh
#
################################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Create test directory
TEST_DIR="/tmp/qc_test_$$"
mkdir -p "$TEST_DIR"

echo -e "${GREEN}QC Scripts Test Suite${NC}"
echo "Test directory: $TEST_DIR"
echo ""

# ============================================================================
# Create test FASTQ data
# ============================================================================

echo -e "${YELLOW}Creating test data...${NC}"

# Create a simple test FASTQ file with correct quality/sequence lengths
cat > "$TEST_DIR/test_input.fastq" << 'EOF'
@sequence_1 sample=test1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@sequence_2 sample=test2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@sequence_3 sample=test3_lowquality
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
###########33333333333333333333333333333333333333333333333333333333333333333333
@sequence_4 sample=test4_short
ATCGATCGATCG
+
IIIIIIII7777
@sequence_5 sample=test5_with_n
ATCGATCGATCGATCGATCGATNNNNNGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIII!!!!!!!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

echo -e "${GREEN}✓ Created test FASTQ file${NC}"

# Create corresponding FASTA and QUAL files
cat > "$TEST_DIR/test_input.fasta" << 'EOF'
>sequence_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>sequence_2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>sequence_3_lowquality
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>sequence_4_short
ATCGATCGATCG
>sequence_5_with_n
ATCGATCGATCGATCGATCGATNNNNNGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > "$TEST_DIR/test_input.qual" << 'EOF'
>sequence_1
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40
>sequence_2
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40
>sequence_3_lowquality
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20
>sequence_4_short
40 40 40 40 40 40 40 40 40 40 40 40
>sequence_5_with_n
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 5 5 5 5 5 5 5 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40
EOF

echo -e "${GREEN}✓ Created test FASTA and QUAL files${NC}"

# ============================================================================
# Test filter_sequences.py
# ============================================================================

echo ""
echo -e "${YELLOW}Testing filter_sequences.py...${NC}"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if ! python3 "$SCRIPT_DIR/filter_sequences.py" "$TEST_DIR/test_input.fastq" \
    --min-length 50 \
    --min-quality 15 \
    --max-n 10 \
    --output "$TEST_DIR/test_filtered.fastq"; then
    echo -e "${RED}✗ filter_sequences.py test failed${NC}"
    exit 1
fi

echo -e "${GREEN}✓ filter_sequences.py test passed${NC}"

if [ -f "$TEST_DIR/test_filtered.fastq" ]; then
    echo -e "${GREEN}✓ Output file created${NC}"
    FILTERED_COUNT=$(grep -c "^@" "$TEST_DIR/test_filtered.fastq" || true)
    echo "  Sequences after filtering: $FILTERED_COUNT"
else
    echo -e "${RED}✗ Output file not created${NC}"
    exit 1
fi

# ============================================================================
# Test trim_quality.py with FASTQ
# ============================================================================

echo ""
echo -e "${YELLOW}Testing trim_quality.py (FASTQ)...${NC}"

if ! python3 "$SCRIPT_DIR/trim_quality.py" "$TEST_DIR/test_input.fastq" \
    --quality 10 \
    --output "$TEST_DIR/test_trimmed_fastq.fastq"; then
    echo -e "${RED}✗ trim_quality.py FASTQ test failed${NC}"
    exit 1
fi

echo -e "${GREEN}✓ trim_quality.py FASTQ test passed${NC}"

if [ -f "$TEST_DIR/test_trimmed_fastq.fastq" ]; then
    echo -e "${GREEN}✓ Output file created${NC}"
    TRIMMED_COUNT=$(grep -c "^@" "$TEST_DIR/test_trimmed_fastq.fastq" || true)
    echo "  Sequences after trimming: $TRIMMED_COUNT"
else
    echo -e "${RED}✗ Output file not created${NC}"
    exit 1
fi

# ============================================================================
# Test trim_quality.py with FASTA+QUAL
# ============================================================================

echo ""
echo -e "${YELLOW}Testing trim_quality.py (FASTA+QUAL)...${NC}"

if ! python3 "$SCRIPT_DIR/trim_quality.py" "$TEST_DIR/test_input.fasta" \
    --qual "$TEST_DIR/test_input.qual" \
    --quality 10 \
    --format fasta \
    --output "$TEST_DIR/test_trimmed_fasta.fasta"; then
    echo -e "${RED}✗ trim_quality.py FASTA test failed${NC}"
    exit 1
fi

echo -e "${GREEN}✓ trim_quality.py FASTA+QUAL test passed${NC}"

if [ -f "$TEST_DIR/test_trimmed_fasta.fasta" ]; then
    echo -e "${GREEN}✓ Output file created${NC}"
    TRIMMED_FASTA_COUNT=$(grep -c "^>" "$TEST_DIR/test_trimmed_fasta.fasta" || true)
    echo "  Sequences after trimming: $TRIMMED_FASTA_COUNT"
else
    echo -e "${RED}✗ Output file not created${NC}"
    exit 1
fi

# ============================================================================
# Test parse_ab1.py (help only - we don't have real AB1 files)
# ============================================================================

echo ""
echo -e "${YELLOW}Testing parse_ab1.py (help check)...${NC}"

if ! python3 "$SCRIPT_DIR/parse_ab1.py" -h > /dev/null 2>&1; then
    echo -e "${RED}✗ parse_ab1.py help test failed${NC}"
    exit 1
fi

echo -e "${GREEN}✓ parse_ab1.py help test passed${NC}"

# ============================================================================
# Test batch_qc.sh (syntax check only)
# ============================================================================

echo ""
echo -e "${YELLOW}Testing batch_qc.sh (syntax check)...${NC}"

if ! bash -n "$SCRIPT_DIR/batch_qc.sh"; then
    echo -e "${RED}✗ batch_qc.sh syntax check failed${NC}"
    exit 1
fi

echo -e "${GREEN}✓ batch_qc.sh syntax check passed${NC}"

# ============================================================================
# Summary
# ============================================================================

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}All tests passed successfully!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Test summary:"
echo "  ✓ filter_sequences.py: Working"
echo "  ✓ trim_quality.py (FASTQ): Working"
echo "  ✓ trim_quality.py (FASTA+QUAL): Working"
echo "  ✓ parse_ab1.py: Syntax OK"
echo "  ✓ batch_qc.sh: Syntax OK"
echo ""
echo "Test results:"
echo "  Original sequences: 5"
echo "  After Q15 filtering (min-length 50): $FILTERED_COUNT"
echo "  After Q10 trimming: $TRIMMED_COUNT"
echo ""
echo "Test data location: $TEST_DIR"
echo "Feel free to explore the output files:"
echo "  cat $TEST_DIR/test_filtered.fastq"
echo "  cat $TEST_DIR/test_trimmed_fastq.fastq"
echo "  cat $TEST_DIR/test_trimmed_fasta.fasta"
echo ""
echo -e "${GREEN}Ready to use QC scripts!${NC}"

exit 0
