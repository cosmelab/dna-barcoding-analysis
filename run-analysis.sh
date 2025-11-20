#!/bin/bash
# DNA Barcoding Analysis - Simple Workflow
# Copy your .ab1 files to data/student_sequences/ then run this script

set -e

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

clear
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║                   DNA BARCODING ANALYSIS                             ║"
echo "║                                                                      ║"
echo "║  This script runs all 5 steps automatically on your sequences       ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""
echo "REQUIREMENTS:"
echo "  1. Your .ab1 files must be in: data/student_sequences/"
echo "  2. You must have completed the tutorial first: ./tutorial.sh"
echo ""
read -p "Press ENTER to start analysis..."

CONTAINER="cosmelab/dna-barcoding-analysis:latest"
OUTPUT_DIR="results/my_analysis"

echo ""
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}STEP 1 of 5: Quality Control${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  $CONTAINER \
  python3 modules/01_quality_control/qc_chromatograms.py \
  data/student_sequences/ \
  $OUTPUT_DIR/01_qc/ \
  --open

echo -e "${GREEN}✓ Step 1 complete!${NC}"
echo ""
read -p "Press ENTER for Step 2..."

echo ""
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}STEP 2 of 5: Create Consensus Sequences${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  $CONTAINER \
  python3 modules/02_consensus/create_consensus.py \
  $OUTPUT_DIR/01_qc/passed_sequences.fasta \
  $OUTPUT_DIR/02_consensus/ \
  --pairs-only \
  --open

echo -e "${GREEN}✓ Step 2 complete!${NC}"
echo ""
read -p "Press ENTER for Step 3..."

echo ""
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}STEP 3 of 5: Combine with Reference Sequences${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

cat $OUTPUT_DIR/02_consensus/consensus_sequences.fasta \
    data/reference_sequences/socal_mosquitoes.fasta \
    > $OUTPUT_DIR/02_consensus/combined_with_references.fasta

NUM_SEQS=$(grep -c "^>" $OUTPUT_DIR/02_consensus/combined_with_references.fasta)
echo "Combined $NUM_SEQS sequences (your samples + references)"

echo -e "${GREEN}✓ Step 3 complete!${NC}"
echo ""
read -p "Press ENTER for Step 4..."

echo ""
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}STEP 4 of 5: Alignment & Phylogenetic Tree${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo "This step takes ~2-3 minutes..."
echo ""

# Alignment
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  $CONTAINER \
  python3 modules/03_alignment/align_sequences.py \
  $OUTPUT_DIR/02_consensus/combined_with_references.fasta \
  $OUTPUT_DIR/03_alignment/

echo ""
# Tree
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  $CONTAINER \
  python3 modules/04_phylogeny/build_tree.py \
  $OUTPUT_DIR/03_alignment/aligned_sequences.fasta \
  $OUTPUT_DIR/04_phylogeny/

echo -e "${GREEN}✓ Step 4 complete!${NC}"
echo ""
read -p "Press ENTER for Step 5 (final step)..."

echo ""
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}STEP 5 of 5: Species Identification (BLAST)${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  $CONTAINER \
  python3 modules/05_identification/identify_species.py \
  $OUTPUT_DIR/02_consensus/consensus_sequences.fasta \
  $OUTPUT_DIR/05_blast/

echo ""
echo -e "${GREEN}╔══════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║                     ✓ ANALYSIS COMPLETE! ✓                          ║${NC}"
echo -e "${GREEN}╚══════════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Your results are in: $OUTPUT_DIR/"
echo ""
echo "  01_qc/          - Quality control report"
echo "  02_consensus/   - Consensus sequences"
echo "  03_alignment/   - Sequence alignment"
echo "  04_phylogeny/   - Phylogenetic tree"
echo "  05_blast/       - Species identification"
echo ""
echo "Open these files:"
echo "  • $OUTPUT_DIR/01_qc/qc_report.html"
echo "  • $OUTPUT_DIR/04_phylogeny/tree.png"
echo "  • $OUTPUT_DIR/05_blast/identification_report.html"
echo ""
