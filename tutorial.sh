#!/bin/bash
# DNA Barcoding Analysis - Simple Workflow
# Copy your .ab1 files to data/test_data/ then run this script

set -e

# No colors for clean logs
GREEN=''
BLUE=''
YELLOW=''
NC=''

# Create output directories automatically
mkdir -p results/tutorial/{01_qc,02_consensus,03_alignment,04_phylogeny,05_blast}

# Start logging
LOG_FILE="results/tutorial/tutorial.log"
exec > >(tee "$LOG_FILE") 2>&1
echo "=== Tutorial Pipeline Started: $(date) ==="

echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║                   DNA BARCODING ANALYSIS                             ║"
echo "║                                                                      ║"
echo "║  This script runs all 5 steps automatically on your sequences       ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""
echo "REQUIREMENTS:"
echo "  1. Your .ab1 files must be in: data/test_data/"
echo "  2. You must have completed the tutorial first: ./tutorial.sh"
echo ""
echo "Starting tutorial..."

CONTAINER="cosmelab/dna-barcoding-analysis:latest"
OUTPUT_DIR="results/tutorial"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 1 of 5: Quality Control"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  $CONTAINER \
  python3 modules/01_quality_control/qc_chromatograms.py \
  data/test_data/ \
  $OUTPUT_DIR/01_qc/

echo "✓ Step 1 complete!"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 2 of 5: Create Consensus Sequences"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  $CONTAINER \
  python3 modules/02_consensus/create_consensus.py \
  $OUTPUT_DIR/01_qc/passed_sequences.fasta \
  $OUTPUT_DIR/02_consensus/ \
  --pairs-only

echo "✓ Step 2 complete!"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 3 of 5: Combine with Reference Sequences"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

cat $OUTPUT_DIR/02_consensus/consensus_sequences.fasta \
    data/reference_sequences/socal_mosquitoes.fasta \
    > $OUTPUT_DIR/02_consensus/combined_with_references.fasta

NUM_SEQS=$(grep -c "^>" $OUTPUT_DIR/02_consensus/combined_with_references.fasta)
echo "Combined $NUM_SEQS sequences (your samples + references)"

echo "✓ Step 3 complete!"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 4 of 5: Alignment & Phylogenetic Tree"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
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

echo "✓ Step 4 complete!"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 5 of 5: Species Identification (BLAST)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  $CONTAINER \
  python3 modules/05_identification/identify_species.py \
  $OUTPUT_DIR/02_consensus/consensus_sequences.fasta \
  $OUTPUT_DIR/05_blast/

echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║                     ✓ ANALYSIS COMPLETE! ✓                          ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
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
echo "=== Tutorial Pipeline Completed: $(date) ==="
echo "Log saved to: $LOG_FILE"
