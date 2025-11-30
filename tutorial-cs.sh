#!/bin/bash
# DNA Barcoding Analysis - Tutorial Workflow (Codespaces Version)
# Learn the pipeline with test data before analyzing your own

set -e

# Suppress Python deprecation warnings (BioPython pairwise2)
export PYTHONWARNINGS="ignore::DeprecationWarning"

OUTPUT_DIR="results/tutorial"

# Create output directories automatically
mkdir -p results/tutorial/{01_qc,02_consensus,03_alignment,04_phylogeny,05_blast}

# Start logging
LOG_FILE="results/tutorial/tutorial.log"
exec > >(tee "$LOG_FILE") 2>&1
echo "=== Tutorial Pipeline Started ==="
echo "☁️  Environment: GitHub Codespaces"
echo "📝 Log file: $LOG_FILE"
echo ""

# Show beautiful welcome banner using Rich
python3 modules/show_welcome.py --welcome --mode tutorial

echo ""
echo "This tutorial uses TEST DATA in data/test_data/"
echo "Learn the workflow here, then run ./run-analysis-cs.sh for your data!"
echo ""
echo "Starting tutorial..."

# Show Step 1 banner
python3 modules/show_welcome.py --step 1

python3 modules/01_quality_control/qc_chromatograms.py \
  data/test_data/ \
  $OUTPUT_DIR/01_qc/

echo "✓ Step 1 complete!"

# Show Step 2 banner
python3 modules/show_welcome.py --step 2

python3 -W ignore::DeprecationWarning modules/02_consensus/create_consensus.py \
  $OUTPUT_DIR/01_qc/passed_sequences.fasta \
  $OUTPUT_DIR/02_consensus/ \
  --pairs-only 2>&1 | grep -v "BiopythonDeprecationWarning\|Bio.pairwise2\|warnings.warn"

echo "✓ Step 2 complete!"

# Show Step 3 banner
python3 modules/show_welcome.py --step 3

cat $OUTPUT_DIR/02_consensus/consensus_sequences.fasta \
    data/reference_sequences/socal_mosquitoes.fasta \
    > $OUTPUT_DIR/02_consensus/combined_with_references.fasta

NUM_SEQS=$(grep -c "^>" $OUTPUT_DIR/02_consensus/combined_with_references.fasta)
echo "Combined $NUM_SEQS sequences (your samples + references)"

echo "✓ Step 3 complete!"

# Show Step 4 banner
python3 modules/show_welcome.py --step 4

echo ""

# Alignment
python3 modules/03_alignment/align_sequences.py \
  $OUTPUT_DIR/02_consensus/combined_with_references.fasta \
  $OUTPUT_DIR/03_alignment/

echo ""
# Tree
python3 modules/04_phylogeny/build_tree.py \
  $OUTPUT_DIR/03_alignment/aligned_sequences.fasta \
  $OUTPUT_DIR/04_phylogeny/

echo "✓ Step 4 complete!"

# Show Step 5 banner
python3 modules/show_welcome.py --step 5

python3 modules/05_identification/identify_species.py \
  $OUTPUT_DIR/02_consensus/consensus_sequences.fasta \
  $OUTPUT_DIR/05_blast/

# Show completion banner using Rich
python3 modules/show_welcome.py --complete

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
echo "=== Tutorial Pipeline Completed ==="

# Clean ANSI color codes from log file for readability in text editors
sed 's/\x1b\[[0-9;]*m//g' "$LOG_FILE" > "${LOG_FILE}.clean" && mv "${LOG_FILE}.clean" "$LOG_FILE"
echo "Log saved to: $LOG_FILE"
