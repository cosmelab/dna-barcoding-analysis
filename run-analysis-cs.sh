#!/bin/bash
# DNA Barcoding Analysis - Simple Workflow (Codespaces Version)
# Copy your .ab1 files to data/student_sequences/ then run this script

set -e

# Suppress Python deprecation warnings (BioPython pairwise2)
export PYTHONWARNINGS="ignore::DeprecationWarning"

OUTPUT_DIR="results/my_analysis"

# Create output directories automatically
mkdir -p results/my_analysis/{01_qc,02_consensus,03_alignment,04_phylogeny,05_blast}
mkdir -p results/lab_analysis

# Start logging
LOG_FILE="results/my_analysis/student.log"
exec > >(tee "$LOG_FILE") 2>&1
echo "=== Student Analysis Pipeline Started ==="
echo "☁️  Environment: GitHub Codespaces"
echo "📝 Log file: $LOG_FILE"
echo ""

# Show beautiful welcome banner using Rich
python3 modules/show_welcome.py --welcome --mode analysis

echo ""
echo "REQUIREMENTS:"
echo "  1. Your .ab1 files must be in: data/student_sequences/"
echo "  2. You must have completed the tutorial first: ./tutorial-cs.sh"
echo ""
echo "Starting analysis..."

# Show Step 1 banner
python3 modules/show_welcome.py --step 1

python3 modules/01_quality_control/qc_chromatograms.py \
  data/student_sequences/ \
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

echo "✓ Step 5 complete!"

# Show Step 6 banner
python3 modules/show_welcome.py --step 6

echo "Generating lab analysis report (extraction, PCR, sequencing results)..."
echo ""

cd modules/06_lab_data_analysis
python3 plot_extraction.py && \
python3 plot_quality.py && \
python3 plot_pcr.py && \
python3 plot_sequencing.py && \
python3 plot_teams.py && \
python3 generate_report.py
cd ../..

echo "✓ Step 6 complete!"

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
echo "Lab analysis report: results/lab_analysis/"
echo ""
echo "Open these files:"
echo "  • $OUTPUT_DIR/01_qc/qc_report.html"
echo "  • $OUTPUT_DIR/04_phylogeny/tree.png"
echo "  • $OUTPUT_DIR/05_blast/identification_report.html"
echo "  • results/lab_analysis/lab_report.html"
echo ""
echo "=== Student Analysis Pipeline Completed ==="

# Clean ANSI color codes from log file for readability in text editors
sed 's/\x1b\[[0-9;]*m//g' "$LOG_FILE" > "${LOG_FILE}.clean" && mv "${LOG_FILE}.clean" "$LOG_FILE"
echo "Log saved to: $LOG_FILE"
