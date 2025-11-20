#!/bin/bash
# DNA Barcoding Tutorial - Learn the Complete Workflow
# This tutorial uses TEST DATA so you can learn without worrying about mistakes!

set -e  # Exit on error

# Color codes for pretty output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                    🧬 DNA BARCODING TUTORIAL 🧬                      ║
║                                                                      ║
║              Learn the 5-Step Workflow with Test Data                ║
╚══════════════════════════════════════════════════════════════════════╝

WELCOME!

This tutorial will teach you how to:
  ✓ Analyze Sanger sequencing chromatograms
  ✓ Create consensus sequences from F+R reads
  ✓ Build phylogenetic trees
  ✓ Identify mosquito species using DNA barcoding

We're using TEST DATA, so you can't break anything. Have fun! 🎉

TIME: About 15-20 minutes

EOF

read -p "Press ENTER to start the tutorial..."

# =============================================================================
# STEP 1: QUALITY CONTROL
# =============================================================================

clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                           STEP 1 of 5                                ║
║                        📊 QUALITY CONTROL 📊                         ║
╚══════════════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHAT THIS STEP DOES:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

When you get DNA sequences back from the UC Genomics Core, they send you
.ab1 files (chromatograms). Not all sequences are perfect! Some might be:

  ✗ Too short (less than 500 base pairs)
  ✗ Low quality (hard to read the bases)
  ✗ Failed completely (no readable sequence)

This step analyzes each chromatogram and tells you which ones are good
enough to use for species identification.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 YOUR TEST DATA:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

You have 8 files in data/test_data/:
  • AT99F_C01_011.ab1    (forward read, sample AT99)
  • AT99R_G01_003.ab1    (reverse read, sample AT99)
  • AT_ROCK_F_D01_009.ab1
  • AT_ROCK_R_H01_001.ab1
  • AT83F_A01_015.ab1
  • AT83R_E01_007.ab1
  • AT94F_B01_013.ab1
  • AT94R_F01_005.ab1

That's 4 samples, each with forward (F) and reverse (R) reads.

EOF

read -p "Ready to check the quality? Press ENTER to run QC..."

echo ""
echo -e "${BLUE}Running QC command...${NC}"
echo ""

docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/01_quality_control/qc_chromatograms.py \
  data/test_data/ \
  results/tutorial/01_qc/ \
  --open

echo ""
echo -e "${GREEN}✓ QC Complete!${NC}"
echo ""
echo "LOOK AT THE HTML REPORT that just opened in your browser."
echo "You should see:"
echo "  • Chromatogram visualizations (the peaks are the DNA signal)"
echo "  • Quality scores for each sequence"
echo "  • PASS/FAIL status for each file"
echo ""
read -p "When you've looked at the report, press ENTER for Step 2..."

# =============================================================================
# STEP 2: CONSENSUS SEQUENCES
# =============================================================================

clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                           STEP 2 of 5                                ║
║                      🎯 CONSENSUS SEQUENCES 🎯                       ║
╚══════════════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHAT THIS STEP DOES:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Why do we sequence both FORWARD and REVERSE?

  🧬 Forward read:  →→→→→→→→→→  (reads 5' to 3' direction)
  🧬 Reverse read:  ←←←←←←←←←←  (reads the other strand)

By combining both directions, we get a CONSENSUS sequence that is more
accurate than either read alone. Think of it like double-checking your work!

The --pairs-only flag means: "Only keep samples where BOTH F and R passed QC"

If a sample only has F or only has R, we skip it. We want the best data!

EOF

read -p "Ready to create consensus sequences? Press ENTER..."

echo ""
echo -e "${BLUE}Running consensus command...${NC}"
echo ""

docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/02_consensus/create_consensus.py \
  results/tutorial/01_qc/passed_sequences.fasta \
  results/tutorial/02_consensus/ \
  --pairs-only \
  --open

echo ""
echo -e "${GREEN}✓ Consensus sequences created!${NC}"
echo ""
echo "LOOK AT THE HTML REPORT."
echo "You should see which samples have complete F+R pairs."
echo ""
read -p "When you've looked at the report, press ENTER for Step 3..."

# =============================================================================
# STEP 3: COMBINE WITH REFERENCES
# =============================================================================

clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                           STEP 3 of 5                                ║
║                  📚 COMBINE WITH REFERENCE SEQUENCES 📚              ║
╚══════════════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHAT THIS STEP DOES:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

To build a phylogenetic tree, we need REFERENCE SEQUENCES - DNA from
mosquitoes we already know the species of.

We'll combine YOUR consensus sequences with 52 reference sequences from
known Southern California mosquito species.

This lets you see WHERE your samples cluster on the tree:
  • Do they cluster with Aedes aegypti?
  • Do they cluster with Culex pipiens?
  • Do they cluster with something else?

EOF

read -p "Ready to combine sequences? Press ENTER..."

echo ""
echo -e "${BLUE}Combining your sequences with reference database...${NC}"
echo ""

cat results/tutorial/02_consensus/consensus_sequences.fasta \
    data/reference_sequences/socal_mosquitoes.fasta \
    > results/tutorial/02_consensus/combined_with_references.fasta

NUM_SEQS=$(grep -c "^>" results/tutorial/02_consensus/combined_with_references.fasta)

echo -e "${GREEN}✓ Combined sequences successfully!${NC}"
echo ""
echo "Created file with ${NUM_SEQS} total sequences"
echo "  • Your consensus sequences"
echo "  • 52 reference sequences from known mosquito species"
echo ""
read -p "Press ENTER for Step 4 (Alignment & Tree)..."

# =============================================================================
# STEP 4: ALIGNMENT AND PHYLOGENETIC TREE
# =============================================================================

clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                           STEP 4 of 5                                ║
║                   📐 ALIGNMENT & PHYLOGENETIC TREE 🌳                ║
╚══════════════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHAT THIS STEP DOES:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

PART 4A: ALIGNMENT (MAFFT)
  • Lines up all sequences so we can compare them
  • Adds gaps (-) where sequences are different lengths
  • Shows conserved vs. variable positions

PART 4B: PHYLOGENETIC TREE (IQ-TREE)
  • Builds an evolutionary tree showing relationships
  • Uses maximum likelihood to find the best tree
  • Adds bootstrap values (confidence scores)

The tree shows you which mosquito species your samples are most closely
related to!

⏱️  This step takes ~2-3 minutes. Time for a coffee break! ☕

EOF

read -p "Ready to build the tree? Press ENTER..."

echo ""
echo -e "${BLUE}Step 4A: Running MAFFT alignment...${NC}"
echo ""

docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/03_alignment/align_sequences.py \
  results/tutorial/02_consensus/combined_with_references.fasta \
  results/tutorial/03_alignment/

echo ""
echo -e "${GREEN}✓ Alignment complete!${NC}"
echo ""
echo -e "${BLUE}Step 4B: Running IQ-TREE (building phylogenetic tree)...${NC}"
echo ""

docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/04_phylogeny/build_tree.py \
  results/tutorial/03_alignment/aligned_sequences.fasta \
  results/tutorial/04_phylogeny/

echo ""
echo -e "${GREEN}✓ Tree construction complete!${NC}"
echo ""
echo "Your tree has been saved as:"
echo "  • results/tutorial/04_phylogeny/tree.png (simple visualization)"
echo "  • results/tutorial/04_phylogeny/tree.treefile (for FigTree)"
echo ""
echo "LOOK AT THE TREE! See where your samples cluster with known species."
echo ""
read -p "When you've looked at the tree, press ENTER for Step 5 (BLAST)..."

# =============================================================================
# STEP 5: SPECIES IDENTIFICATION (BLAST)
# =============================================================================

clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                           STEP 5 of 5                                ║
║                   🏷️  SPECIES IDENTIFICATION (BLAST) 🏷️             ║
╚══════════════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHAT THIS STEP DOES:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

BLAST (Basic Local Alignment Search Tool) compares your DNA sequences to
the NCBI GenBank database - millions of sequences from around the world.

It will tell you:
  • What species matches your sequence
  • How similar they are (% identity)
  • How confident we can be in the match

INTERPRETING RESULTS:
  • >98% identity → Same species
  • 95-98% identity → Same genus, possibly different species
  • <95% identity → Different genus or poor quality sequence

⏱️  This takes ~30 seconds (NCBI limits how fast we can query)

EOF

read -p "Ready to identify species? Press ENTER..."

echo ""
echo -e "${BLUE}Running BLAST queries...${NC}"
echo ""

docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/05_identification/identify_species.py \
  results/tutorial/02_consensus/consensus_sequences.fasta \
  results/tutorial/05_blast/

echo ""
echo -e "${GREEN}✓ Species identification complete!${NC}"
echo ""
echo "LOOK AT THE HTML REPORT that shows:"
echo "  • Top 5 species matches for each sample"
echo "  • % identity scores"
echo "  • Scientific names and common names"
echo ""
read -p "When you've looked at the BLAST results, press ENTER to finish..."

# =============================================================================
# TUTORIAL COMPLETE!
# =============================================================================

clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                    🎉 TUTORIAL COMPLETE! 🎉                          ║
╚══════════════════════════════════════════════════════════════════════╝

CONGRATULATIONS! You've learned the complete DNA barcoding workflow:

  ✓ Step 1: Quality Control
  ✓ Step 2: Consensus Sequences
  ✓ Step 3: Combine with References
  ✓ Step 4: Alignment & Phylogenetic Tree
  ✓ Step 5: Species Identification (BLAST)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 YOUR RESULTS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

All your results are saved in: results/tutorial/

  📊 QC report:        results/tutorial/01_qc/qc_report.html
  🎯 Consensus report: results/tutorial/02_consensus/consensus_report.html
  📐 Alignment report: results/tutorial/03_alignment/alignment_report.html
  🌳 Tree files:       results/tutorial/04_phylogeny/
  🏷️  BLAST results:   results/tutorial/05_blast/identification_report.html

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 NEXT STEPS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Now you're ready to analyze YOUR OWN sequences!

1. Put your .ab1 files in: data/student_sequences/

2. Run the same 5 steps (replace data/test_data with data/student_sequences)

3. Compare your tree to the tutorial tree - do you have similar species?

4. Fill in the assignment worksheet with your results

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 NEED HELP?
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

• Read: docs/pipeline_workflow.md (visual guide)
• Read: docs/iqtree_guide.md (understanding trees)
• Ask your TA or instructor

Good luck with your analysis! 🧬🦟

EOF

echo ""
echo -e "${GREEN}Tutorial complete!${NC}"
echo ""
