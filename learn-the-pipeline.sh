#!/bin/bash
# DNA Barcoding Tutorial - Learn with Test Data
# This script walks you through the entire pipeline step-by-step

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Clear screen and show welcome
clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                                                                      ║
║                    🧬 DNA BARCODING TUTORIAL 🧬                      ║
║                                                                      ║
║              Learning the Pipeline with Test Data                   ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝
EOF

echo ""
echo "Welcome! This tutorial will teach you how to:"
echo ""
echo "  📊 Check sequence quality"
echo "  🧬 Align DNA sequences"
echo "  🌳 Build phylogenetic trees"
echo "  🏷️  Identify species with BLAST"
echo ""
echo "We'll use TEST DATA first so you can see how everything works."
echo "Then you'll analyze YOUR OWN DATA from the lab."
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  This takes about 15 minutes. Relax and learn! ☕"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
read -p "Press ENTER when you're ready to start..."

# ============================================================================
# STEP 1: QUALITY CONTROL
# ============================================================================
clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                         STEP 1 of 4                                  ║
║                      📊 QUALITY CONTROL 📊                           ║
╚══════════════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHAT THIS STEP DOES:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DNA sequencing isn't perfect! Some sequences might be:
  ✗ Too short (less than 500 base pairs)
  ✗ Low quality (hard to read the bases)
  ✗ Lots of ambiguous bases (N's instead of A/T/G/C)

This step checks each sequence and marks it as PASS or FAIL.

Only sequences that PASS will be used in the next steps.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHY IT MATTERS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Bad sequences → Wrong species identification
Good sequences → Confident results

It's NORMAL for some sequences to fail. Even professional labs get
50-80% success rates from Sanger sequencing.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
EOF

echo ""
read -p "Press ENTER to run Quality Control on 8 test samples..."

echo ""
echo "Running QC..."
python3 modules/01_quality_control/qc_chromatograms.py data/test_data results/tutorial/qc --open

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "${GREEN}✓ STEP 1 COMPLETE!${NC}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📊 The QC report should have opened in your web browser."
echo ""
echo "WHAT TO LOOK FOR IN THE REPORT:"
echo "  • How many sequences PASSED? (should be 4)"
echo "  • How many sequences FAILED? (should be 4)"
echo "  • Why did they fail? (click on red rows to see reasons)"
echo "  • Look at the chromatograms - see the colored peaks?"
echo "  • Notice that forward and reverse reads are grouped together"
echo ""
echo "QUESTIONS TO THINK ABOUT:"
echo "  • Why might a sequence fail quality control?"
echo "  • Would you trust a sequence with average quality score of 15?"
echo "  • What does the chromatogram tell you?"
echo ""
read -p "When you've looked at the report, press ENTER for Step 2..."

# ============================================================================
# STEP 2: SEQUENCE ALIGNMENT
# ============================================================================
clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                         STEP 2 of 4                                  ║
║                    🧬 SEQUENCE ALIGNMENT 🧬                          ║
╚══════════════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHAT THIS STEP DOES:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Your sequences aren't all exactly the same length. Alignment "lines them up"
so we can compare them position-by-position.

BEFORE alignment:                AFTER alignment:
>Sample1                         >Sample1
ATCGATCGATCG                     ATCGAT---CGATCG
>Sample2                         >Sample2
ATCGATCGATCG                     ATCGATCGATCG---
>Sample3                         >Sample3
ATCGATCG                         ATCGAT---CG----

The dashes (---) are gaps added to make everything line up.

We use MAFFT, a tool that's very good at finding the best alignment.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHY IT MATTERS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

To compare DNA sequences, we need to know which position in one sequence
corresponds to which position in another sequence. Alignment does this!

After alignment, you can see:
  • Which positions are CONSERVED (same in all sequences)
  • Which positions are VARIABLE (different between species)
  • Where insertions or deletions occurred

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
EOF

echo ""
read -p "Press ENTER to align the 4 sequences that passed QC..."

echo ""
echo "Running MAFFT alignment..."
python3 modules/02_alignment/align_sequences.py results/tutorial/qc/passed_sequences.fasta results/tutorial/alignment --open

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "${GREEN}✓ STEP 2 COMPLETE!${NC}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📊 The alignment report should have opened in your browser."
echo ""
echo "WHAT TO LOOK FOR IN THE REPORT:"
echo "  • What's the alignment length? (should be ~785 bp)"
echo "  • In the visual alignment:"
echo "    - UPPERCASE letters = conserved positions (same in all sequences)"
echo "    - lowercase letters = variable positions (different between sequences)"
echo "  • Notice the color coding for different bases (A/T/G/C)"
echo "  • See where gaps (---) were added"
echo ""
echo "QUESTIONS TO THINK ABOUT:"
echo "  • Why are some positions conserved and others variable?"
echo "  • What do the gaps represent biologically?"
echo "  • Which species are most similar (have fewest differences)?"
echo ""
read -p "When you've looked at the alignment, press ENTER for Step 3..."

# ============================================================================
# STEP 3: PHYLOGENETIC TREE
# ============================================================================
clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                         STEP 3 of 4                                  ║
║                   🌳 PHYLOGENETIC TREE 🌳                            ║
╚══════════════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHAT THIS STEP DOES:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

IQ-TREE builds an evolutionary "family tree" showing how your sequences
are related to each other.

Think of it like a genealogy tree:
  • Closely related organisms = close together on the tree
  • Distantly related organisms = far apart on the tree

Example tree:
                ┌─ Your_Sample_1
      ┌─100%───┤
      │         └─ Your_Sample_2  ← These are closely related
  ────┤
      │         ┌─ Your_Sample_3
      └─95%────┤
                └─ Your_Sample_4  ← These are also close

The numbers (100%, 95%) are BOOTSTRAP VALUES = how confident we are.
  • >90% = Very confident
  • 70-90% = Moderately confident
  • <70% = Not very confident

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHY IT MATTERS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The tree helps you understand:
  • Which samples are from the same species (cluster together)
  • Which samples are from different species (separate branches)
  • How confident we are in these relationships (bootstrap values)

When you add reference sequences (from GenBank), the tree shows:
  • Which reference species your samples cluster with
  • This helps identify your unknown samples!

NOTE: We only have 4 sequences from YOUR samples. To really identify
species, we'll use BLAST in Step 4. But the tree is still useful to see
relationships!

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
EOF

echo ""
read -p "Press ENTER to build the phylogenetic tree with IQ-TREE..."

echo ""
echo "Running IQ-TREE (this takes ~30 seconds)..."
python3 modules/03_phylogeny/build_tree.py results/tutorial/alignment/aligned_sequences.fasta results/tutorial/phylogeny --open

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "${GREEN}✓ STEP 3 COMPLETE!${NC}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📊 The tree report should have opened in your browser."
echo ""
echo "WHAT TO LOOK FOR IN THE REPORT:"
echo "  • tree.png - A simple visualization of your tree"
echo "  • Which samples cluster together?"
echo "  • What are the bootstrap values? (should be >70% for reliable branches)"
echo "  • Branch lengths show how different sequences are"
echo ""
echo "ADVANCED (OPTIONAL):"
echo "  • Download FigTree: http://tree.bio.ed.ac.uk/software/figtree/"
echo "  • Open results/tutorial/phylogeny/tree.treefile"
echo "  • Make a publication-quality figure!"
echo ""
echo "QUESTIONS TO THINK ABOUT:"
echo "  • Do your forward and reverse reads cluster together? (they should!)"
echo "  • What does a long branch vs. short branch mean?"
echo "  • Why might bootstrap values be low with only 4 sequences?"
echo ""
read -p "When you've looked at the tree, press ENTER for Step 4..."

# ============================================================================
# STEP 4: SPECIES IDENTIFICATION
# ============================================================================
clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                         STEP 4 of 4                                  ║
║                🏷️  SPECIES IDENTIFICATION 🏷️                        ║
╚══════════════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHAT THIS STEP DOES:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

BLAST (Basic Local Alignment Search Tool) compares your sequences to
millions of sequences in GenBank (NCBI's database).

It finds the closest matches and tells you:
  • What species they match
  • How good the match is (% identity)
  • How confident we are (e-value)

Example results:
  Your_Sample_1 → Aedes aegypti (98.5% identity)
  Your_Sample_2 → Culex quinquefasciatus (99.2% identity)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 HOW TO INTERPRET RESULTS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

% Identity Guidelines:
  • 98-100%    = Almost certainly the same species
  • 95-97%     = Probably the same species (could be subspecies)
  • 90-95%     = Same genus, different species
  • 85-90%     = Related genus
  • <85%       = Distantly related (be cautious!)

IMPORTANT: Always check multiple top hits!
  • If top 3 hits are all the same species → confident ID
  • If top hits are different species → less confident

BLAST searches the internet, so this takes ~30 seconds per sequence.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
EOF

echo ""
read -p "Press ENTER to identify species with BLAST (takes ~2 minutes)..."

echo ""
echo "Running BLAST searches against GenBank..."
echo "(This searches the internet, so it takes a bit longer)"
python3 modules/04_identification/identify_species.py results/tutorial/qc/passed_sequences.fasta results/tutorial/blast --open

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "${GREEN}✓ STEP 4 COMPLETE!${NC}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📊 The BLAST results should have opened in your browser."
echo ""
echo "WHAT TO LOOK FOR IN THE REPORT:"
echo "  • What species is each sequence?"
echo "  • What's the % identity for the top hit?"
echo "  • Do the top 3 hits agree on the species?"
echo "  • Are the results believable for Southern California mosquitoes?"
echo ""
echo "QUESTIONS TO THINK ABOUT:"
echo "  • Why might forward and reverse reads get slightly different matches?"
echo "  • What % identity would you need to confidently call it a species?"
echo "  • What if your top hit was 85%? Would you trust it?"
echo ""
read -p "Press ENTER to see the final summary..."

# ============================================================================
# TUTORIAL COMPLETE
# ============================================================================
clear
cat << "EOF"
╔══════════════════════════════════════════════════════════════════════╗
║                                                                      ║
║                  🎉 TUTORIAL COMPLETE! 🎉                            ║
║                                                                      ║
║          You just learned the DNA barcoding pipeline!                ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 WHAT YOU LEARNED:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

 ✓ STEP 1: Quality Control
   → Separate good sequences from bad ones
   → Check length, quality, ambiguous bases

 ✓ STEP 2: Sequence Alignment
   → Line up sequences for comparison
   → Find conserved and variable regions

 ✓ STEP 3: Phylogenetic Tree
   → Show evolutionary relationships
   → Cluster similar sequences together

 ✓ STEP 4: Species Identification
   → BLAST against GenBank database
   → Find closest matching species

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 YOUR TUTORIAL RESULTS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

All results are saved in: results/tutorial/

  📊 QC Report:       results/tutorial/qc/qc_report.html
  🧬 Alignment:       results/tutorial/alignment/alignment_report.html
  🌳 Tree:            results/tutorial/phylogeny/tree.png
  🏷️  Species IDs:    results/tutorial/blast/identification_report.html

You can re-open these anytime to review!

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 NEXT STEP: Analyze YOUR real data!
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1. Get your .ab1 files from UC Genomics Core

2. Put them in: data/student_sequences/

3. Run this command:
   ./analyze-my-data.sh

4. Answer the assignment questions using your results!

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 NEED HELP?
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  • Re-run this tutorial anytime: ./learn-the-pipeline.sh
  • Read the IQ-TREE guide: docs/iqtree_guide.md
  • Check troubleshooting: docs/troubleshooting.md
  • Ask on Canvas discussion board

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Good luck with your mosquito identification! 🦟

EOF
