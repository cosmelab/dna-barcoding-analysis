# Week 8 Lab: DNA Barcoding Analysis

**Due**: End of Week 8
**Points**: 100

---

## 🎯 Objectives

By the end of this lab, you will:
- Analyze Sanger sequencing chromatograms (forward and reverse reads)
- Perform quality control on DNA sequences
- Create consensus sequences from F+R pairs
- Build phylogenetic trees with reference sequences
- Identify mosquito species using COI barcodes

---

## 🖥️ Choose Your Environment

### Option A: GitHub Codespaces (Recommended - No Installation)

1. Click the green **"Code"** button on GitHub
2. Select **"Codespaces"** tab
3. Click **"Create codespace on main"**
4. Wait for the environment to load (~2-3 minutes)
5. You're ready! Use the `-cs` scripts (see below)

### Option B: Local Docker (Advanced)

If you have Docker installed locally, you can run the analysis on your own computer.
See the README for Docker installation instructions.

---

## ⚠️ STEP 0: Complete the Tutorial FIRST (REQUIRED)

**Before analyzing the class data, you MUST complete the interactive tutorial:**

**In Codespaces:**
```bash
./tutorial-cs.sh
```

**With local Docker:**
```bash
./tutorial.sh
```

**Why this is mandatory:**
- ✓ Teaches you all 5 steps of the workflow
- ✓ Uses test data (you can't break anything)
- ✓ Shows you what results should look like
- ✓ Takes only 15-20 minutes
- ✓ Makes the actual assignment much easier!

**Do NOT skip this!** Students who skip the tutorial get confused and make mistakes.

---

## Part 1: Run the Analysis (60 points)

### Class Data

**Important:** Everyone will analyze the **same class dataset**. These are the pooled sequences from all students this semester. Sample names are initials only.

The `.ab1` chromatogram files are in: `data/student_sequences/`

The dataset includes:
- Forward reads (sample names ending in F)
- Reverse reads (sample names ending in R)
- Example: AT-HV1F and AT-HV1R are a pair

---

### 💡 How to Run the Analysis

**OPTION 1 (RECOMMENDED): Use the automated script**

**In Codespaces:**
```bash
./run-analysis-cs.sh
```

**With local Docker:**
```bash
./run-analysis.sh
```

This runs all 6 steps below automatically. **This is the easiest way!**

**OPTION 2: Run each step individually (shown below)**

The commands below show you what happens in each step. You can run them one by one if you want to understand the process better, or if you need to re-run just one step.

**Note:** The individual commands below use Docker. In Codespaces, you can run Python directly (e.g., `python3 modules/01_quality_control/qc_chromatograms.py ...`).

---

### Step 1: Quality Control (10 points)

Check which sequences are good enough to use:

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/01_quality_control/qc_chromatograms.py \
  data/student_sequences/ \
  results/my_analysis/qc/ \
  --open
```

**Look at the HTML report that opens. Count:**
- How many sequences PASSED QC?
- How many samples have BOTH F and R that passed?

### Step 2: Create Consensus Sequences (10 points)

Combine forward and reverse reads:

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/02_consensus/create_consensus.py \
  results/my_analysis/qc/passed_sequences.fasta \
  results/my_analysis/consensus/ \
  --pairs-only \
  --open
```

**The `--pairs-only` flag means:** Only keep samples where BOTH F and R passed QC. This ensures high-quality consensus sequences.

### Step 3: Combine with Reference Sequences (5 points)

Add the class consensus sequences to the database of known SoCal mosquitoes:

```bash
cat results/my_analysis/consensus/consensus_sequences.fasta \
    data/reference_sequences/socal_mosquitoes.fasta \
    > results/my_analysis/consensus/combined_with_references.fasta
```

This creates a file with the CLASS sequences + 52 reference sequences.

### Step 4: Alignment (10 points)

Line up all sequences so we can compare them:

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/03_alignment/align_sequences.py \
  results/my_analysis/consensus/combined_with_references.fasta \
  results/my_analysis/alignment/
```

### Step 5: Phylogenetic Tree (15 points)

Build an evolutionary tree showing relationships:

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/04_phylogeny/build_tree.py \
  results/my_analysis/alignment/aligned_sequences.fasta \
  results/my_analysis/phylogeny/
```

**This takes ~2-3 minutes.** Be patient!

### Step 6: Species Identification (BLAST) (10 points)

Compare your sequences to GenBank database:

```bash
docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest \
  python3 modules/05_identification/identify_species.py \
  results/my_analysis/consensus/consensus_sequences.fasta \
  results/my_analysis/blast/
```

---

## Part 2 & 3: Answer Questions Interactively (40 points)

**After completing Steps 1-6 above, run the interactive question script:**

```bash
python3 answer_assignment.py
```

**What this script does:**
- Asks you questions about your analysis results
- Guides you through the HTML reports (QC, BLAST, phylogeny)
- Collects your answers in a structured format
- Saves to `answers.json` for automatic grading

**Why use this script?**
- ✓ No formatting errors (it creates perfect JSON)
- ✓ Interactive and easy to use
- ✓ Automatically graded when you push to GitHub
- ✓ Immediate feedback on correctness

**You will answer questions about:**
1. Species identification (BLAST results)
2. Quality control statistics
3. Phylogenetic tree interpretation
4. Mosquito diversity assessment

---

## (ALTERNATIVE) Part 2: Results Table (20 points)

**NOTE:** If you prefer to fill in answers manually instead of using the interactive script, you can fill in the tables and questions below. However, the interactive script (`python3 answer_assignment.py`) is recommended!

---

Fill in this table with the BLAST results for the class dataset:

| Sample | Species Identified | % Identity | Common Name |
|--------|-------------------|------------|-------------|
|        |                   |            |             |
|        |                   |            |             |
|        |                   |            |             |
|        |                   |            |             |

**Instructions:**
- Only include samples with consensus sequences (had both F and R pass QC)
- Use the **top BLAST hit** from each sample
- Species names must be in *italics*: *Genus species*
- Get % identity from BLAST HTML report

---

## Part 3: Analysis Questions (20 points)

### Question 1: Quality Control (5 points)

**a)** How many of the 30 class sequences (.ab1 files) passed QC?
**b)** How many samples had BOTH forward AND reverse reads pass?
**c)** Why is it important to have both F and R reads?

**Your answer:**
```
a)
b)
c)




```

### Question 2: Phylogenetic Tree (7 points)

Look at the phylogenetic tree (`results/my_analysis/phylogeny/tree.png`).

**a)** Do the class samples cluster together, or are they spread across different parts of the tree?
**b)** Which reference species are the class samples most closely related to?
**c)** What does this tell you about mosquito diversity in the class sampling locations?

**Your answer:**
```
a)
b)
c)




```

### Question 3: Species Identification (8 points)

**a)** What mosquito species did the class identify? List all unique species found.
**b)** Do the BLAST results (% identity) agree with where samples clustered on the tree?
**c)** Are these species known to occur in Southern California? (You may need to Google this!)
**d)** How confident are you in the species identifications? (Consider % identity scores)

**Your answer:**
```
a)
b)
c)
d)




```

---

## Part 4: Submission (GitHub)

**Submit your work by committing and pushing to GitHub:**

```bash
# Add your completed assignment and results
git add answers.json results/

# Commit your work
git commit -m "Complete DNA barcoding analysis and assignment"

# Push to GitHub
git push origin main
```

**Important:** Make sure to add `answers.json` (created by `python3 answer_assignment.py`)

**Verify your submission:**

1. Go to your GitHub repository
2. Click the **"Actions"** tab
3. Check that **"Auto-Grading"** workflow passed ✅
4. If it failed ❌, read the error message and fix missing files

**What auto-grading checks:**
- ✅ Tutorial completed (`results/tutorial/` has all reports)
- ✅ Analysis completed (`results/my_analysis/` has all reports)
- ✅ Assignment file exists (`assignment.md`)
- ✅ Answers are correct (since everyone has the same data, answers should match!)

**Auto-grading checks your:**
- Species identifications (BLAST results table)
- QC statistics (number of sequences that passed)
- Written answers (keyword matching for concepts)

---

## 📚 Grading Rubric

| Component | Points | Criteria |
|-----------|--------|----------|
| Part 1: Commands executed correctly | 60 | All 6 steps completed, files generated |
| Part 2: Results table | 20 | Accurate data, correct format |
| Part 3: Question 1 | 5 | Complete, accurate answers |
| Part 3: Question 2 | 7 | Thoughtful tree interpretation |
| Part 3: Question 3 | 8 | Complete species analysis |
| **Total** | **100** | |

---

## 🆘 Need Help?

1. **Re-run the tutorial**: `./tutorial-cs.sh` (Codespaces) or `./tutorial.sh` (Docker)
2. **Read the visual workflow**: `docs/pipeline_workflow.md`
3. **Understand IQ-TREE**: `docs/iqtree_guide.md`
4. **Check your work**: Compare your results to the tutorial results
5. **Codespaces issues**: See `.devcontainer/README.md`
6. **Ask your instructor**

---

## 🎓 Learning Goals

This assignment teaches you:
- How DNA barcoding identifies species
- Why quality control matters in sequencing
- How consensus sequences improve accuracy
- How to interpret phylogenetic trees
- How to use bioinformatics tools for real research

**Remember:** The goal is to LEARN the workflow, not just get it done. Take your time, look at the results, and think about what they mean!

Good luck! 🧬🦟
