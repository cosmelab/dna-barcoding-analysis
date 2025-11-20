# START HERE - DNA Barcoding Analysis

## For Students: 3 Simple Steps

### STEP 0: Setup (one time only)

1. **Make sure Docker Desktop is running**

2. **Login to Docker:**
   ```bash
   docker login
   ```

3. **Pull the container:**
   ```bash
   docker pull cosmelab/dna-barcoding-analysis:latest
   ```

4. **Get this repository:**
   ```bash
   git clone https://github.com/cosmelab/dna-barcoding-analysis.git
   cd dna-barcoding-analysis
   ```

---

### STEP 1: Learn the Workflow (REQUIRED - 15 minutes)

**Run the tutorial with test data:**
```bash
./tutorial.sh
```

This teaches you all 5 steps using example data. **DO NOT SKIP THIS!**

---

### STEP 2: Analyze Your Own Mosquito Sequences

**Put your .ab1 files in the data folder:**
```bash
cp ~/Desktop/*.ab1 data/student_sequences/
```

**Run the complete analysis:**
```bash
./run-analysis.sh
```

This script runs all 5 steps for you automatically:
1. Quality Control
2. Consensus Sequences
3. Combine with References
4. Alignment & Tree
5. Species ID (BLAST)

Your results will be in: `results/my_analysis/`

---

### STEP 3: Fill Out Assignment

Open these files and answer the questions in `assignment.md`:

- `results/my_analysis/01_qc/qc_report.html` - How many sequences passed?
- `results/my_analysis/04_phylogeny/tree.png` - What species cluster together?
- `results/my_analysis/05_blast/identification_report.html` - What species did you find?

---

## That's It!

**The workflow in 3 commands:**
```bash
./tutorial.sh              # Learn (15 min)
./run-analysis.sh          # Analyze your data (5 min)
# Fill out assignment.md
```

**Need help?**
- Re-run tutorial: `./tutorial.sh`
- Read: `docs/pipeline_workflow.md`
- Ask your TA
