# START HERE - DNA Barcoding Analysis

## For Students: 3 Simple Steps

### STEP 0: Setup (one time only)

1. **Make sure Docker Desktop is running**

2. **Login to Docker Hub:**
   ```bash
   docker login
   ```
   Enter your Docker Hub username and password.

3. **Pull the container:**
   ```bash
   docker pull cosmelab/dna-barcoding-analysis:latest
   ```

4. **Get your assignment repository:**

   **GitHub Classroom students - IMPORTANT:**

   a. **Click the assignment link from Canvas**

   b. **Accept the assignment** - this creates YOUR personal repo

   c. **Refresh the page** to see your new repository

   d. **Copy the clone URL** from the green Code button

   e. **Clone your repo:**
   ```bash
   git clone YOUR-REPO-URL
   cd dna-barcoding-analysis-YOUR-GITHUB-USERNAME
   ```

   Git will automatically prompt you to login to GitHub.

   **Or if using the template directly (not GitHub Classroom):**
   ```bash
   git clone https://github.com/cosmelab/dna-barcoding-analysis.git
   cd dna-barcoding-analysis
   ```

---

### OPTIONAL: Interactive Terminal (Advanced)

**Want a fancy terminal with colorful output?** The container includes a beautiful Zsh shell setup!

Instead of running `./tutorial.sh` directly, you can run commands inside the container interactively:

```bash
docker run --rm -it -v $(pwd):/workspace -w /workspace \
  cosmelab/dna-barcoding-analysis:latest zsh
```

**Inside the container, try these commands:**
- `ls` - Colorful file listings with icons
- `ll` - Detailed view with file sizes, dates, and git status
- `lt` - Tree view of directories
- `./tutorial.sh` - Run the tutorial
- `exit` - Leave the container

**What makes it fancy?**
- **Zsh** - Modern shell with smart tab completion
- **Dracula theme** - Professional dark color scheme
- **Git integration** - See git status in your prompt
- **Colorful output** - File types shown with different colors and icons

**For beginners:** Just use `./tutorial.sh` and `./run-analysis.sh` normally - they work perfectly without the fancy terminal!

---

### STEP 1: Learn the Workflow (REQUIRED - 15 minutes)

**Run the tutorial with test data:**
```bash
./tutorial.sh
```

This teaches you all 5 steps using example data. **DO NOT SKIP THIS!**

---

### STEP 2: Analyze the Class Mosquito Sequences

**Important:** Everyone analyzes the **same class dataset**. The .ab1 files are already in `data/student_sequences/`

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

### STEP 3: Answer Assignment Questions

**Run the interactive question script:**
```bash
python3 answer_assignment.py
```

This script will:
- Guide you through your analysis results
- Ask questions about QC, BLAST, and phylogeny
- Save your answers to `answers.json`

It's interactive and easy - just follow the prompts!

---

### STEP 4: Submit to GitHub

**Commit and push your work:**
```bash
git add answers.json results/
git commit -m "Complete DNA barcoding analysis"
git push origin main
```

**Check auto-grading:**
1. Go to your GitHub repository
2. Click the "Actions" tab
3. Look for ✅ (passed) or ❌ (failed)
4. If failed, click on it to see feedback

---

## That's It!

**The complete workflow:**
```bash
./tutorial.sh              # STEP 1: Learn (15 min)
./run-analysis.sh          # STEP 2: Analyze class data (5 min)
python3 answer_assignment.py  # STEP 3: Answer questions (10 min)
git add answers.json results/
git commit -m "Complete assignment"
git push origin main       # STEP 4: Submit (auto-graded!)
```

**Need help?**
- Re-run tutorial: `./tutorial.sh`
- Read: `docs/pipeline_workflow.md`
- Ask your TA
