# START HERE - DNA Barcoding Analysis

## For Students: 3 Simple Steps

### STEP 0: Setup (one time only)

**Complete these steps IN ORDER:**

#### 1. Open VSCode

Open Visual Studio Code on your computer (Mac or Windows).

#### 2. Open Terminal in VSCode

- **Mac:** Press `Control + ~` (tilde) or go to Terminal → New Terminal
- **Windows:** Press `Ctrl + ~` or go to Terminal → New Terminal

You should see a terminal window at the bottom of VSCode.

#### 3. Navigate to Where You Want the Project

```bash
# Example: Go to your Documents folder
cd ~/Documents

# Or create a course folder first
mkdir -p ~/Documents/ENTM201L
cd ~/Documents/ENTM201L
```

#### 4. Check Docker is Running

**IMPORTANT:** Make sure Docker Desktop is running BEFORE continuing!

```bash
# Check Docker is installed and running
docker --version
```

**Expected output:** `Docker version 24.x.x` or similar

**If you get an error:**
- Open Docker Desktop application
- Wait for it to fully start (whale icon stops animating)
- Try `docker --version` again

#### 5. Login to Docker Hub

```bash
docker login
```

Enter your Docker Hub username and password when prompted.

**Expected output:** `Login Succeeded`

#### 6. Pull the Analysis Container

```bash
docker pull cosmelab/dna-barcoding-analysis:latest
```

This downloads the container (~2.5GB). Wait for it to complete.

#### 7. Get Your Assignment Repository

**GitHub Classroom students - IMPORTANT:**

a. **Click the assignment link from Canvas**

b. **Accept the assignment** - this creates YOUR personal repo

c. **Refresh the page** to see your new repository

d. **Copy the clone URL** from the green Code button

e. **Clone your repo in the terminal:**
   ```bash
   git clone YOUR-REPO-URL
   cd dna-barcoding-analysis-YOUR-GITHUB-USERNAME
   ```

   Replace `YOUR-REPO-URL` with the URL you copied!
   Replace `YOUR-GITHUB-USERNAME` with your actual GitHub username!

   Git will automatically prompt you to login to GitHub.

**Or if using the template directly (not GitHub Classroom):**
```bash
git clone https://github.com/cosmelab/dna-barcoding-analysis.git
cd dna-barcoding-analysis
```

#### 8. Verify You're in the Right Place

```bash
# Check you're in the project directory
ls
```

**You should see:** `tutorial.sh`, `run-analysis.sh`, `data/`, `modules/`, etc.

**Now you're ready to run the analysis!**

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

## ⚠️ IMPORTANT: Where to Run Commands

**Run all commands ON YOUR COMPUTER (Mac/Windows/Linux), NOT inside Docker!**

✅ **CORRECT:**
- Open terminal/VSCode on your Mac or Windows computer
- Navigate to your project folder
- Run `./tutorial.sh` from there
- The script automatically uses Docker for you

❌ **WRONG:**
- Don't try to run commands inside the Docker container
- Don't manually start Docker and work inside it
- The scripts handle Docker automatically!

**Think of Docker like a tool:** The scripts call Docker for you. You don't need to "go inside" Docker.

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

**How it works:**
1. Script tells you which HTML file to open (e.g., QC report, BLAST results, tree)
2. You open the file and look at the results
3. Script asks you questions about what you see
4. You type your answers
5. Script saves everything to `answers.json`

**This is interactive!** The script guides you step-by-step through your results. Just follow the prompts and answer based on what you see in the HTML reports.

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
- Ask your instructor
