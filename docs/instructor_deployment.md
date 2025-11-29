# Instructor Deployment Guide

Complete setup instructions for deploying the DNA Barcoding Analysis assignment via GitHub Classroom with Codespaces support.

---

## Pre-Deployment Checklist

- [ ] Template repository is public (or accessible to GitHub Classroom)
- [ ] Docker image `cosmelab/dna-barcoding-analysis:latest` is published
- [ ] `.devcontainer/devcontainer.json` is configured
- [ ] Test data in `data/test_data/` works with tutorial
- [ ] Class sequence data in `data/student_sequences/` is ready
- [ ] Answer key values are finalized

---

## Step 1: Set Up the Answer Key Secret

The autograding workflow uses a secret to store correct answers (hidden from students).

### 1.1 Prepare the Answer Key

1. Open `.github/grading/answer_key.json`
2. Update values based on your class data:
   - `q1a_sequences_passed`: Number of sequences passing QC
   - `q1b_pairs_passed`: Number of F/R pairs
   - `num_samples`: Number of unique samples
   - `species_table`: Expected species identifications
   - `q2a_clustering`: Expected phylogeny interpretation
   - `q3a_species_list`: List of species found
   - `q3b_confidence`: Expected confidence level

### 1.2 Add Secret to Repository

1. Go to your **template repository** on GitHub
2. Navigate to: **Settings** > **Secrets and variables** > **Actions**
3. Click **"New repository secret"**
4. Name: `ANSWER_KEY_JSON`
5. Value: Paste the entire JSON content (without the `_comment` and `_instructions` fields)
6. Click **"Add secret"**

Example secret value:
```json
{
  "q1a_sequences_passed": "12",
  "q1b_pairs_passed": "4",
  "num_samples": "4",
  "species_table": [
    {"sample": "AT-HV1", "species": "Aedes albopictus", "percent_identity": "99.55"},
    {"sample": "AT-HV3", "species": "Culex pipiens", "percent_identity": "98.12"},
    {"sample": "AT-JM2", "species": "Culex pipiens", "percent_identity": "99.25"},
    {"sample": "AT-WL2", "species": "Culex pipiens", "percent_identity": "98.67"}
  ],
  "q2a_clustering": "Spread across different parts (multiple species)",
  "q3a_species_list": ["Aedes albopictus", "Culex pipiens"],
  "q3b_confidence": "Confident (97-99% identity)"
}
```

> **Note:** Secrets are NOT copied when GitHub Classroom creates student repos. You must add the secret to the **Classroom organization** or use a **grading workflow** that fetches from a secure location.

---

## Step 2: Enable Codespaces for Your Classroom

### 2.1 Organization Requirements

- Your organization must use **GitHub Team** (free for education)
- You must be an **organization owner** and **classroom admin**
- Apply for [GitHub Education benefits](https://education.github.com/discount_requests/application) if not already verified

### 2.2 Enable Codespaces

**Option A: When Creating a New Classroom**
1. Go to [classroom.github.com](https://classroom.github.com)
2. Click **"New classroom"**
3. Select your organization
4. Under "Codespaces in your Classroom," click **"Enable"**
5. Complete classroom creation

**Option B: For an Existing Classroom**
1. Go to [classroom.github.com](https://classroom.github.com)
2. Select your classroom
3. Click **"Settings"**
4. Under "GitHub Codespaces," click **"Enable"**

> **Important:** Enabling Codespaces applies to ALL repositories in the organization.

---

## Step 3: Create the Assignment

### 3.1 Assignment Settings

1. Go to your classroom on GitHub Classroom
2. Click **"New assignment"**
3. Configure:
   - **Title:** DNA Barcoding Analysis (or your preferred name)
   - **Deadline:** Set your due date (configured in workflow as Dec 5, 2025 07:59 UTC)
   - **Repository visibility:** Private (recommended)
   - **Grant admin access:** No (students don't need admin)

### 3.2 Starter Code

1. Select **"Add a template repository"**
2. Choose your template repo: `cosmelab/dna-barcoding-analysis`
3. **Important:** The repository is snapshotted at this moment

### 3.3 Enable Codespaces

1. Under "Add a supported editor," select **"GitHub Codespaces"**
2. This adds the "Open in GitHub Codespaces" button for students

### 3.4 Autograding (Optional)

GitHub Classroom's built-in autograding is separate from our workflow. Our custom autograding runs via GitHub Actions when students push to `results/`.

You can leave Classroom's autograding disabled since we use `.github/workflows/autograding.yml`.

### 3.5 Create Assignment

Click **"Create assignment"** and share the invitation link with students.

---

## Step 4: Student Workflow

Students will:

1. Accept the assignment via the invitation link
2. Click **"Open in GitHub Codespaces"** (or clone locally with Docker)
3. Run `./tutorial-cs.sh` (Codespaces) or `./tutorial.sh` (local Docker)
4. Run `./run-analysis-cs.sh` or `./run-analysis.sh`
5. Answer questions: `python3 answer_assignment.py`
6. Commit and push results

The autograding workflow triggers automatically when they push to `results/`.

---

## Updating After Assignment Creation

**GitHub Classroom takes a snapshot** of your template repository when you create an assignment. Changes to the template after this point do NOT affect existing student repositories.

### If You Need to Update:

**Option 1: Minor Fixes (Manual)**
- Students can pull changes manually:
  ```bash
  git remote add upstream https://github.com/cosmelab/dna-barcoding-analysis.git
  git fetch upstream
  git merge upstream/main --allow-unrelated-histories
  ```
- Communicate changes via announcement

**Option 2: Critical Fixes (New Assignment)**
1. Make fixes to your template repository
2. Create a **new assignment** in GitHub Classroom
3. Have students accept the new assignment
4. Students can copy their `results/` folder to the new repo

**Option 3: Script Distribution**
- For script-only updates, provide a download link
- Students can replace specific files without full re-clone

### What CAN Be Updated Without New Assignment:

- Docker image (`cosmelab/dna-barcoding-analysis:latest`) - pulled fresh each run
- Answer key secret - update in organization/repo settings
- Autograding workflow deadline - students pull workflow changes with git

---

## Grading and Results

### View Student Progress

1. Go to your classroom on GitHub Classroom
2. Click on the assignment
3. View submission status for each student

### Access Detailed Grading

1. Go to student's repository
2. Click **"Actions"** tab
3. Find the latest "Auto-Grading" workflow run
4. Download the **"grading-details"** artifact for detailed results

### Grading Artifacts Include:

- Question-by-question breakdown
- Partial credit calculations
- Late submission flags
- Detailed error messages (instructor only)

---

## Troubleshooting

### Codespaces Won't Start

- Check that Docker image is accessible: `docker pull cosmelab/dna-barcoding-analysis:latest`
- Verify `.devcontainer/devcontainer.json` syntax
- Check organization Codespaces quota

### Autograding Not Running

- Workflow triggers on pushes to `results/**`
- Check Actions tab for workflow errors
- Verify `ANSWER_KEY_JSON` secret is set

### Students Can't See Results

- Detailed grading results are in artifacts (instructor access)
- Students see only pass/fail summary in workflow output
- This prevents answer leakage

### Docker Image Issues

- Rebuild: `./scripts/build-container.sh`
- Push: `./scripts/push-container.sh`
- Students will get new image on next Codespace creation

---

## Timeline Checklist

### 1 Week Before Release
- [ ] Finalize class sequence data
- [ ] Test full pipeline end-to-end
- [ ] Set answer key secret
- [ ] Verify Docker image is current

### Day of Release
- [ ] Create assignment in GitHub Classroom
- [ ] Test student workflow with a test account
- [ ] Send invitation link to students

### During Assignment Period
- [ ] Monitor Actions tab for common errors
- [ ] Respond to student issues promptly
- [ ] Do NOT modify template (breaks student repos)

### After Deadline
- [ ] Download grading artifacts
- [ ] Review late submission flags
- [ ] Export final grades

---

## Resources

- [GitHub Classroom Documentation](https://docs.github.com/en/education/manage-coursework-with-github-classroom)
- [Codespaces with Classroom](https://docs.github.com/en/education/manage-coursework-with-github-classroom/integrate-github-classroom-with-an-ide/using-github-codespaces-with-github-classroom)
- [GitHub Education Benefits](https://education.github.com/benefits)
