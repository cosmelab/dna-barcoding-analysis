# Instructor Setup: Automated Grading

## How to Set Up the Answer Key (CRITICAL!)

The answer key is stored as a **GitHub Secret** so students can't see it.

### Step 1: Generate the Answer Key

Run the analysis yourself to get the correct answers:

```bash
# Run the analysis on the class data
./run-analysis.sh

# Answer the questions interactively
python3 answer_assignment.py

# This creates answers.json with the CORRECT answers
```

### Step 2: Copy the Answer Key JSON

```bash
# Copy the contents of answers.json
cat answers.json
```

**Copy the entire JSON output** (you'll need this in the next step)

### Step 3: Add to GitHub Secrets

**For the Template Repository:**

1. Go to: https://github.com/cosmelab/dna-barcoding-analysis/settings/secrets/actions
2. Click **"New repository secret"**
3. Name: `ANSWER_KEY_JSON`
4. Value: **Paste the entire JSON from answers.json**
5. Click **"Add secret"**

**For GitHub Classroom (ALL student repos):**

You need to add the secret to EACH student repository, OR use GitHub Classroom's secret sync feature:

1. Go to your GitHub Classroom assignment settings
2. Under "Secrets", add the secret there
3. It will sync to all student repositories

**Alternative: Manual per-repo setup**

If you don't use Classroom secrets:

```bash
# Use GitHub CLI to add secret to all student repos
gh secret set ANSWER_KEY_JSON --body "$(cat answers.json)" -R cosmelab/dna-barcoding-analysis-STUDENT1
gh secret set ANSWER_KEY_JSON --body "$(cat answers.json)" -R cosmelab/dna-barcoding-analysis-STUDENT2
# ... repeat for all students
```

Or write a script:

```bash
#!/bin/bash
ANSWER_KEY=$(cat answers.json)

for repo in $(gh repo list cosmelab --json name -q '.[] | select(.name | startswith("dna-barcoding-analysis-")) | .name'); do
  echo "Adding secret to $repo..."
  gh secret set ANSWER_KEY_JSON --body "$ANSWER_KEY" -R cosmelab/$repo
done
```

### Step 4: Verify It Works

1. Make a test commit to a student repo
2. Check the Actions tab
3. The grading step should run and check answers

## How the Grading Works

1. Student runs `python3 answer_assignment.py` → creates `answers.json`
2. Student commits and pushes `answers.json`
3. GitHub Actions workflow runs:
   - Checks tutorial complete
   - Checks analysis complete
   - **Grades answers.json against secret answer key**
   - Returns feedback and score
4. Student sees ✅ or ❌ in Actions tab

## Updating the Answer Key

If you need to update the answer key (e.g., fix a typo):

1. Update your local `answers.json`
2. Copy the new JSON
3. Go to repository secrets
4. Edit `ANSWER_KEY_JSON`
5. Paste the new JSON
6. Save

The next time students push, they'll be graded against the new key.

## Security

- ✅ Answer key stored as GitHub Secret (encrypted)
- ✅ Students cannot access secrets
- ✅ Answer key not in git history
- ✅ Answer key created in /tmp and deleted after grading
- ✅ GitHub Classroom can sync secrets to all student repos

## Troubleshooting

**Error: "ANSWER_KEY_JSON secret not found"**
- Make sure you added the secret to the repository
- Check the secret name is exactly `ANSWER_KEY_JSON`
- For Classroom, make sure secret sync is enabled

**Students getting wrong grades**
- Check your answer key JSON is valid
- Run the analysis yourself to verify correct answers
- Look at the Actions log for grading feedback

**Need to disable auto-grading temporarily?**
- Remove the secret (students will get "answers_graded=false")
- Or comment out the grading step in `.github/workflows/autograding.yml`
