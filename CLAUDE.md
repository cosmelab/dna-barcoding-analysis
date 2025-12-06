# DNA Barcoding Analysis

## MANDATORY STARTUP CHECKLIST

Read these files IN ORDER:

| # | File | Contains |
|---|------|----------|
| 1 | `tracking/state.yaml` | This project's state, container build instructions |
| 2 | `~/Projects/tracking-hub/tracking/state.yaml` | ALL APIs, credentials, security |

### In hub state.yaml, find:
- `all_apis_available.github` → gh CLI commands
- `credentials.security` → Security hooks (don't expose credentials)

### Container build (no GUI needed):
```bash
gh workflow run docker-build.yml
gh run watch
```

DO NOT run `ls`. DO NOT guess. READ THE TRACKING FILES.

---

## Hub Connection

This repo is tracked by the central **tracking-hub** at:
```
~/Projects/tracking-hub/
```

## Rules

Follow rules in ~/Projects/tracking-hub/CLAUDE.md

## HOW TO USE APIs

### Step 0: Credentials
If any API fails with 401/auth error, ask user: "Please load credentials"

**Perplexity (web search - USE INSTEAD OF HALLUCINATING):**
```bash
python3 ~/Projects/tracking-hub/tools/perplexity/perplexity.py "Your question" --model sonar-pro
```
Models: `sonar` (fast/cheap), `sonar-pro` (complex queries)

**NCBI/CrossRef (verify citations):**
```bash
cd /Users/lucianocosme/Projects/entm201l-fall2025/tools/reference-verification
python3 verify_citation.py "Author" "Year" "keywords"
python3 verify_citation.py --doi "10.xxxx/xxxxx"
```

**Google APIs (Gmail, Sheets, Calendar) - READY:**
- Auth: `gcloud auth application-default login` (already done)
- Credentials: `~/.config/gcloud/application_default_credentials.json`
- Python: `from googleapiclient.discovery import build`
- APIs enabled: Gmail, Sheets, Calendar, Drive, Forms

**Workflow:**
1. Load credentials (Step 0)
2. Perplexity (fast) → Get literature, DOIs
3. NCBI/CrossRef (verify) → Confirm citations are real

Full API docs: `~/Projects/tracking-hub/tracking/state.yaml` → `all_apis_available`

**HPC Note:** On HPC, use Singularity container. Paths inside container may differ.

### GitHub Authentication

Always use `gh auth login --web` - NEVER use local tokens
