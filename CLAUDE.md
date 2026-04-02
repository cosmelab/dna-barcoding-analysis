# DNA Barcoding Analysis

## ⚠️ WEB SEARCH: ONE TOOL PER QUERY — NO DUPLICATES ⚠️

**Pick ONE search tool per query. NEVER call both for the same question.**

| Need | Tool | Notes |
|------|------|-------|
| Quick answer | Perplexity `web_search` model=**sonar** | DEFAULT. $1/1K. **NEVER use sonar-pro unless user explicitly asks.** |
| Find paper URLs | Exa `research_papers` | Returns links + highlights |
| Find code examples | Exa `code_search` | GitHub/SO/docs |
| Fetch URL content | crawl4ai `fetch_page` | Free, saves to RAG. Exa `crawl` as fallback only. |
| Deep literature review | Perplexity `research_search` model=sonar | Upgrade to sonar-pro ONLY if sonar results are insufficient |

**RULES:**
- **sonar-pro costs 15x more than sonar.** Do NOT use it for general web searches.
- **Do NOT run Perplexity AND Exa for the same query.** That wastes money and tokens.
- **crawl4ai first** for fetching URLs (free). Exa crawl only if crawl4ai fails.

---

## MANDATORY - READ FIRST

```
~/Projects/tracking-hub/tracking/rules.yaml   # AGENT INSTRUCTIONS, credentials, common errors
tracking/state.yaml                            # THIS repo's current state
```

After EVERY prompt: Update tracking/state.yaml with what you did and if it worked.

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

## SESSION TRACKING (MANDATORY)

**Always use `/start` with explicit project argument. Never rely on auto-detect.**

```bash
/start dna-barcoding session-name
```

**Workflow:**
1. Claude starts -> SessionStart hook captures session ID automatically
2. `/start dna-barcoding <name>` -> locks project
3. Work -> `/wrapup` at end -> stores handoff
4. Resume: `claude --resume <session-id>` (shown by /pickup)

**Task ID prefix:** `DBA`

**Use MCP tools for session management:**
- `mcp__session-tracker__create_task` for task tracking
- `mcp__session-tracker__log_event` for documenting findings
- NEVER write YAML session files manually

## Rules

Follow rules in ~/Projects/tracking-hub/CLAUDE.md

## USE RAG TO FIND CONTEXT

Before starting any task, search RAG for relevant history:
```bash
python3 ~/Projects/tracking-hub/scripts/51_search_rag.py "your query"
python3 ~/Projects/tracking-hub/scripts/51_search_rag.py "topic" --project dna-barcoding-analysis
```

## SEND EMAILS (Cosme Lab Template v21)

Read template instructions first:
```bash
cat ~/Projects/tracking-hub/templates/email/agent-instructions.md
```

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

### Authentication Commands (MEMORIZE THESE)

**GitHub CLI:**
```bash
gh auth login
# Opens browser, follow prompts
```

**Docker (GHCR container registry):**
```bash
docker login
# Uses saved credentials (Username: cosmelab)
```

**Google Cloud:**
```bash
gcloud auth application-default login
```

NEVER expose tokens - use these simple commands.
