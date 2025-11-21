# Project Development Tracking System

This directory contains project management and development files using a hybrid YAML tracking system. This repository was developed with agentic AI assistance following structured tracking protocols.

## Overview

This tracking system ensures:
- ✅ Continuity across sessions (AI agents can resume work)
- ✅ Design decisions are documented and never re-debated
- ✅ Clear project state and progress tracking
- ✅ Transparent development process (open science)

## File Structure

### Core Tracking Files

- **`state.yaml`** - Current project state (phase, task, status, blockers)
- **`decisions.yaml`** - All design decisions (D001-D018) with rationale
- **`tasks.yaml`** - Task definitions organized by phase
- **`log.md`** - Chronological work log (structured entries)

### AI Assistant Configuration

- **`rules.md`** - Comprehensive assistant rules with session protocol (39 rules)
- **`assistant_rules.md`** - Basic assistant rules (17 core rules)

**For AI agents:** Read `rules.md` for full protocols. The session start protocol (lines 170-206) ensures continuity across sessions.

### Documentation Files

- **`student_data_results.md`** - Example workflow results on class dataset
- **`README.md`** - This file (tracking system overview)

### Utility Scripts

- **`monitor.py`** - Progress tracking script
- **`visual_progress.py`** - Generate progress reports
- **`visual_progress.txt`** - Generated progress snapshot
- **`master_state.json`** - Legacy state tracking (deprecated, use state.yaml)

### Archive

- **`archive/`** - Archived old files (no longer active)

## How to Use This System with Agentic AI

### Session Start Protocol (for AI Agents)

Before any work in a session, AI agents MUST:

1. **Read these 4 files in order:**
   - `tracking/rules.md` (comprehensive rules)
   - `tracking/state.yaml` (current state)
   - `tracking/decisions.yaml` (past decisions)
   - `tracking/tasks.yaml` (task definitions)

2. **Print status summary:**
   ```
   📍 Current: [phase] → [task]
   ✅ Last completed: [task]
   🔄 Active jobs: [list]
   📋 Next tasks: [list]
   ⚠️  Key decisions: [list relevant]
   ```

3. **Reference past decisions:**
   - If user asks something already decided, cite decision ID (e.g., "Per D015: ...")
   - NEVER re-ask or reconsider documented decisions
   - Add new decisions to `decisions.yaml` with next ID

4. **Update after work:**
   - Mark tasks complete in `tasks.yaml`
   - Append structured entry to `log.md`
   - Update `state.yaml` if phase/task changed

### Key Principles

- **Single source of truth:** `state.yaml` + `decisions.yaml` (not .md files)
- **No rework:** Documented decisions are final unless explicitly revisited
- **Transparency:** All design choices visible and justified
- **Continuity:** Any agent can resume work by reading tracking files

## For Students: What's Useful Here?

### Educational Files (interesting to explore):

1. **`decisions.yaml`** - See why we chose these tools, this workflow, this structure
   - D001-D005: Tool selection (MAFFT, IQ-TREE, BLAST)
   - D006-D010: Workflow design (consensus sequences, paired reads)
   - D011-D014: Repository structure and Docker strategy
   - D015-D018: GitHub Classroom integration and public release

2. **`log.md`** - Development timeline showing how the project evolved

3. **`student_data_results.md`** - Example results from running the complete workflow
   - Shows what output to expect
   - Explains why some samples fail QC
   - Sample success rates and species identified

### Internal Files (can ignore):

- State tracking files (`state.yaml`, `tasks.yaml`, `master_state.json`)
- Scripts (`monitor.py`, `visual_progress.py`)
- AI configuration (`rules.md`, `assistant_rules.md`)

## Tracking System Structure

### state.yaml Schema

```yaml
project: project-name
current_phase: phase-name
current_task: task-id
last_updated: ISO-8601-timestamp
active_jobs: [list]
blocking_issues: [list]
completed_today: [list]
current_status: [list]
pending_decisions: [list]
next_up: [list]
```

### decisions.yaml Schema

```yaml
- id: D###
  date: YYYY-MM-DD
  question: "What needs to be decided?"
  decision: "Clear decision statement"
  rationale: "Why this decision?"
  status: active | superseded | deprecated
  superseded_by: D### (if applicable)
  tags: [list]
```

### log.md Format

```markdown
## YYYY-MM-DD – [agent] – [task-id]
actions:
  - [what was done]
results:
  - [what happened]
next:
  - [what's next]
```

## Why This Approach?

This tracking system enables:

1. **Reproducibility** - Every decision documented with rationale
2. **Efficiency** - AI agents don't re-debate past decisions
3. **Continuity** - Any session can resume from current state
4. **Transparency** - Open science principles in practice
5. **Learning** - Students see real software development process

## Questions?

### About the tracking system:
- See `rules.md` for full protocols
- See `decisions.yaml` for design rationale

### About the DNA barcoding workflow:
- See `README.md` in root directory
- See `docs/pipeline_workflow.md` for detailed steps
- Ask during lab or office hours
