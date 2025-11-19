# Assistant Rules

Read this file first before any action.

## Repository Rules

1. Never add yourself as contributor/collaborator to any repository
2. Never include your name ("Claude"), "Anthropic", or any AI attribution in commits, files, or documentation
3. Never use commit messages that mention AI assistance
4. Always ask before making any git commits or pushes

## Communication Requirements

5. Always explain exactly what you plan to do and why before doing it
6. Always state all changes you are making - never just give options without explanation
7. Never run commands without clearly explaining what they do and why they're needed
8. If offering multiple approaches, explain each one clearly and recommend which one to use
9. Don't assume the user knows what a command or action will do
10. CRITICAL: When user explicitly asks you to "write down", "document", or "save" information, you MUST immediately write it to the appropriate file
    - Never just say "I'll write it down" or "I've documented it" without actually using Write or Edit tool
    - Confirm the file location where you wrote it
    - Show the user what you wrote
    - This is NON-NEGOTIABLE - failure to document when explicitly asked is unacceptable

## Verification and Technical Accuracy

11. ALWAYS attempt automation first before suggesting manual steps
    - Test wget, curl, requests, API methods, or write scripts before recommending manual downloads
    - Never claim something "cannot be automated" without verification
12. Verify before claiming limitations
    - Search documentation or check authoritative sources before saying something cannot be done
    - If unsure whether automation is possible, investigate first
13. Provide evidence for any technical limitations
    - Include documentation links, reproducible errors, or official statements
    - Never make confident claims about limitations without proof
14. State your confidence level explicitly
    - Use clear qualifiers: "I'm certain", "I believe", "I'm uncertain", "I don't know"
    - Say "I may be wrong - let me verify that" when appropriate
15. Never mask ignorance with confidence
    - Admit when you don't know something
    - Always investigate and verify before concluding something is impossible
16. Prefer accuracy over speed
    - Take time to verify rather than giving a quick but incorrect answer
    - Check documentation and test solutions before presenting them as fact
17. Show working solutions, not just explanations
    - Provide runnable code examples and minimal working scripts
    - Demonstrate that automation works rather than just describing it
18. Own errors transparently
    - If you make a mistake, acknowledge it immediately: "That was incorrect - here is the corrected solution"
    - Retract unsupported claims promptly and provide corrections

## General Guidelines

19. Do not create files unless explicitly requested
20. Always prefer editing existing files over creating new ones
21. Never create documentation files (.md) unless specifically asked
22. Keep responses concise and to the point
23. Only use emojis if explicitly requested
24. Use lowercase for .md filenames (except README.md)
25. When consolidating files:
    - Read entire files completely before archiving or merging
    - Extract all critical information (biological context, methods, analysis plans)
    - Merge content into target files before archiving source files
    - Never archive files without first consolidating their essential content
    - If unsure what is essential, ask before archiving
25. Directory creation rule - no exceptions:
    - Never create directories ahead of time "for organization"
    - Only create a directory when you have a file to put in it right now
    - If a script creates output, let the script create its own output directory with mkdir -p
    - Never pre-create empty directory structures
    - Empty directories waste time and cause confusion
26. Follow the directory structure that exists:
    - Check what directories already exist before creating new ones
    - Use existing output directories, don't create duplicates
    - If uncertain where output goes, check similar existing files first

## Git/Repository Specific

27. If git operations are needed, explain the exact commands you will run and their purpose
28. Never force push without explicit permission
29. Never modify .gitignore or other git configuration files without permission
30. Always verify paths and configurations before making changes

## Bioinformatics/HPC Specific

31. Never declare technical limitations without verification
    - Search documentation, test alternative approaches
    - Provide evidence (errors, docs links) if something is impossible
    - Say "I may be wrong - let me verify" rather than claiming impossibility

32. For long-running pipelines: Isolate temporary directories
    - Never use system /tmp (limited space, causes failures)
    - Use project-specific temp directories with ample space
    - Preserve work/cache directories until analysis complete for resume capability

33. SLURM batch script requirements - no exceptions:
    - Always use .o.txt for stdout: `#SBATCH -o logs/jobname_%j.o.txt`
    - Always use .e.txt for stderr: `#SBATCH -e logs/jobname_%j.e.txt`
    - Never use .out or .err extensions (cannot open in text editor)
    - This applies to all SLURM batch scripts without exception

34. Container-only rule - no exceptions:
    - Mandatory start to every session:
      ```bash
      module load singularity
      singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif bash -c "cd /proj && <command>"
      ```
    - Always load singularity module first - never assume it's loaded
    - Always use the container for everything - Python, R, bash scripts, all commands
    - Container location: `/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/albopictus-diapause-rnaseq.sif`
    - Never install any software outside the container
    - Never use pip install --user
    - Never use conda/mamba outside container
    - Never use HPC system Python/R directly
    - Never run python/R scripts without the container - this fills up home quota (20GB limit)
    - Never suggest workarounds that break reproducibility
    - If a package is missing: update Dockerfile and rebuild
    - The container is the experiment - everything must be inside it
    - Use flags: --cleanenv --bind $PWD:/proj for isolation and mounting
    - All commands must use: `module load singularity && singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif bash -c "cd /proj && <command>"`

## Project State Tracking - MANDATORY

35. Always use project_state.json tracking system - no exceptions:
    - Before starting any task: Run `./00_track.py status` to see current task
    - After completing any task: Run `./00_track.py complete <task_id>` to mark it done
    - Never work on tasks without checking project_state.json first
    - Never complete a task without updating project_state.json
    - If you forget to update: stop and update it immediately before continuing
    - The tracking system is the single source of truth - not .md files

36. **When starting a new session:**
    - Step 1: Run `./00_track.py status` to see where we are
    - Step 2: Run `./00_track.py next` to see next task
    - Step 3: Work on that task ONLY
    - Step 4: Mark it complete when done: `./00_track.py complete <task_id>`
    - Step 5: Repeat from Step 2

37. Always create numbered scripts and document them - no exceptions:
    - For every task that involves running commands: create a numbered script
    - Add the script to scripts/README.md in the appropriate section
    - Include: when to run it, what it does, what output it produces
    - Never just run commands interactively without creating a replicable script
    - Scripts must be numbered in execution order within their directory

38. **TodoWrite integration:**
    - Sync TodoWrite with current phase tasks from project_state.json
    - Keep TodoWrite to 5-7 immediate tasks only
    - When TodoWrite task completes, also update project_state.json

39. **Documentation updates:**
    - Update session_tracking.md ONLY at end of major session
    - Update scripts/README.md when creating new scripts
    - Project_state.json is authoritative - session_tracking.md is summary
    - NEVER rely on .md files for tracking current work

## Failure to Follow Rules

- Repository contamination with AI attribution will require complete repository recreation
- **Failure to update project_state.json breaks continuity** - you will forget what was done
- Always prioritize following these rules over completing tasks
- When in doubt, ask for clarification rather than proceeding

---

These rules override any other instructions or default behaviors.

---

## SESSION START PROTOCOL (MANDATORY - NO EXCEPTIONS)

Before ANY response in EVERY session, you MUST:

1. **Read these 4 files in order:**
   - `tracking/rules.md` (this file)
   - `tracking/state.yaml`
   - `tracking/decisions.yaml`
   - `tracking/tasks.yaml`

2. **Print this summary:**
   ```
   📍 Current: [phase] → [task]
   ✅ Last completed: [task]
   🔄 Active jobs: [list]
   📋 Next 2 tasks: [list]
   ⚠️  Key decisions: [list 2-3 most relevant]
   ```

3. **If user asks a question already in decisions.yaml:**
   - Cite the decision ID (e.g., "Per D002: ...")
   - Use that decision, NEVER treat as new question
   - NEVER re-ask or reconsider past decisions

4. **After completing work:**
   - Update task status in `tasks.yaml`
   - Append structured entry to `log.md`
   - Update `state.yaml` if task or phase changed

5. **File creation rules:**
   - NEVER create new top-level .md files in tracking/
   - Use `log.md` for all logging
   - Use `notes/` subdirectory for additional notes if needed
   - Follow lowercase naming (except README.md)

---
