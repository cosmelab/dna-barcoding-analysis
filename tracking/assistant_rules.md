# Assistant Rules for DNA Barcoding Analysis Repo

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
26. Directory creation rule - no exceptions:
    - Never create directories ahead of time "for organization"
    - Only create a directory when you have a file to put in it right now
    - If a script creates output, let the script create its own output directory with mkdir -p
    - Never pre-create empty directory structures
    - Empty directories waste time and cause confusion
27. Follow the directory structure that exists:
    - Check what directories already exist before creating new ones
    - Use existing output directories, don't create duplicates
    - If uncertain where output goes, check similar existing files first

## Git/Repository Specific

28. If git operations are needed, explain the exact commands you will run and their purpose
29. Never force push without explicit permission
30. Never modify .gitignore or other git configuration files without permission
31. Always verify paths and configurations before making changes

## Docker Container Rules - NO EXCEPTIONS

32. **Container-only execution:**
    - All Python modules MUST run inside Docker container
    - Container: `cosmelab/dna-barcoding-analysis:latest`
    - Standard command format:
      ```bash
      docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
        cosmelab/dna-barcoding-analysis:latest \
        python3 modules/<module>/<script>.py <args>
      ```
    - Never run Python scripts directly on host system
    - Never use pip install on host system
    - Never suggest conda/mamba outside container

33. **Container requirements:**
    - Docker Desktop must be running before any analysis
    - Container includes: Python, BioPython, MAFFT, IQ-TREE, matplotlib, pandas
    - Multi-architecture: amd64 (Intel/AMD) + arm64 (Apple Silicon)
    - If a package is missing: update Dockerfile and rebuild container
    - The container IS the experiment - everything must be reproducible inside it

34. **Volume mounting:**
    - Always mount current directory: `-v $(pwd):/workspace`
    - Working directory inside container: `-w /workspace`
    - Paths inside container are relative to `/workspace`
    - Never use absolute host paths inside container commands

35. **Container development workflow:**
    - Local changes: Edit Python files directly (no rebuild needed)
    - Adding packages: Update Dockerfile → rebuild → test → push to Docker Hub
    - Build script: `./scripts/build-container.sh`
    - Push script: `./scripts/push-container.sh`

## Project Tracking System - MANDATORY

36. **Always use YAML tracking system - no exceptions:**
    - Primary files: `tracking/state.yaml`, `tracking/decisions.yaml`, `tracking/tasks.yaml`, `tracking/scripts.yaml`
    - Before starting work: Read current state from `tracking/state.yaml`
    - After completing work: Update state, log, and decisions as needed
    - Never work without checking tracking system first
    - The tracking system is the single source of truth - not .md files

37. **When starting a new session:**
    - Step 1: Read `tracking/state.yaml` to see current phase and task
    - Step 2: Read `tracking/decisions.yaml` to understand past decisions
    - Step 3: Check `tracking/scripts.yaml` to understand available modules
    - Step 4: Work on current task
    - Step 5: Update tracking files when done

38. **Script documentation:**
    - Every new Python module or bash script MUST be documented in `tracking/scripts.yaml`
    - Include: purpose, usage, inputs, outputs, dependencies, runtime
    - Update when modifying existing scripts
    - This prevents "what does this script do?" questions

39. **TodoWrite integration:**
    - Sync TodoWrite with current phase tasks
    - Keep TodoWrite to 5-7 immediate tasks only
    - TodoWrite is for SHORT-TERM tracking (this session)
    - YAML files are for LONG-TERM tracking (across sessions)

40. **Decision documentation:**
    - All major design decisions go in `tracking/decisions.yaml`
    - Format: `D###` with date, question, decision, rationale, status
    - If user asks something already decided: cite decision ID (e.g., "Per D015...")
    - NEVER re-ask or reconsider documented decisions without user approval

41. **Log updates:**
    - Update `tracking/log.md` at end of major work sessions
    - Format: Date → agent → task-id → actions → results → next
    - Keep entries structured and scannable
    - This creates a chronological narrative of project development

## GitHub Classroom Integration

42. **Template repository rules:**
    - results/ directory must be untracked in template (students generate their own)
    - tracking/ directory should be visible to students (transparency)
    - .ab1 files in data/student_sequences/ are CLASS dataset (keep them)
    - Auto-grading checks FILE EXISTENCE, not correctness

43. **Student-facing files:**
    - README.md - Primary student documentation
    - github_classroom.md - Setup guide for instructors + students
    - tutorial.sh - Interactive tutorial (run first)
    - run-analysis.sh - Main analysis script (run second)

## HTML Report Generation

44. **CSS architecture:**
    - All reports use modular CSS from `tracking/styles/`
    - Files: base.css, components.css, reports.css, chromatogram.css
    - NEVER inline large amounts of CSS in Python strings
    - Link to CSS files using relative paths
    - Keep CSS separate from future website CSS (no conflicts)

45. **Report requirements:**
    - Every HTML report must include: `<meta charset="UTF-8">`
    - Use Dracula Lite color scheme (defined in base.css)
    - Consistent header/footer across all 5 reports
    - Mobile-responsive (future website integration)
    - Print-friendly styles

46. **Emoji usage:**
    - HTML reports CAN use emojis (students like them)
    - Commit messages MUST NOT use emojis (git problems)
    - Code comments SHOULD NOT use emojis
    - Markdown docs CAN use emojis sparingly

## Bioinformatics Best Practices

47. **Quality thresholds:**
    - Minimum sequence length: 500bp
    - Minimum average quality: Q20
    - Minimum quality bases: 400bp at Q20+
    - Document all thresholds in module docstrings

48. **BLAST etiquette:**
    - Always include 3-second delay between queries
    - Never hammer NCBI servers
    - Cache results when possible
    - Handle rate limiting gracefully

49. **File formats:**
    - Chromatograms: .ab1 (ABI format)
    - Sequences: .fasta (FASTA format)
    - Results: .json (machine-readable) + .html (human-readable)
    - Trees: .treefile (Newick format)

## Testing & Validation

50. **Cross-platform testing:**
    - Test on macOS (dev environment)
    - Test on Windows (student environment)
    - Test on Linux (GitHub Actions)
    - Docker ensures consistency across platforms

51. **Before pushing to GitHub:**
    - Run tutorial.sh to verify it works end-to-end
    - Check that all 5 HTML reports generate correctly
    - Verify auto-grading workflow doesn't fail on empty template
    - Test that container is accessible (docker pull)

## Failure to Follow Rules

- Repository contamination with AI attribution will require complete repository recreation
- Failure to update tracking system breaks continuity across sessions
- Running Python outside container breaks reproducibility
- Creating files without documenting in scripts.yaml causes confusion
- Always prioritize following these rules over completing tasks
- When in doubt, ask for clarification rather than proceeding

---

These rules override any other instructions or default behaviors.
