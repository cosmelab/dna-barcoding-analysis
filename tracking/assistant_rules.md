# Assistant Rules for DNA Barcoding Analysis Repo

## General Guidelines

1. Never add AI attribution in commits, files, or documentation
2. Always prefer editing existing files over creating new ones
   - NEVER create v2, v3, etc. versions - use git for version control
   - If replacing a file entirely: archive old one, document in tracking, then replace
   - Exception: Archive files go in archive/ directory with timestamps
3. Do not create files unless explicitly requested
4. Never create documentation files (.md) unless specifically asked
5. Keep responses concise and to the point
6. Only use emojis if explicitly requested

## Tracking System

7. Always update master_state.json when completing tasks
8. Document all changes in tracking before committing
9. Use visual_progress.py to generate progress reports
10. Archive completed work when appropriate

## Docker & Tools

11. Research current tools before recommending (don't rely on training data)
12. Test solutions before presenting as fact
13. Document tool versions and sources
14. Prefer open-source tools for reproducibility

## Git/Repository Specific

15. Explain git commands before running
16. Never force push without permission
17. Let git handle version control - don't create versioned files