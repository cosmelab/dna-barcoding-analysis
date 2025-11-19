# Project Activity Log

## Template
```
## YYYY-MM-DD – [agent] – [task-id]
actions:
  - [what was done]
results:
  - [what happened]
next:
  - [what's next]
```

---

## 2025-11-17 – migration – system
actions:
  - created new hybrid YAML tracking system for dna-barcoding-analysis repo
  - created state.yaml, decisions.yaml, tasks.yaml, log.md
results:
  - tracking system ready for project start
  - 5 initial setup tasks defined
  - 3 initial decisions documented (tool choices, dual pathway)
next:
  - create Dockerfile with current 2025 tools
  - write GUI and CLI pathway instructions

## 2025-11-17 – claude-sonnet – verify-tool-versions
actions:
  - searched for current best GUI tools for Sanger sequencing (2025)
  - searched for current best CLI tools for sequence analysis (2025)
  - verified Tracy actively maintained (updated Mar 2025)
  - confirmed MAFFT, IQ-TREE2, ClipKIT are current standards
  - documented findings in decisions.yaml (D005)
results:
  - all recommended tools verified as current and actively maintained
  - Tracy: https://github.com/gear-genomics/tracy (updated Mar 2025)
  - Alternative GUI tools identified: Jalview, AliView
  - Base image decided: Ubuntu 22.04 LTS (D003)
  - Perplexity API available for literature review (D004)
next:
  - create Dockerfile with verified tool versions
  - write GUI pathway instructions (Tracy web, Jalview/AliView)
  - write CLI pathway instructions (Docker usage)
  - create example workflow with sample data
  - set up GitHub Classroom assignment

## 2025-11-18 – claude-sonnet – strategy-revision
actions:
  - reviewed albopictus-diapause-rnaseq Dockerfile and terminal setup
  - analyzed Hoque et al 2022 paper on mosquito COI barcoding
  - extracted GenBank accessions from Tables 2-3 (85 sequences, 19 species)
  - reviewed GitHub Actions workflow for container builds
  - brainstormed student-friendly automation approach
  - updated tracking system with new strategy
results:
  - MAJOR PIVOT: Use micromamba base (not Ubuntu) based on RNA-seq success
  - Container strategy: <800MB, minimal tools, beautiful terminal (zsh/Dracula/colorls)
  - Reference data: Hoque et al 2022 - 85 mosquito COI sequences with improved AUCOS primers
  - Output format: HTML dashboard (not PDF) with interactive plotly visualizations
  - Workflow: CLI-only but ultra-simple (one command: analyze-sequences)
  - Decisions added: D006-D011 (container, data, visualization, UX strategy)
  - Tasks restructured: 12 new tasks covering container, scripts, automation, testing
next:
  - create minimal Dockerfile based on RNA-seq template
  - create Python script to fetch Hoque sequences from GenBank
  - write automation scripts (QC, alignment viz, BLAST, master pipeline)
  - create HTML dashboard generator
  - set up GitHub Actions for container builds
  - create GitHub Classroom template
  - test with real chromatograms when user returns
