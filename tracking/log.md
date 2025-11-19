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

## 2025-11-18 – claude-sonnet – implementation-sprint
actions:
  - created minimal Dockerfile (168 lines, 50% reduction from RNA-seq)
  - created fetch_reference_sequences.py and fetched 52/53 sequences from GenBank
  - created qc_chromatograms.py (BioPython .ab1 parser, HTML reports)
  - created identify_species.py (BLAST via BioPython, HTML reports)
  - created master_pipeline.py (QC -> BLAST workflow)
  - created GitHub Actions workflow (.github/workflows/docker-build.yml)
  - fixed ALL CAPS filename violations (renamed 3 files to lowercase)
  - committed and pushed to trigger container build
results:
  - Core functionality complete: QC + Species ID pipeline ready
  - Container build triggered on GitHub Actions (building now)
  - 6 major tasks completed in single session
  - Reference data prepared (52 Southern California mosquito COI sequences)
  - Modular structure created (modules/01_quality_control, modules/04_identification)
  - Docker container will be available at ghcr.io/cosmelab/dna-barcoding-analysis:latest
next:
  - wait for container build to complete (check GitHub Actions)
  - optional: create alignment visualization script (modules/02_alignment)
  - optional: create phylogeny script with Bio.Phylo (modules/03_phylogeny)
  - create module README files for students
  - create GitHub Classroom template repository structure
  - test end-to-end with real chromatograms when user provides samples

## 2025-11-18 – claude-sonnet – testing-and-visualization
actions:
  - received real test data from UC genomics core (8 .ab1 files from Nov 13)
  - tested QC script with real chromatograms
  - created Module 0: F/R sequence assembly (merge_forward_reverse.py)
  - created Module 2: MAFFT alignment wrapper (align_sequences.py)
  - created Module 3: IQ-TREE phylogeny with Bio.Phylo (build_tree.py)
  - improved chromatogram visualization based on user feedback:
    - added sequence overlay on traces (color-coded by quality)
    - show middle region (bases 50-200) avoiding poor quality ends
    - collapsible HTML sections for better layout
    - full sequence display with horizontal scroll
    - quality legend: green (Q≥30), orange (Q≥20), red (Q<20)
  - updated tracking system with all completed tasks
results:
  - All 5 modules complete (00-04): Assembly, QC, Alignment, Phylogeny, Identification
  - Real data test results:
    - 8 samples analyzed
    - 4 passed QC: AT99F/R (93% high quality), AT_ROCK_F/R (83% high quality)
    - 4 failed QC: AT83F/R, AT94F/R (2-10% high quality)
    - AT_ROCK F+R successfully merged into consensus
  - Improved HTML layout: professional, clean, student-friendly
  - Chromatogram visualization shows actual base calls with quality scores
next:
  - create README files for each module (student documentation)
  - update master pipeline to integrate all 5 modules
  - create GitHub Classroom template repository
  - check GitHub Actions container build status
  - optional: create alignment variation visualization (large/small letters)

## 2025-11-18 – claude-sonnet – dockerfile-debugging
actions:
  - debugged repeated GitHub Actions build failures (6+ attempts)
  - investigated actual error logs: colorls gem native extension compilation failing
  - identified root cause: missing C/C++ compilers (cmake, make, gcc, gxx)
  - systematically compared RNA-seq Dockerfile vs DNA-barcoding Dockerfile
  - added missing dependencies:
    - compilers: cmake, make, gcc, gxx (conda)
    - system libs: software-properties-common, dirmngr, lsb-release, gnupg2 (apt)
  - removed unnecessary packages (datamash, openjdk - not needed for DNA barcoding)
  - verified tool availability in conda channels (mafft, iqtree, blast all available)
results:
  - Fixed Dockerfile now matches RNA-seq structure for critical dependencies
  - Commit ca42519 pushed with compiler fixes
  - Build #7 triggered on GitHub Actions
  - Learned: colorls requires native C extension compilation (clocale gem)
next:
  - monitor build #7 for success
  - if successful: test container locally, deploy agents for module READMEs
  - if failed: investigate missing Ruby development headers
  - update master pipeline
  - create GitHub Classroom template

## 2025-11-18 – claude-sonnet – finalization
actions:
  - BUILD SUCCESS! Container finally built on GitHub Actions
  - created 5 wrapper scripts for easy commands:
    - analyze-sequences (full pipeline)
    - qc-sequences, align-sequences, build-tree, blast-identify
  - updated Dockerfile to install wrapper scripts and copy modules
  - created docker-compose.yml for GitHub Classroom deployment
  - updated main README with correct docker-compose commands
  - created GitHub Classroom template:
    - .github/classroom/autograding.json (100 points, 5 auto-graded tests)
    - ASSIGNMENT.md (complete lab instructions for students)
    - .gitignore for results and temp files
  - deployed 4 agents to populate tutorial script directories:
    - scripts/quality_control/ - 7 files (parse_ab1, trim, filter, batch_qc)
    - 06_alignment/scripts/ - 5 files (MAFFT wrapper, stats, visualization, trimming)
    - scripts/phylogeny/ - 7 files (IQ-TREE, tree viz, rooting, distances)
    - scripts/utilities/ - 11 files (FASTA tools, stats, reverse-comp, converter, subsample)
results:
  - Container: ✓ Built and available at ghcr.io/cosmelab/dna-barcoding-analysis:latest
  - Modules: ✓ All 5 modules complete with READMEs (Assembly, QC, Alignment, Phylogeny, ID)
  - Master Pipeline: ✓ Integrated all modules with command-line flags
  - Wrapper Scripts: ✓ 5 simple commands for students
  - Docker Setup: ✓ docker-compose.yml for one-command deployment
  - GitHub Classroom: ✓ Autograding + assignment instructions ready
  - Tutorial Scripts: ✓ 38 files, 14,244 lines, heavily commented for learning
  - Documentation: ✓ README updated, ASSIGNMENT.md created, STATUS.md summary
  - Total commits: 7 (ca42519 → 91c8978)
next:
  - user to test container with: docker-compose up -d
  - user to create GitHub Classroom assignment from this repo
  - user to test full pipeline with real student data
  - optional: test with Hoque reference sequences
  - ready for Week 8 lab deployment!
