# DNA Barcoding Analysis - Architecture Design Document

**Version**: 1.0
**Date**: 2025-11-07
**Designer**: architecture-designer-agent
**Purpose**: Document all architectural decisions and rationale for the dna-barcoding-analysis repository

---

## Executive Summary

This repository provides a comprehensive, modular educational platform for DNA barcoding analysis. It combines three critical components:

1. **Command-line tutorials** - Linux, Python, R basics for beginners
2. **Analysis workflows** - Quality control through species identification
3. **Container infrastructure** - Reproducible computing environment

**Target Users**:
- ENTM201L students (UC Riverside) - first exposure to bioinformatics
- Lab members - learning computational biology
- Broader community - anyone doing DNA barcoding

**Key Design Principles**:
- Modularity (use any module independently)
- Beginner-friendly (assume zero command-line experience)
- Reproducibility (containerized environment)
- Scalability (GitHub Classroom ready)

---

## Repository Structure Rationale

### Numbered Module Approach (00-08)

**Decision**: Use numbered directories (00_introduction, 01_linux_basics, etc.)

**Rationale**:
- Clear learning progression
- Students know where to start (00)
- Logical workflow matches actual analysis pipeline
- Easier to reference ("go to Module 05")
- Self-documenting order

**Alternative Considered**: Semantic names only (introduction/, linux/, etc.)
**Rejected Because**: No clear starting point, harder to sequence learning

### Module Independence

**Decision**: Each module is self-contained with its own README, data, exercises

**Rationale**:
- Instructors can assign specific modules
- Students can focus on weak areas
- Advanced users skip basics
- Easier to maintain (changes isolated)
- Can extract modules for other courses

**Implementation**:
```
XX_module_name/
├── README.md       # Complete module documentation
├── data/           # Practice datasets
├── exercises/      # Hands-on problems
├── solutions/      # Answer keys
└── [module-specific content]
```

### Tutorials First, Analysis Second (Modules 01-03)

**Decision**: Teach computing skills before bioinformatics analysis

**Rationale**:
- Students need CLI skills for *all* bioinformatics
- Better to learn grep/awk with simple examples than fail on real data
- Reduces frustration (learn one thing at a time)
- Python/R skills transferable to other projects

**Learning Curve Management**:
- Module 00: No skills required (just reading)
- Module 01: Complete beginner command line
- Modules 02-03: Programming basics (build on Module 01)
- Modules 05-08: Apply skills to real problems

---

## Linux Tutorial Design (Module 01)

### Why So Comprehensive?

**Decision**: Create 8 detailed lessons + cheatsheet + 20+ exercises

**Rationale**:
1. **Students have zero CLI experience**
   - Cannot assume any prior knowledge
   - Need hand-holding initially
   - Build confidence gradually

2. **CLI is fundamental to all bioinformatics**
   - HPC clusters are Linux
   - Most tools are command-line only
   - Industry standard

3. **Pipes and streams are critical**
   - DNA barcoding = text file processing
   - grep, sed, awk are daily tools
   - Workflow automation requires this knowledge

### Lesson Structure

**Decision**: Each lesson follows this pattern:
- Introduction (why learn this?)
- Command syntax
- Practical examples
- Bioinformatics use cases
- Common errors and solutions
- Exercises

**Rationale**:
- Motivation first (students need "why")
- Show don't tell (examples > explanation)
- Anticipate mistakes (reduce frustration)
- Practice makes perfect (exercises essential)

### Streams and Redirection Emphasis

**Decision**: Dedicated 45-minute lesson on pipes (06_streams_redirection.md)

**Rationale**:
- **This is the power tool** of Linux for bioinformatics
- Students struggle with this concept
- Once mastered, productivity 10x
- Enables complex workflows

**Examples Focused on Biology**:
```bash
# Not just "cat file.txt | grep pattern"
# But actual use cases:
grep "^>" sequences.fasta | sed 's/>//' | sort > ids.txt
```

### Cheatsheet Strategy

**Decision**: Comprehensive cheatsheet (500+ lines)

**Rationale**:
- Reference during exercises
- Reduces googling (faster learning)
- Organized by task (navigation, searching, etc.)
- Bioinformatics-specific section
- Print and keep handy

---

## Python Tutorial Design (Module 02)

### Jupyter Notebooks vs Scripts

**Decision**: Use Jupyter notebooks for learning, scripts for automation

**Rationale**:
- **Notebooks for beginners**:
  - Interactive (see results immediately)
  - Mix code, text, plots
  - Can experiment safely
  - Great for teaching

- **Scripts for production**:
  - Reusable
  - Version control friendly
  - Automation ready
  - HPC compatible

**Structure**:
```
02_python_basics/
├── notebooks/        # Learning (interactive)
└── scripts/          # Production (reusable)
```

### BioPython Focus

**Decision**: Teach BioPython from the start, not generic Python

**Rationale**:
- Students want to process sequences NOW
- More motivating than "Hello World"
- Learn Python syntax in context
- Immediate applicability to their data

**Example Approach**:
- Lesson 1: Python basics (BUT using Bio.Seq objects)
- Lesson 2: BioPython explicitly
- Lesson 3: Real FASTA parsing

---

## R Tutorial Design (Module 03)

### ggtree Emphasis

**Decision**: Focus on ggtree for tree visualization, not base R plotting

**Rationale**:
- ggtree is publication-quality out of the box
- Consistent with modern R (tidyverse ecosystem)
- Easier for beginners than base graphics
- What professionals actually use

### Parallel to Python

**Decision**: Similar structure to Python module (scripts, exercises, data)

**Rationale**:
- Consistency reduces cognitive load
- Students know what to expect
- Choose Python OR R (not mandatory to do both)
- Both lead to same analysis capability

---

## Data Organization (Module 04)

### Three-Directory Structure

**Decision**:
```
04_data/
├── reference_sequences/    # Published data
├── student_sequences/      # User data
└── metadata/               # Documentation
```

**Rationale**:
1. **reference_sequences/**: Immutable, shared
   - Don't want students editing this
   - Source of truth for comparisons

2. **student_sequences/**: Mutable, personal
   - Students add their samples here
   - Instructor can review this directory
   - .gitkeep ensures directory exists

3. **metadata/**: Critical for reproducibility
   - Specimen info, GPS, dates
   - Template provided
   - Teaches data management

### Naming Conventions

**Decision**: Enforce `LastName_SpecimenID_COI.fasta` format

**Rationale**:
- No spaces (breaks scripts)
- No special characters (safer)
- Identifiable (know whose sample)
- Sortable (alphabetical by student)
- Searchable (grep for student name)

---

## Analysis Modules (05-08) Design

### Pipeline Order Matches Reality

**Decision**: Modules follow actual analysis workflow:
05 (QC) → 06 (Alignment) → 07 (Phylogeny) → 08 (ID)

**Rationale**:
- Matches how students actually work
- Output of module N is input to module N+1
- Can't skip steps (dependencies clear)
- Reinforces proper workflow

### Tool Selection Rationale

| Step | Tool | Why This Tool? |
|------|------|----------------|
| QC | BioPython | Standard for .ab1 parsing |
| Alignment | MAFFT | Fast, accurate, beginner-friendly |
| Phylogeny | IQ-TREE | Modern, automatic, best support |
| ID | BLAST + BOLD | Industry standards |

**Alternatives Provided**:
- Alignment: MUSCLE, ClustalW (for comparison)
- Phylogeny: RAxML, FastTree (learn trade-offs)

**Decision Rationale**:
- Primary tool = best practice
- Alternatives = educational value
- Students learn when to use each

### Scripts vs Containers

**Decision**: Provide both raw commands AND containerized workflows

**Rationale**:
- **Raw commands**: Understand what's happening
- **Container**: Reproducibility, ease of use
- Students learn both approaches
- Transition to HPC clusters later (need both)

---

## Container Architecture

### Base Image Selection

**Decision**: Use mambaorg/micromamba (not conda, not Ubuntu+pip)

**Rationale**:
- Faster than conda
- Lighter than full Anaconda
- Reliable package management
- Easy to add bioinformatics tools (bioconda)

**Borrowed from**: entm201l-fall2025/Dockerfile (proven in production)

### Multi-Stage Install Strategy

**Decision**: Install in order:
1. System dependencies (apt)
2. Core tools (micromamba)
3. Bioinformatics (bioconda)
4. Python packages (pip)
5. R packages (install.packages)
6. Shell configuration (zsh, starship)

**Rationale**:
- Layered for caching (faster rebuilds)
- System → language → tools
- Dependencies resolved before use
- Configuration last (changes often)

### Directory Structure Inside Container

**Decision**:
```
/workspace  # Student mounts their data here
/course     # (not used in this repo, kept for compatibility)
```

**Rationale**:
- `/workspace` is intuitive
- Students mount project root
- Can access all modules
- Consistent with VS Code Dev Containers

---

## GitHub Classroom Integration

### Template Repository Design

**Decision**: Repository is a template, students fork via Classroom

**Rationale**:
- Each student gets personal copy
- Can commit their work
- Instructor can review all repos
- Students learn Git workflow
- Grade via GitHub (Classroom tools)

### What Students Commit

**Design**:
```
04_data/student_sequences/      # Their sequences (COMMIT)
05_quality_control/cleaned/     # Cleaned data (COMMIT)
06_alignment/                   # Alignments (COMMIT)
07_phylogeny/                   # Trees (COMMIT)
08_identification/              # Results (COMMIT)
analysis_report.md              # Write-up (COMMIT)
```

**What's Gitignored**:
```
Large reference data (>100MB)
Build artifacts
Temporary files
Personal credentials
```

**Rationale**:
- Students commit their work (analysis outputs)
- Don't commit data that's already in template
- Reasonable repo size for GitHub Classroom
- Instructor sees their analysis, not template

---

## Documentation Strategy

### Markdown for Everything

**Decision**: All documentation in Markdown (.md files)

**Rationale**:
- Readable as plain text
- Renders nicely on GitHub
- Can convert to PDF/HTML if needed
- Students familiar with Markdown
- Version control friendly

### README Hierarchy

**Decision**:
```
README.md                  # Repository overview
XX_module/README.md        # Module documentation
XX_module/lesson.md        # Specific lessons
XX_module/exercises.md     # Problem sets
```

**Rationale**:
- Top-down navigation
- Each level self-contained
- Students can find info quickly
- Instructors can link directly to modules

### Code Examples Style

**Decision**: Every example is runnable, not pseudocode

**Rationale**:
```bash
# BAD (not actually runnable)
command --option file.txt

# GOOD (copy-paste ready)
mafft --auto sequences.fasta > aligned.fasta
```

- Students can copy-paste
- Reduces errors (they see exactly what to type)
- Teaches correct syntax
- Can test examples before publishing

---

## VS Code Integration

### Dev Container Configuration

**Decision**: Provide .devcontainer/ with full setup

**Rationale**:
- Students click "Reopen in Container"
- Everything just works
- No installation headaches
- Consistent environment (no "works on my machine")

### Extensions Included

**Decision**: Pre-install Python, Jupyter, R extensions

**Rationale**:
- Students don't know which extensions to install
- Reduces setup friction
- Quality of life improvements
- Still works without VS Code (not required)

---

## Tracking System Integration

### Minimal Tracking for This Repo

**Decision**: Only `.tracking/system.json`, no proliferation

**Rationale**:
- This is NEW repo (clean start)
- Follow tracking_system_design.md rules
- NO new tracking files
- NO ALL CAPS filenames
- Agents update system.json only

**Structure**:
```json
{
  "tracking_version": "3.0",
  "project": "dna-barcoding-analysis",
  "milestones": {
    "repo_created": "completed",
    "architecture_designed": "completed",
    ...
  }
}
```

---

## Design Trade-offs

### Completeness vs Simplicity

**Trade-off**: Make tutorials comprehensive OR keep them minimal?

**Decision**: Comprehensive (erring on side of too much detail)

**Rationale**:
- Beginners need hand-holding
- Better to skip sections than be stuck
- Can skim if experienced
- Reference value (come back later)

**Cost**: More maintenance, longer files

### Container Size vs Convenience

**Trade-off**: Minimal container (small, fast) OR batteries-included (large, slow)?

**Decision**: Batteries-included (~2GB container)

**Rationale**:
- One-time download cost
- Students don't want to install tools individually
- Instructor doesn't want support requests
- Disk space is cheap, time is expensive

**Cost**: Longer initial download

### GitHub Classroom vs Standalone

**Trade-off**: Optimize for classroom OR individual learners?

**Decision**: Both (default to standalone, Classroom-compatible)

**Rationale**:
- Most users will be self-learners
- But ENTM201L needs Classroom
- Design for standalone, add Classroom features
- Documented in main README

**Implementation**: Instructions for both use cases

---

## Future-Proofing

### Modular Structure

**Design**: Each module can be used independently

**Benefit**: If new modules needed (e.g., "10_population_genetics"), easy to add

### Container Updates

**Design**: Dockerfile uses version pinning with comments

```dockerfile
# IQ-TREE 2.2+ (update when newer version available)
iqtree=2.2.0
```

**Benefit**: Easy to update, know what versions are running

### Tool Alternatives

**Design**: Show multiple tools for same task

**Benefit**: When tool X is deprecated, students know tool Y

---

## Accessibility Considerations

### Beginner Friendly

- Zero assumptions about prior knowledge
- Explain every term (e.g., "Phred score" defined)
- Screenshots where helpful
- Error messages explained

### Multiple Learning Styles

- Text lessons (for readers)
- Code examples (for hands-on learners)
- Exercises (for practice)
- Visualizations (for visual learners)

### Self-Paced

- No hard deadlines mentioned in materials
- Optional advanced sections
- "Skip if experienced" notes
- Multiple difficulty paths

---

## Success Metrics

**How to know if architecture is good?**

1. **Student Feedback**:
   - Can complete analyses without getting stuck?
   - Find documentation helpful?
   - Understand what commands do?

2. **Adoption**:
   - Other instructors use it?
   - Lab members reference it?
   - GitHub stars/forks?

3. **Reproducibility**:
   - Same results across machines?
   - Container works on all platforms?
   - Scripts don't break?

4. **Maintenance**:
   - Easy to update tools?
   - Can add modules without breaking existing?
   - Issues get resolved quickly?

---

## Architectural Principles Summary

1. **Modularity**: Each module stands alone
2. **Clarity**: Assume zero prior knowledge
3. **Reproducibility**: Container ensures consistency
4. **Practicality**: Real data, real workflows
5. **Completeness**: Cover full pipeline (QC → ID)
6. **Flexibility**: Multiple tools, multiple paths
7. **Documentation**: Everything explained
8. **Scalability**: Works for 1 student or 100
9. **Maintainability**: Easy to update and extend
10. **Accessibility**: Designed for beginners

---

## Files Created by Architecture Design

### Documentation (17 files)
- README.md (root)
- ARCHITECTURE.md (this file)
- 00_introduction/README.md
- 01_linux_basics/README.md
- 01_linux_basics/01_navigation.md
- 01_linux_basics/06_streams_redirection.md
- 01_linux_basics/cheatsheet.md
- 02_python_basics/README.md
- 03_r_basics/README.md
- 04_data/README.md
- 05_quality_control/README.md
- 06_alignment/README.md
- 07_phylogeny/README.md
- 08_identification/README.md
- container/README.md
- .gitignore

### Directory Structure (41 directories)
All numbered module directories with subdirectories (data, exercises, solutions, etc.)

### Configuration
- .tracking/system.json (pre-existing, to be updated)
- .devcontainer/ (to be created by container setup)

---

## Implementation Notes

**What Remains To Be Done**:

1. **Detailed lessons** for 01_linux_basics (lessons 02-05, 07-08)
2. **Jupyter notebooks** for 02_python_basics
3. **R scripts** for 03_r_basics
4. **Example datasets** for all modules
5. **Exercise files** with solutions
6. **Dockerfile** (can copy from entm201l-fall2025)
7. **VS Code devcontainer.json**
8. **Helper scripts** in scripts/ directory
9. **Sample data** (reference sequences, student templates)

**Priority for Next Phase**:
1. Copy Dockerfile from entm201l-fall2025 ✓
2. Create example datasets (small, for practice)
3. Complete Linux lessons (most critical)
4. Create first Jupyter notebook (Python basics)
5. Test container build
6. Create first complete exercise with solution

---

## Conclusion

This architecture provides a solid foundation for teaching DNA barcoding analysis from first principles. The modular design allows flexibility while the comprehensive documentation ensures accessibility for beginners. The container infrastructure guarantees reproducibility. The GitHub Classroom integration enables scalable teaching.

**Design Philosophy**: "Make the easy things easy, and the hard things possible."

---

**Document Version History**:
- v1.0 (2025-11-07): Initial architecture design

**Designer**: architecture-designer-agent
**Approved by**: (pending user review)
**Status**: Design Complete, Implementation In Progress
