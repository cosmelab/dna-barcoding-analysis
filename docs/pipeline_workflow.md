# DNA Barcoding Pipeline - Visual Workflow

## Overview: From Chromatograms to Species Identification

```
┌─────────────────────────────────────────────────────────────────────┐
│                    DNA BARCODING WORKFLOW                            │
│                   (5 Steps to Species ID)                           │
└─────────────────────────────────────────────────────────────────────┘

YOUR SAMPLES                          WHAT HAPPENS
(from UC Genomics)                    (This pipeline does it)
─────────────────                     ──────────────────────


📊 Raw Data
─────────────
.ab1 files
(chromatograms)                       ┌──────────────────────┐
  • AT-HV1F.ab1 ────┐                │  STEP 1: QC          │
  • AT-HV1R.ab1 ────┼──────────────> │  Quality Control     │
  • AT-HV3F.ab1 ────┤                └──────────────────────┘
  • AT-HV3R.ab1 ────┘                         │
  • (30 total files)                          │ Checks quality & length
                                              ▼
                                       ✓ PASSED: 12 sequences
                                       ✗ FAILED: 18 sequences
                                              │
                                              ▼
🧬 Passed Sequences
──────────────────
Individual F & R reads                ┌──────────────────────┐
>AT-HV1F                             │  STEP 2: CONSENSUS   │
ATCGATCGATCG...       ──────────────> │  Combine F+R pairs   │
>AT-HV1R                             └──────────────────────┘
ATCGATCGATCG...                              │
>AT-HV3F                                     │ Pairs F+R reads
ATCGATCGATCG...                              │ Creates consensus
>AT-HV3R                                     ▼
ATCGATCGATCG...                  ✓ 4 consensus sequences created
                                 ✗ 4 samples missing F or R
                                              │
                                              ▼
🎯 Consensus Sequences
─────────────────────
Best sequence from F+R            ┌──────────────────────────┐
>AT-HV1 (consensus)               │  STEP 3: COMBINE         │
ATCGATCGATCG...       ────────────>│  Add reference sequences │
>AT-HV3 (consensus)               └──────────────────────────┘
ATCGATCGATCG...                            │
>AT-JM2 (consensus)                        │ Adds 52 known
ATCGATCGATCG...                            │ SoCal mosquitoes
>AT-WL2 (consensus)                        ▼
ATCGATCGATCG...               56 sequences (4 yours + 52 refs)
                                           │
                                           ▼
📐 Combined Dataset
──────────────────
Your samples + references         ┌──────────────────────┐
(56 sequences total)              │  STEP 4: ALIGN+TREE  │
                  ──────────────> │  Compare all         │
                                  └──────────────────────┘
                                           │
                                           │ MAFFT aligns
                                           │ IQ-TREE builds tree
                                           ▼
                                  Tree shows where YOUR
                                  samples cluster with
                                  known species!
                                           │
                                           ▼
🌳 Phylogenetic Tree                ┌──────────────────────┐
────────────────────               │  STEP 5: BLAST       │
Your 4 consensus     ──────────────>│  Identify species    │
sequences                          └──────────────────────┘
                                           │
                                           │ Compares to
                                           │ GenBank database
                                           ▼
🏷️  Species Identified!
────────────────────
AT-HV1  → Aedes albopictus (99.55%) - Asian tiger mosquito
AT-HV3  → Culex pipiens (98.12%)     - Northern house mosquito
AT-JM2  → Culex pipiens (99.25%)     - Northern house mosquito
AT-WL2  → Culex pipiens (98.67%)     - Northern house mosquito


📊 Final Results: HTML Reports + Tree Figures + Summary Table
```

## What You'll Learn at Each Step

### Step 1: Quality Control (5 minutes)
**QUESTION:** Are my sequences good enough?

**YOU LEARN:**
- What makes a good Sanger sequence
- Why some sequences fail (low quality, too short)
- How to read chromatograms
- Why you need BOTH forward and reverse reads

**YOU DO:**
- Run one simple command
- Look at HTML report with chromatogram visualizations
- Count how many F and R reads passed

**KEY CONCEPT:** Not all sequences from the core are usable - that's normal! We need high quality data for accurate species ID.

---

### Step 2: Consensus Sequences (3 minutes)
**QUESTION:** How do we combine forward and reverse reads?

**YOU LEARN:**
- Why forward and reverse reads are sequenced
- How consensus sequences improve accuracy
- What happens if only F or only R passes QC

**YOU DO:**
- Run one simple command with `--pairs-only` flag
- View HTML report showing which samples have complete pairs
- See how many samples make it to the next step

**KEY CONCEPT:** Consensus from F+R is more accurate than a single read! Only complete pairs are used.

---

### Step 3: Alignment & Tree (10 minutes)
**QUESTION:** How do my samples compare to known mosquito species?

**YOU LEARN:**
- What "alignment" means (lining up DNA sequences)
- How phylogenetic trees show evolutionary relationships
- How to interpret your samples clustering with references

**YOU DO:**
- Run alignment + tree commands (combined with 52 reference sequences)
- View tree showing where YOUR 4 samples cluster
- Identify which reference species are closest to yours

**KEY CONCEPT:** The tree shows you which known species your samples are most similar to!

---

### Step 4: Species Identification (5 minutes)
**QUESTION:** What species are my samples?

**YOU LEARN:**
- How BLAST works (compares to GenBank database)
- What % identity means (>98% = same species)
- How to interpret top hits and write scientific names correctly

**YOU DO:**
- Run BLAST command on your consensus sequences
- Read HTML report with top matches
- Fill in summary table with species names and % identity

**KEY CONCEPT:** BLAST confirms what the tree suggests! >98% match = probably the same species.

---

## Time Investment

**First time (with tutorial):** ~30 minutes
- 10 min learning workflow
- 20 min running test data

**Your real data:** ~10 minutes
- Just run the commands
- Focus on interpreting results

**Writing up results:** ~30 minutes
- Answer assignment questions
- Include figures in report
