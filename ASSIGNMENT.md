# Week 8 Lab: DNA Barcoding Analysis

**Due**: End of Week 8
**Points**: 100

---

## 🎯 Objectives

By the end of this lab, you will:
- Analyze Sanger sequencing chromatograms
- Perform quality control on DNA sequences
- Build phylogenetic trees
- Identify mosquito species using COI barcodes

---

## 📋 Instructions

### 1. Setup (10 points)

1. Accept this GitHub Classroom assignment
2. Clone your repository to your computer
3. Verify Docker is installed:
   ```bash
   docker --version
   docker-compose --version
   ```

**Deliverable**: Repository cloned with `.ab1` files in `data/my_sequences/`

### 2. Upload Your Chromatograms (10 points)

Place your `.ab1` chromatogram files from the UC Genomics Core in the `data/my_sequences/` directory.

**Required**:
- Minimum 2 samples (Forward + Reverse for at least 1 mosquito specimen)
- Files must be `.ab1` format
- Follow naming convention: `SampleID_F.ab1` and `SampleID_R.ab1`

### 3. Run the Analysis (60 points)

Start the container and run the complete pipeline:

```bash
# Start container
docker-compose up -d
docker-compose exec dna-barcoding zsh

# Inside container
analyze-sequences

# Exit
exit
docker-compose down
```

**Deliverables** (auto-graded):
- QC Report (20 points): `results/run_*/01_qc/qc_report.html`
- Passed Sequences (20 points): `results/run_*/01_qc/passed_sequences.fasta`
- Phylogenetic Tree (25 points): `results/run_*/03_phylogeny/tree.png`
- Species ID (25 points): `results/run_*/04_identification/identification_report.html`

### 4. Interpret Results (20 points)

Create a file `RESULTS.md` in your repository answering these questions:

1. **Quality Control** (5 points)
   - How many sequences passed QC?
   - What was the average quality score?
   - Did any sequences fail? Why?

2. **Phylogenetic Analysis** (5 points)
   - Which species are most closely related to your samples?
   - What is the bootstrap support for the placement of your samples?

3. **Species Identification** (10 points)
   - What species did BLAST identify your samples as?
   - What was the % identity?
   - Do the BLAST results agree with the phylogenetic tree?
   - Are you confident in the species ID? Why or why not?

---

## 📤 Submission

1. Commit all results to your repository:
   ```bash
   git add data/ results/ RESULTS.md
   git commit -m "Complete Week 8 DNA barcoding analysis"
   git push
   ```

2. Verify your submission on GitHub - you should see:
   - ✅ Your `.ab1` files in `data/my_sequences/`
   - ✅ Results in `results/run_*/`
   - ✅ `RESULTS.md` with your interpretations

3. Autograding will run automatically and assign points

---

## 🆘 Troubleshooting

**No sequences passed QC?**
- Check your chromatogram quality in the QC report
- Low quality is common at sequence ends - this is normal
- Contact the instructor if all sequences fail

**BLAST returns no results?**
- Your sequence may be from a species not in the database
- Check for contamination
- Try aligning to the reference sequences manually

**Docker errors?**
- Make sure Docker Desktop is running
- Try: `docker-compose down` then `docker-compose up -d`
- Restart Docker Desktop if needed

**Need help?**
- Office hours: [Schedule]
- Email: [Instructor email]
- GitHub Issues: Use for technical problems

---

## 📚 Grading Rubric

| Component | Points | Criteria |
|-----------|--------|----------|
| Setup | 10 | Repository cloned, chromatograms uploaded |
| QC Report | 20 | HTML report generated with quality metrics |
| Passed Sequences | 20 | FASTA file with QC-passed sequences |
| Phylogeny | 25 | Tree built and visualized |
| Species ID | 25 | BLAST identification complete |
| Interpretation | 20 | RESULTS.md with thoughtful answers |
| **Total** | **100** | |

---

## 🎓 Learning Goals

This lab teaches you:
- Real-world bioinformatics workflows
- Interpreting sequence quality
- Phylogenetic relationships
- Species identification methods
- Reproducible research practices (Docker, Git)

Good luck! 🦟🧬
