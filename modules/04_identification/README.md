# Module 04: Species Identification (BLAST)

## What is BLAST?

BLAST stands for **Basic Local Alignment Search Tool**. Think of it like Google for DNA sequences!

Just like you type a search query into Google to find similar web pages, BLAST takes your DNA sequence and searches through a huge database (NCBI GenBank) to find similar sequences from known species. When it finds matches, it tells you which species your sequence is most similar to.

## How Does Species Identification Work?

Here's the simple process:

1. **You provide** a DNA sequence (from your quality-controlled reads)
2. **BLAST searches** the NCBI GenBank database (millions of known sequences)
3. **BLAST finds matches** and ranks them by similarity
4. **You get results** showing which species your sequence most likely came from

---

## How to Use It

### Command

```bash
python identify_species.py <input_fasta> [output_directory]
```

### Example

```bash
python identify_species.py results/passed_sequences.fasta results/
```

This command will:
- Read quality-controlled sequences from `passed_sequences.fasta`
- Send each sequence to NCBI for BLAST searching
- Create a nice HTML report and JSON file with results
- Show you the top 5 matches for each sequence

---

## Input Data

**Required format:** FASTA file containing DNA sequences

Each sequence should be:
- From your quality control step (Module 05)
- In FASTA format (header line starting with `>`, followed by sequence)
- Well-trimmed and quality-filtered

Example input:
```
>sample1
ATCGATCGATCGATCGATCG
>sample2
GCTAGCTAGCTAGCTAGCTA
```

---

## Output Data

The script creates two output files:

### 1. HTML Report (`identification_report.html`)
A beautiful, visual report showing:
- The top species match for each sequence
- A table of the top 5 matches for each sequence
- Color-coding for match quality:
  - **Green** = Good match (>97% identity)
  - **Yellow** = OK match (>90% identity)
  - **Red** = Poor match (<90% identity)

Open this file in your web browser to view it!

### 2. JSON Results (`identification_results.json`)
A machine-readable file with detailed numbers for each sequence:
- Accession numbers (GenBank ID)
- Exact identity percentages
- E-values (measure of statistical significance)
- Alignment information

---

## Understanding Your Results

### What Do These Numbers Mean?

#### **Percent Identity (%)**
This tells you how similar your sequence is to the match in the database.

- **97-100%** = Almost identical, very confident match
- **90-97%** = Very similar, good match for genus level
- **80-90%** = Reasonably similar, less certain
- **<80%** = Weak match, probably not the right species

**Rule of thumb:**
- `>97%` = You can confidently say it's this **species**
- `>90%` but `<97%` = You can say it's this **genus** but not the exact species
- `<90%` = Results are unclear, be cautious

#### **E-value (Expect Value)**
This is a statistical measure. Smaller is better!

- **Very small numbers** (like `1e-50` or `2e-100`) = Very significant match, extremely rare by chance
- **Medium numbers** (like `1e-10`) = Good match, unlikely to be by chance
- **Large numbers** (like `0.05` or `1`) = Not a good match

**Simple rule:** If your top hit has a small E-value and high identity, it's a good match!

#### **Coverage**
This tells you what percentage of your sequence matched the database sequence.

- **High coverage** (>90%) = Your whole sequence was similar
- **Low coverage** (<50%) = Only a small part matched, be suspicious

---

## What Makes a Good Match?

### Species-Level Identification (Very Confident)
✓ Identity **>97%**
✓ E-value **very small** (like `1e-50`)
✓ Coverage **>90%**
✓ Top hit is significantly better than second hit

**Example:** 98.5% identity to *Homo sapiens* - You can confidently say the sample is human!

### Genus-Level Identification (Confident)
✓ Identity **>90%** but **<97%**
✓ E-value **small** (like `1e-20`)
✓ Coverage **>80%**

**Example:** 92% identity to *Danio* species - You can say it's likely a zebrafish, but not which exact species.

### Questionable Results (Be Careful!)
⚠ Identity **<90%**
⚠ Low E-value or inconsistent results
⚠ Multiple species with similar scores
⚠ Low coverage

**What to do:** These results are unclear. Maybe your sample wasn't well-preserved, or it's a new/rare species not in the database.

---

## Where Do These Results Come From?

All BLAST searches are run against **NCBI GenBank**, which is:
- A massive public database of DNA sequences
- Run by the U.S. National Institutes of Health
- Updated continuously as scientists discover new species
- The gold standard for genetic databases worldwide

This means your results are compared against millions of sequences from known, published research.

---

## Troubleshooting

### Problem: "No BLAST hits found"

**What it means:** Your sequence didn't match anything in GenBank, or matched very poorly.

**Possible reasons:**
1. Your sequence is from a species that's not in the database yet
2. The sequence is damaged or has too many errors
3. The sequence is from something unusual (contamination, plasmid, etc.)

**What to do:**
- Check your sequence quality again (go back to Module 05)
- Manually check a small part of the sequence in NCBI's BLAST tool
- Ask your instructor - it might be a new species!

---

### Problem: Ambiguous or conflicting results

**What it means:** The top matches are very close in identity, or you get different species in the top 5.

**Possible reasons:**
1. Your species hasn't been fully sequenced/documented yet
2. The DNA region you're using is very similar across multiple species
3. The sample might be a hybrid or new species

**What to do:**
- Look at ALL top 5 matches - do they all belong to the same genus?
- If yes, you can confidently identify genus
- If no, the results are unclear - report this in your lab notes!

---

### Problem: Very high identity (99-100%) to unexpected species

**What it means:** Your sample is nearly identical to a different species than expected.

**Possible reasons:**
1. Your sample was contaminated with a different species
2. It's actually a different species than you thought
3. Lab mix-up or mislabeling

**What to do:**
- Re-check your lab notes and sample handling
- Check if the species actually make sense (e.g., are they found in the same location?)
- Ask your instructor for guidance

---

## Next Steps

After identifying your sequences:

1. **Record your results** in your lab notebook
2. **Summarize the matches** (which species, which confidence level)
3. **Move to the next step** (usually phylogenetic analysis)
4. **See the master pipeline** (`master_pipeline.py`) to learn how all modules work together

For advanced analysis, you might want to:
- Run phylogenetic trees with your identified sequences
- Compare multiple samples to see evolutionary relationships
- Search for additional barcoding regions for ambiguous results

---

## Quick Reference

| Metric | Good Match | OK Match | Poor Match |
|--------|-----------|----------|-----------|
| Identity | >97% | 90-97% | <90% |
| E-value | Very small | Small | Large or absent |
| Coverage | >90% | >80% | <80% |
| Confidence | Species | Genus | Unreliable |

---

## Questions?

- Check the code comments in `identify_species.py`
- Look at example output files in the `examples/` directory
- Ask your instructor or teaching assistant
- See the main README.md for the full course overview
