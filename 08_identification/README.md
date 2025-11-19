# Module 08: Species Identification

**Duration**: 2 hours
**Prerequisites**: Modules 06-07

---

## Overview

Identify unknown specimens using BLAST searches and phylogenetic placement.

---

## Methods

### 1. BLAST Against NCBI

```bash
# Online BLAST (web interface)
# https://blast.ncbi.nlm.nih.gov/

# Command-line BLAST
blastn -query unknown_sample.fasta \
       -db nt \
       -remote \
       -outfmt 6 \
       -out blast_results.txt

# Output format 6 columns:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

### 2. BOLD Identification

**Web Interface**: http://www.boldsystems.org/index.php/IDS_OpenIdEngine

**API**:
```bash
curl -F "sequences=@sequences.fasta" \
     http://www.boldsystems.org/index.php/API_Public/identify
```

### 3. Phylogenetic Placement

Place unknown sequence on reference tree and identify based on clade membership.

---

## Interpreting Results

### BLAST Results

| Percent Identity | Interpretation |
|-----------------|----------------|
| >99% | Same species (high confidence) |
| 97-99% | Likely same species |
| 95-97% | Same genus, possibly different species |
| <95% | Different genus or poor quality |

**Caution**: High similarity doesn't always mean same species!

### BOLD Results

- **Species match**: High confidence if >98% similarity
- **Genus match**: Moderate confidence
- **Family match**: Low confidence for species ID

---

## Example Workflow

```python
#!/usr/bin/env python3
"""
Identify species using BLAST
"""

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

def blast_identify(fasta_file):
    """BLAST sequence and report top hit."""
    record = SeqIO.read(fasta_file, "fasta")

    print(f"BLASTing {record.id}...")

    # Run BLAST
    result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)

    # Parse results
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)

    if blast_record.alignments:
        # Get top hit
        top_hit = blast_record.alignments[0]
        hsp = top_hit.hsps[0]

        percent_id = (hsp.identities / hsp.align_length) * 100

        print(f"Top hit: {top_hit.title}")
        print(f"Percent identity: {percent_id:.2f}%")
        print(f"E-value: {hsp.expect}")

        # Interpret
        if percent_id >= 99:
            print("Identification: HIGH CONFIDENCE")
        elif percent_id >= 97:
            print("Identification: MODERATE CONFIDENCE")
        else:
            print("Identification: LOW CONFIDENCE")
    else:
        print("No BLAST hits found")

if __name__ == "__main__":
    import sys
    blast_identify(sys.argv[1])
```

---

## Confidence Assessment

**High Confidence ID Requires**:
1. >98% sequence similarity
2. Multiple reference sequences from same species
3. Phylogenetic placement in species clade
4. Geographic range makes sense
5. Morphology consistent (if available)

**Report Identification As**:
```
Aedes aegypti (99.2% COI similarity, n=15 reference sequences)
```

---

## Dealing with Ambiguity

### Scenario 1: Multiple Top Hits
```
98.5% - Species A
98.3% - Species B
98.1% - Species C
```

**Action**: Report as "Species complex" or "cf. Species A"

### Scenario 2: No Close Matches
```
Best hit: 92% similarity
```

**Action**:
- Possibly novel species
- Check for sequencing errors
- Try different gene region
- Morphological examination required

### Scenario 3: Contradictory Results
```
BLAST: Species A (99%)
BOLD: Species B (97%)
Phylogeny: Places with Species C
```

**Action**:
- Check database errors
- Re-sequence sample
- Use multiple genes
- Consult taxonomic expert

---

## Best Practices

1. **Always Use Multiple Lines of Evidence**
   - BLAST + BOLD + Phylogeny
   - Morphology when possible

2. **Check Reference Quality**
   - Are reference sequences verified?
   - Multiple voucher specimens?

3. **Consider Geography**
   - Does species occur in collection location?
   - Check invasive species databases

4. **Document Everything**
   - Save BLAST results
   - Record accession numbers
   - Note percent similarities

---

## Final Report Template

```
Species Identification Report
=============================

Sample ID: Mosquito_001
Collection: UCR Campus, 2025-01-15

BLAST Results:
- Top hit: Aedes aegypti (GenBank: KX987654)
- Percent identity: 99.2%
- E-value: 0.0

BOLD Results:
- Match: Aedes aegypti
- Similarity: 99.1%

Phylogenetic Placement:
- Clades with A. aegypti reference clade
- Bootstrap support: 95%

Conclusion:
HIGH CONFIDENCE identification as Aedes aegypti
```

---

## Resources

- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/)
- [BOLD Systems](http://www.boldsystems.org/)
- [Mosquito Barcoding Guide](https://doi.org/10.1371/journal.pone.0029329)

---

**Analysis Complete! Now write up your results.** ✓
