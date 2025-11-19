# Module 04: Reference Data and Datasets

**Duration**: 1 hour
**Prerequisites**: Module 01

---

## Overview

Learn how to download, organize, and manage reference sequence data for DNA barcoding analysis.

---

## Directory Structure

```
04_data/
├── README.md                   # This file
├── reference_sequences/        # Published data
│   ├── socal_insects/         # Southern California dataset
│   ├── ncbi_coi/              # Downloaded from GenBank
│   └── bold_public/           # Downloaded from BOLD
├── student_sequences/          # Your sequences go here
│   ├── .gitkeep
│   └── metadata_template.csv
└── metadata/                   # Specimen information
    └── data_dictionary.md
```

---

## Data Sources

### 1. Southern California Insect Dataset

**Citation**: [To be added - manuscript in preparation]

**Contents**:
- COI barcode sequences from multiple insect orders
- Southern California biodiversity focus
- Voucher specimen data
- GPS coordinates
- Collection dates

**Location**: `reference_sequences/socal_insects/`

### 2. NCBI GenBank

**Download COI Sequences**:
```bash
# Using Entrez Direct tools
esearch -db nucleotide -query "COI[Gene] AND Aedes[Organism]" |
efetch -format fasta > aedes_coi.fasta

# Download specific accessions
efetch -db nucleotide -id "accession1,accession2" -format fasta > sequences.fasta
```

### 3. BOLD Systems

**Web Interface**: http://www.boldsystems.org/

**API Download**:
```bash
# Download sequences for specific taxon
wget "http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Aedes" \
  -O aedes_bold.fasta
```

---

## Student Data Integration

### File Naming Convention

```
LastName_SpecimenID_COI.fasta
LastName_SpecimenID_COI.ab1
```

**Examples**:
```
Cosme_Mosquito01_COI.fasta
Cosme_Mosquito01_COI.ab1
Smith_Beetle05_COI.fasta
```

### Metadata Template

**File**: `student_sequences/metadata_template.csv`

```csv
student_name,specimen_id,collection_date,location,latitude,longitude,habitat,notes
Cosme,Mosquito01,2025-01-15,UCR Campus,33.9737,-117.3281,Urban,Found near fountain
```

---

## Data Organization Best Practices

1. **Never Edit Raw Data** - Always work on copies
2. **Use Descriptive Filenames** - Include date, sample ID
3. **Document Everything** - Metadata is crucial
4. **Version Control** - Use Git for tracking
5. **Backup Regularly** - Data loss is catastrophic

---

## Example: Download Reference Data

```bash
#!/bin/bash
# Download COI sequences for common mosquito species

# Create output directory
mkdir -p reference_sequences/ncbi_coi

# Download Aedes aegypti
esearch -db nucleotide -query "Aedes aegypti[Organism] AND COI[Gene] AND 600:700[Sequence Length]" |
efetch -format fasta > reference_sequences/ncbi_coi/aedes_aegypti_coi.fasta

# Download Culex pipiens
esearch -db nucleotide -query "Culex pipiens[Organism] AND COI[Gene] AND 600:700[Sequence Length]" |
efetch -format fasta > reference_sequences/ncbi_coi/culex_pipiens_coi.fasta

# Combine all
cat reference_sequences/ncbi_coi/*.fasta > reference_sequences/ncbi_coi/all_mosquito_coi.fasta

echo "Download complete!"
```

---

## Next Steps

Proceed to **Module 05: Quality Control** to learn how to process your sequences.

---
