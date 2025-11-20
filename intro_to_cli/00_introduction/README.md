# Module 00: Introduction to DNA Barcoding

**Duration**: 30 minutes
**Prerequisites**: None

---

## What is DNA Barcoding?

DNA barcoding is a method of species identification using a short, standardized genetic marker. Think of it as a molecular barcode scanner for biodiversity.

**The Concept**:
- Each species has unique DNA sequences
- A standard gene region acts as the "barcode"
- Sequencing this region identifies unknown specimens
- Like scanning a product barcode at a store

---

## The COI Gene: Nature's Barcode

### Why Cytochrome c Oxidase Subunit I (COI)?

**For Animals (especially insects)**:
- Universal across animal kingdom
- Highly variable between species
- Conserved within species
- Easy to amplify with PCR
- ~650 base pairs (perfect length)

**The "Folmer Region"**:
- Named after Folmer et al. (1994)
- Standard primers: LCO1490 and HCO2198
- Works across diverse animal groups
- Over 30 years of sequence data available

---

## Applications in Entomology

### 1. Species Identification
- Identify mosquito species
- Distinguish cryptic species (look identical)
- Identify larvae (can't use morphology)

### 2. Biodiversity Surveys
- Rapid assessment of insect communities
- Monitor invasive species
- Track seasonal changes

### 3. Food Safety
- Detect insect contamination
- Verify species in products

### 4. Forensics
- Link suspects to locations (insects as evidence)
- Determine time of death (carrion insects)

---

## The DNA Barcoding Workflow

```
Specimen Collection
        ↓
DNA Extraction ---------------→ Module 01-02 (ENTM201L)
        ↓
PCR Amplification ------------→ Module 07 (ENTM201L)
        ↓
Sanger Sequencing ------------→ Module 09 (ENTM201L)
        ↓
Quality Control --------------→ Module 05 (This repository)
        ↓
Sequence Alignment -----------→ Module 06 (This repository)
        ↓
Phylogenetic Analysis --------→ Module 07 (This repository)
        ↓
Species Identification -------→ Module 08 (This repository)
```

---

## Real-World Example

**Scenario**: You caught a mosquito in your backyard. Is it *Aedes aegypti* (disease vector)?

1. **Extract DNA** from one leg
2. **PCR amplify** COI gene
3. **Sequence** the amplicon
4. **Compare** to database (BOLD, GenBank)
5. **Result**: 99.8% match to *Aedes aegypti*
6. **Conclusion**: Yes, it's a disease vector!

---

## Key Databases

### BOLD Systems (Barcode of Life Data Systems)
- http://www.boldsystems.org/
- Dedicated to DNA barcoding
- Curated sequences with specimen data
- Identification engine
- Over 10 million barcode sequences

### GenBank (NCBI)
- https://www.ncbi.nlm.nih.gov/genbank/
- General sequence database
- BLAST search tool
- Broader coverage but less curated

---

## COI Gene Structure

```
Mitochondrial DNA (circular)
     |
     ├── COI gene (~1500 bp)
     |       |
     |       └── Folmer region (~650 bp) ← THIS IS THE BARCODE
     ├── COII gene
     ├── Cytochrome b
     └── Other genes
```

**Properties**:
- Mitochondrial (high copy number = easy to amplify)
- Protein-coding (less sequencing errors than non-coding)
- Maternal inheritance (no recombination)
- Fast evolution (species-specific variation)

---

## What Makes a Good Barcode?

1. **Universal Primers** - Work across many species
2. **Variability** - Different between species
3. **Conserved** - Similar within species
4. **Short** - Can sequence with Sanger
5. **Easy to Amplify** - Reliable PCR

COI checks all boxes for animals!

---

## Limitations

### When DNA Barcoding Doesn't Work Well

1. **Very Recent Speciation**
   - Species too similar genetically

2. **Hybridization**
   - Hybrid specimens have mixed signals

3. **Incomplete Databases**
   - Species not yet sequenced and deposited

4. **Poor DNA Quality**
   - Degraded specimens (old museum samples)

5. **Symbionts/Parasites**
   - Might amplify wrong organism

**Solution**: Use multiple approaches (morphology + barcoding)

---

## This Repository's Focus

We'll take you from **raw Sanger sequences** to **species identification**:

1. **Quality Control** - Clean up sequencing artifacts
2. **Alignment** - Compare multiple sequences
3. **Phylogeny** - Build evolutionary tree
4. **Identification** - Match to known species

All using open-source tools and reproducible workflows!

---

## Key References

**Original DNA Barcoding Paper**:
Hebert, P.D.N., Cywinska, A., Ball, S.L., & deWaard, J.R. (2003). Biological identifications through DNA barcodes. *Proceedings of the Royal Society B*, 270(1512), 313-321. [DOI: 10.1098/rspb.2002.2218](https://doi.org/10.1098/rspb.2002.2218)

**COI Primers**:
Folmer, O., Black, M., Hoeh, W., Lutz, R., & Vrijenhoek, R. (1994). DNA primers for amplification of mitochondrial cytochrome c oxidase subunit I from diverse metazoan invertebrates. *Molecular Marine Biology and Biotechnology*, 3(5), 294-299.

**Insect Barcoding Review**:
Meier, R., Shiyang, K., Vaidya, G., & Ng, P.K.L. (2006). DNA barcoding and taxonomy in Diptera: A tale of high intraspecific variability and low identification success. *Systematic Biology*, 55(5), 715-728.

---

## Let's Begin!

Now that you understand DNA barcoding, proceed to:

**Module 01: Linux Basics** - Learn command-line skills needed for bioinformatics

---

**Questions? Proceed to Module 01!** →
