# COI Reference Sequences for Southern California Mosquitoes

## Overview

This directory contains COI (Cytochrome c Oxidase subunit I) reference sequences for mosquito species, with a focus on Southern California genera (Aedes, Culex, Anopheles).

**Dataset:** 52 COI sequences from 35 mosquito species
**Source:** GenBank accessions from Hoque et al. 2022
**Gene region:** Mitochondrial COI barcode region (~600-1500 bp)
**Primary file:** `socal_mosquitoes.fasta`

## Citation

All GenBank accessions were extracted from:

**Hoque MM, Valentine MJ, Kelly PJ, Barua S, Barrantes Murillo DF, Wang C.** (2022)
*Modification of the Folmer primers for the cytochrome c oxidase gene facilitates identification of mosquitoes.*
Parasites & Vectors 15:437.
DOI: [10.1186/s13071-022-05494-2](https://doi.org/10.1186/s13071-022-05494-2)

### Paper Summary

The Hoque et al. 2022 study developed improved COI primers (AUCOS - Auburn University COI System) for mosquito identification, which showed superior performance compared to the traditional Folmer primers (FCOS). The AUCOS system:
- Successfully amplified 67.5% of mosquito samples vs 16.7% with FCOS
- Provided higher quality sequence data
- Correctly identified Aedes, Anopheles, Deinocerites, and Uranotaenia species
- Had challenges with Culex species identification (common issue in mosquito barcoding)

## Southern California Priority Species

The following species are of particular relevance to Southern California vector surveillance:

### Aedes (5 species, 11 sequences)
- **Aedes aegypti** (Yellow fever mosquito) - 5 sequences
  - Vector: Dengue, Zika, Chikungunya, Yellow fever
  - GenBank: MN299002.1, MK300221.1, MN298992.1, MN298993.1, MK300224.1
- **Aedes albopictus** (Asian tiger mosquito) - 2 sequences
  - Vector: Dengue, Zika, Chikungunya
  - GenBank: KF211505, MN513368
- **Aedes taeniorhynchus** (Black salt marsh mosquito) - 1 sequence
  - Coastal pest species
  - GenBank: MN626442.1

### Culex (7 species, 16 sequences)
- **Culex pipiens** (Northern house mosquito) - 4 sequences
  - Vector: West Nile virus, Western equine encephalitis
  - GenBank: KP293422.1, MK714012.1, KP293425.1, MK714001.1
- **Culex quinquefasciatus** (Southern house mosquito) - 4 sequences
  - Vector: West Nile virus, St. Louis encephalitis
  - GenBank: MW509603, MH463059.1, MN389462.1, MN005046.1
- **Culex tarsalis** (Western encephalitis mosquito) - 2 sequences
  - Primary vector: West Nile virus in California
  - GenBank: NC_036006.1, AF425847
- **Culex erraticus** - 3 sequences
  - GenBank: MH129001, MH128999.1, MN389459
- **Culex nigripalpus**, **C. usquatissimus**, **C. sitiens**, **C. tritaeniorhynchus**, **C. coronator** - 1 sequence each

### Anopheles (2 species, 3 sequences)
- **Anopheles quadrimaculatus** (Common malaria mosquito) - 2 sequences
  - Historical malaria vector in USA
  - GenBank: L04272, L04272.1
- **Anopheles punctipennis** - 1 sequence
  - Former malaria vector
  - GenBank: KR666470.1

## Additional Reference Species

The dataset includes comparative sequences from other genera:

### Other Aedes species (5 species, 7 sequences)
- Aedes busckii, A. japonicus, A. tortilis, A. vexans, A. triseriatus

### Other Anopheles species (5 species, 5 sequences)
- A. crucians, A. albimanus, A. funestus, A. gambiae, A. pseudopunctipennis, A. stephensi

### Other genera (8 species, 10 sequences)
- **Deinocerites magnus** (1 seq)
- **Psorophora** (5 species: howardii, ferox, pygmaea, cingulata, confinis) - 6 sequences
- **Uranotaenia sapphirina** (1 seq)
- **Culiseta annulata** (1 seq)
- **Mansonia annulata** (1 seq)
- **Ochlerotatus taeniorhynchus** (1 seq)

## File Format

**File:** `socal_mosquitoes.fasta`
**Format:** FASTA nucleotide sequences
**Headers:** `>Genus_species_AccessionNumber [Organism] [Location]`

Example:
```
>Aedes_aegypti_MN299002.1 Aedes aegypti
CATTGGTCAACAAATCATAAAGATATTGGAACTTTATATTTCATTTTTGGAGTATGATCTGGAATAGTCGGAACTTCTCT...
```

## Sequence Statistics

| Statistic | Value |
|-----------|-------|
| Total sequences | 52 |
| Total species | 35 |
| Southern California priority species | 9 |
| Sequence length range | 639-2306 bp |
| Average sequence length | ~1100 bp |
| Failed retrievals | 1 (KR653634.100) |
| Success rate | 98.1% |

### Sequences by Genus

| Genus | Species | Sequences | SoCal Priority |
|-------|---------|-----------|----------------|
| Aedes | 8 | 13 | 3 species |
| Anopheles | 7 | 8 | 2 species |
| Culex | 7 | 16 | 4 species |
| Psorophora | 5 | 6 | 0 |
| Deinocerites | 1 | 1 | 0 |
| Uranotaenia | 1 | 1 | 0 |
| Culiseta | 1 | 1 | 0 |
| Mansonia | 1 | 1 | 0 |
| Ochlerotatus | 1 | 1 | 0 |

## Usage

This reference dataset can be used for:

1. **BLAST searches** - Identify unknown mosquito COI sequences
2. **Phylogenetic analysis** - Build trees showing evolutionary relationships
3. **Primer design** - Design species-specific or genus-specific primers
4. **Barcode gap analysis** - Assess intra- vs inter-species variation
5. **Vector surveillance** - Molecular identification of field-collected specimens

### Example BLAST Usage

```bash
# Create BLAST database
makeblastdb -in socal_mosquitoes.fasta -dbtype nucl -parse_seqids

# Search unknown sequence
blastn -query unknown_mosquito.fasta -db socal_mosquitoes.fasta -outfmt 6
```

### Example Multiple Sequence Alignment

```bash
# Using MAFFT
mafft --auto socal_mosquitoes.fasta > socal_mosquitoes_aligned.fasta

# Using Clustal Omega
clustalo -i socal_mosquitoes.fasta -o socal_mosquitoes_aligned.fasta --threads=4
```

## Data Collection

**Collection date:** November 2025
**Retrieval method:** BioPython Entrez API
**Script:** `scripts/fetch_reference_sequences.py`

All sequences were downloaded directly from GenBank using accession numbers documented in Tables 2 and 3 of Hoque et al. 2022.

## Notes and Limitations

1. **Culex species identification:** As noted in the original paper, COI barcoding has challenges distinguishing between closely related Culex species, particularly within the Cx. pipiens complex (pipiens, quinquefasciatus, and their hybrids).

2. **Geographic variation:** These sequences represent global samples. Local California populations may show some sequence variation.

3. **Sequence lengths:** Sequences vary in length (639-2306 bp) as some represent:
   - Partial COI barcode region (~650 bp)
   - Standard barcode region (~700 bp)
   - Complete COI gene (~1500 bp)
   - Extended COI sequences (>2000 bp)

4. **Missing accession:** One accession (KR653634.100) failed to download due to a malformed accession number in the original paper.

5. **Taxonomic updates:** Mosquito taxonomy is actively being revised. Some species names or classifications may have been updated since original sequence submission.

## Quality Control

- All sequences were retrieved directly from GenBank
- Sequences verified to match accessions in Hoque et al. 2022
- COI gene regions extracted when annotated in GenBank records
- FASTA headers include organism and location metadata when available

## Future Updates

This dataset may be expanded to include:
- Additional California-specific specimens
- More recent sequence submissions
- Sequences from other relevant COI studies
- Population-level variants from Southern California

## Contact

For questions about this dataset or the DNA barcoding pipeline, please refer to the main project documentation.

## References

1. Hoque MM, Valentine MJ, Kelly PJ, Barua S, Barrantes Murillo DF, Wang C. (2022) Modification of the Folmer primers for the cytochrome c oxidase gene facilitates identification of mosquitoes. Parasites & Vectors 15:437. https://doi.org/10.1186/s13071-022-05494-2

2. Folmer O, Black M, Hoeh W, Lutz R, Vrijenhoek R. (1994) DNA primers for amplification of mitochondrial cytochrome c oxidase subunit I from diverse metazoan invertebrates. Molecular Marine Biology and Biotechnology 3:294-299.

3. Hebert PDN, Cywinska A, Ball SL, DeWaard JR. (2003) Biological identifications through DNA barcodes. Proceedings of the Royal Society B 270:313-321.

---

**Last updated:** November 2025
**Dataset version:** 1.0
**License:** Data from public GenBank records
