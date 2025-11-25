# Reference Sequence Trimming

## Why Do We Trim Reference Sequences?

### The COI Barcode Region

When we use DNA barcoding to identify mosquito species, we sequence a specific ~700 base pair (bp) region of the **cytochrome c oxidase I (COI)** gene. This region is called the "barcode region" because it's:

1. **Variable enough** to distinguish between species
2. **Conserved enough** to amplify with universal primers
3. **Short enough** to sequence completely in one reaction

### The Problem with Untrimmed References

GenBank (NCBI's DNA database) contains COI sequences of varying lengths:

- **Barcode sequences:** ~640-750 bp (just the barcode region)
- **Partial COI genes:** ~1,500 bp (larger portion of the gene)
- **Complete COI genes:** ~2,300 bp (entire gene)

When we downloaded reference sequences from GenBank, we got a **mix of all three types**.

### What Happens When References Are Too Long?

If you try to align a 700bp student sequence with 2,300bp reference sequences:

```
Student sequence:    ATCG... (700bp)
Long reference:      ATCG..........................................ATCG (2,300bp)
```

MAFFT (the alignment program) adds **gaps** to make everything line up:

```
Aligned student:     ATCG-------------------------- (lots of gaps!)
Aligned reference:   ATCG...full sequence...ATCG
```

**Result:** Our alignment was **2,334 bp** with ~1,000 gaps! This is bad because:

- Phylogenetic trees become unreliable with too many gaps
- BLAST works better with sequences of similar length
- Processing takes longer
- Results are harder to interpret

### The Solution: Trim to the Barcode Region

We created a Python script (`data/reference_sequences/trim_references_to_barcode.py`) that:

1. **Finds the barcode region** in each reference using the AUCOS primer sequences from [Hoque et al. 2022](https://doi.org/10.1186/s13071-022-05494-2)
2. **Trims** references longer than 800bp down to ~712bp (the AUCOS amplicon size)
3. **Keeps** references that are already the right size (640-800bp)
4. **Removes** references that are too short (<640bp)

### Results After Trimming

| Metric | Before Trimming | After Trimming |
|--------|----------------|----------------|
| Reference sequences | 52 | 50 |
| Average length | 1,123 bp | 700 bp |
| Length range | 639 - 2,306 bp | 657 - 750 bp |
| **Alignment length** | **2,334 bp** | **784 bp** |
| Gaps introduced | ~1,000 | ~80 |

### Scientific Basis

According to **Hoque et al. 2022** (Parasites & Vectors):
- The AUCOS primers amplify a **712 bp** COI barcode region
- This is the same variable region used by the original Folmer primers
- GenBank sequences can be used for identification if trimmed to this region

### For Students

**You don't need to do any trimming!** We've already:
- ✅ Trimmed all reference sequences to the barcode region
- ✅ Verified they're the correct length (~700bp)
- ✅ Tested the alignment (784bp, perfect!)

Your student sequences from UC Genomics are already the right length (~700bp) because the PCR primers only amplify the barcode region.

### Technical Details

The trimming script:
- Uses the AUCOS primers: `AU-COI-F` and `AU-COI-R` (Table 1, Hoque et al. 2022)
- Allows up to 3 mismatches for degenerate base matching
- Handles IUPAC degenerate codes (W, R, Y, etc.)
- Preserves sequence metadata (species names, GenBank IDs)

### Files

- **Original references:** `data/reference_sequences/socal_mosquitoes_original.fasta` (backup)
- **Trimmed references:** `data/reference_sequences/socal_mosquitoes.fasta` (active)
- **Trimming script:** `data/reference_sequences/trim_references_to_barcode.py`

## Reference

Hoque MM, Valentine MJ, Kelly PJ, Barua S, Murillo DFB, Wang C. (2022). Modification of the Folmer primers for the cytochrome c oxidase gene facilitates identification of mosquitoes. *Parasites & Vectors* 15:437. https://doi.org/10.1186/s13071-022-05494-2
