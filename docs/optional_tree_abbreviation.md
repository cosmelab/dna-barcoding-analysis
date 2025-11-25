# Optional: Tree Name Abbreviation Tool

## What Is This?

An optional tool to shorten long species names in phylogenetic trees for cleaner visualization.

**Example:**
- Original: `Culex_pipiens_KP293422.1` (25 characters)
- Abbreviated: `cu_pip_KP29` (11 characters)
- **Result:** 57% shorter names, cleaner tree visualizations

## When To Use

**Use abbreviations if:**
- Tree labels are crowded or overlap
- Creating publication figures
- Printing trees where space is limited
- You want a cleaner appearance

**Skip abbreviations if:**
- Current tree is readable (default is usually fine!)
- This is your first analysis (keep it simple)
- You need full species names for reference

## How To Use

**Step 1: Run normal analysis first**
```bash
./run-analysis.sh
# Creates: results/your_analysis/04_phylogeny/tree.treefile
```

**Step 2: (OPTIONAL) Create abbreviated version**
```bash
python scripts/abbreviate_tree_names.py results/your_analysis/04_phylogeny/tree.treefile
```

**Output:**
- `tree_abbreviated.treefile` - Shortened names (open in FigTree)
- `abbreviation_mapping.csv` - Reference table (full name → abbreviation)
- `abbreviation_mapping.md` - Human-readable reference

## Abbreviation Format

**Pattern:** `GG_sss_AAAA`
- `GG` = First 2 letters of genus (lowercase)
- `sss` = First 3 letters of species (lowercase)
- `AAAA` = First 4 characters of accession number

**Examples:**
| Original | Abbreviated | Reduction |
|----------|-------------|-----------|
| Culex_pipiens_KP293422.1 | cu_pip_KP29 | 56% shorter |
| Aedes_aegypti_MN299002.1 | ae_aeg_MN299 | 56% shorter |
| Anopheles_gambiae_MG753769 | an_gam_MG753 | 62% shorter |

**Your samples preserved:**
| Original | Abbreviated | Notes |
|----------|-------------|-------|
| AT-HV1 | AT-HV1 | Unchanged (student sample) |
| AT99 | AT99 | Unchanged (student sample) |
| AT_ROCK_ | AT_ROCK_ | Unchanged (student sample) |

## Viewing Abbreviated Trees

**In FigTree (recommended):**
1. Download FigTree: http://tree.bio.ed.ac.uk/software/figtree/
2. Open `tree_abbreviated.treefile`
3. Refer to `abbreviation_mapping.csv` for full names

**Reference table:**
```csv
Original,Abbreviated,Genus,Species,Accession
Culex_pipiens_KP293422.1,cu_pip_KP29,Culex,pipiens,KP293422.1
...
```

## Tips

- **First time?** Skip abbreviations - default tree is fine!
- **Publication?** Use abbreviated version for cleaner figures
- **Confused?** Stick with original full names - they're self-explanatory

## Technical Details

For developers and instructors, see:
- `scripts/abbreviate_tree_names.py` - Source code
- Automatic generation ensures uniqueness
- Maintains all student sample names unchanged
- 100% reversible (mapping table provided)

---

**Bottom line:** This is an optional enhancement. The default tree with full names works great for most students!
