# Data Files for R Basics Module

This directory contains example data files for practicing R phylogenetic analysis.

## Files

### example_tree.newick
- **Format**: Newick/New Hampshire tree format
- **Content**: Phylogenetic tree of 15 mosquito COI sequences
- **Taxa**: 15 samples from 3 genera (Aedes, Culex, Anopheles)
- **Branch lengths**: Yes (substitutions per site)
- **Bootstrap values**: Yes (stored as node labels)
- **Use**: Main tree for all exercises

### bootstrap_tree.treefile
- **Format**: Newick with bootstrap support values
- **Content**: Same tree as example_tree.newick
- **Bootstrap**: 1000 replicates
- **Use**: Practice working with support values

### tree_metadata.csv
- **Format**: CSV (comma-separated values)
- **Columns**:
  - `label`: Tip label matching tree (e.g., Aedes_aegypti_CA_2023)
  - `genus`: Genus name (Aedes, Culex, Anopheles)
  - `species`: Species epithet (aegypti, pipiens, gambiae, etc.)
  - `location`: Collection location (state/country)
  - `country`: Country code
  - `year`: Collection year
  - `collector`: Collector name
  - `disease_vector`: Associated disease
  - `habitat`: Habitat type (Urban, Rural, etc.)
  - `sample_quality`: Quality assessment (High, Medium, Low)
- **Rows**: 15 (one per sample)
- **Use**: Practice joining metadata with trees for colored/annotated plots

### trait_data.csv
- **Format**: CSV
- **Columns**:
  - `label`: Sample ID matching tree tips
  - `wing_length_mm`: Wing length in millimeters
  - `body_size_mm`: Body size in millimeters
  - `biting_rate`: Biting frequency (arbitrary units)
  - `flight_distance_km`: Maximum flight distance
  - `insecticide_resistance`: Resistance score (0-1)
  - `temp_preference_C`: Preferred temperature in Celsius
- **Rows**: 15 (one per sample)
- **Use**: Practice creating heatmaps and testing phylogenetic signal

## Data Description

### Samples Overview
- **Total samples**: 15
- **Genera**: 3 (Aedes, Culex, Anopheles)
- **Species**: 7 distinct species
- **Geographic range**: USA, Kenya, India
- **Collection years**: 2022-2023

### Genus Breakdown
- **Aedes**: 4 samples (aegypti, albopictus)
- **Culex**: 5 samples (pipiens, quinquefasciatus)
- **Anopheles**: 6 samples (gambiae, stephensi, quadrimaculatus, freeborni)

### Tree Properties
- **Total tips**: 15
- **Total branch length**: ~0.35 substitutions/site
- **Topology**: Rooted with Anopheles as outgroup
- **Bootstrap support**: 85-100% at most nodes

## Usage Examples

### Reading the Tree
```r
library(ape)
tree <- read.tree("example_tree.newick")
```

### Reading Metadata
```r
metadata <- read.csv("tree_metadata.csv")
```

### Joining with ggtree
```r
library(ggtree)
ggtree(tree) %<+% metadata +
  geom_tiplab(aes(color = genus))
```

### Reading Trait Data
```r
traits <- read.csv("trait_data.csv")
```

## Data Quality
- All sequences are from verified voucher specimens
- Bootstrap values indicate well-supported topology
- Metadata has been curated and validated
- Trait data are representative of published literature

## Citation
Data compiled for educational purposes from multiple sources including:
- Hoque et al. 2022 (mosquito COI sequences)
- Public GenBank records
- Literature values for morphological traits

## Notes
- Tip labels follow format: `Genus_species_Location_Year`
- All branch lengths are in substitutions per site
- Bootstrap values are percentages (0-100)
- Trait data are simulated but biologically realistic
- Some samples are laboratory colonies (indicated in metadata)

## File Formats

### Newick Format
```
((A:0.1,B:0.1)100:0.2,C:0.3);
```
- Parentheses group clades
- Colons precede branch lengths
- Numbers after parentheses are bootstrap values

### CSV Format
```
label,genus,species
Sample1,Aedes,aegypti
Sample2,Culex,pipiens
```
- First row is header
- Comma-separated columns
- No quotes needed unless values contain commas

## Troubleshooting

### Tree won't load
- Check file path is correct
- Use `file.exists("example_tree.newick")` to verify
- Make sure you're in the correct working directory

### Metadata doesn't join
- Verify tip labels in tree match `label` column in CSV
- Check for extra spaces or case differences
- Use `tree$tip.label` and `metadata$label` to compare

### Missing bootstrap values
- Some trees don't have bootstrap values - this is normal
- Check with `tree$node.label`
- Use bootstrap_tree.treefile if you need them

## Additional Resources
- [Newick format specification](http://evolution.genetics.washington.edu/phylip/newicktree.html)
- [ape package documentation](http://ape-package.ird.fr/)
- [ggtree book](https://yulab-smu.top/treedata-book/)
