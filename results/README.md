# Results Directory

This directory is where your analysis outputs are saved.

## Structure

When you run the analysis, you'll get numbered subdirectories:

```
results/
├── tutorial/           # Tutorial results (test data)
│   ├── 01_qc/
│   ├── 02_consensus/
│   ├── 03_alignment/
│   ├── 04_phylogeny/
│   └── 05_blast/
│
└── my_analysis/        # Your real analysis results
    ├── 01_qc/
    ├── 02_consensus/
    ├── 03_alignment/
    ├── 04_phylogeny/
    └── 05_blast/
```

## What Creates These?

- **tutorial/** - Created when you run `./tutorial.sh`
- **my_analysis/** - Created when you run `./run-analysis.sh`

The numbered subdirectories (01, 02, 03, 04, 05) match the 5 steps of the workflow.

## Important Files to Check

After running your analysis, look at:

- `my_analysis/01_qc/qc_report.html` - Which sequences passed QC?
- `my_analysis/04_phylogeny/tree.png` - Phylogenetic tree
- `my_analysis/05_blast/identification_report.html` - Species identification
