# Lab Data Analysis Module

**Statistical analysis of DNA barcoding lab data**

This module analyzes class lab data to identify factors affecting DNA extraction quality, PCR success, and sequencing outcomes.

---

## Features

### Statistical Analyses
- **Extraction Methods**: Compare HMW (magnetic beads) vs Column extraction
  - Welch's t-test for mean comparison
  - Box plots with confidence intervals

- **DNA Quality**: NanoDrop purity assessment
  - 260/280 ratio (protein contamination)
  - 260/230 ratio (organic contamination)
  - Quality zones visualization

- **PCR Success**: Correlation with DNA quality
  - One-way ANOVA by quality grade
  - Bar charts with error bars

- **Sequencing Success**: QC pass rates and consensus generation
  - Success metrics by student
  - Species identification results

- **Summary Heatmap**: Overall quality scores across all samples

### Visualizations
- Interactive plotly plots (zoom, pan, hover tooltips)
- Statistical annotations (p-values, significance)
- Color-coded quality grades
- Responsive HTML report with tabs

---

## Usage

### Running the Analysis

```bash
# From the module directory
python3 analyze_lab_data.py

# Or from the repository root
python3 modules/06_lab_data_analysis/analyze_lab_data.py
```

### Input Data

The module automatically loads data from `../../data/class_data/`:
- `raw_measurements.csv` - HMW-DNA extraction measurements
- `column_extraction.csv` - Column-based extraction results
- `pcr_products_raw.csv` - PCR product concentrations
- `pcr_products_summary.csv` - PCR statistics by student
- `nanodrop_quality.csv` - DNA purity ratios
- `extraction_quality_assessment.csv` - Quality grades
- `coi_pcr_gel_results.csv` - Gel electrophoresis results
- `wolbachia_pcr.csv` - Wolbachia detection
- `sequencing_success_summary.csv` - Sanger sequencing QC

### Output Files

Results are saved to `../../results/lab_analysis/`:

**Plots (HTML)**:
- `plots/extraction_comparison.html` - Extraction method comparison
- `plots/quality_scatter.html` - DNA quality scatter plot
- `plots/pcr_success.html` - PCR yields by student
- `plots/sequencing_success.html` - Sequencing success metrics
- `plots/quality_heatmap.html` - Quality score heatmap

**Statistics (JSON)**:
- `statistics/statistical_tests.json` - All test results with p-values

**Report**:
- `lab_analysis_report.html` - Complete interactive report

---

## Running in Docker

```bash
# Start container with mounted data
docker run -it --rm \
  -v "$(pwd)/data:/workspace/data" \
  -v "$(pwd)/results:/workspace/results" \
  -v "$(pwd)/modules:/workspace/modules" \
  cosmelab/dna-barcoding-analysis:latest

# Inside container
cd modules/06_lab_data_analysis
python3 analyze_lab_data.py

# View report
# results/lab_analysis/lab_analysis_report.html
```

---

## Requirements

**Python Packages** (all included in container):
- pandas - Data manipulation
- numpy - Numerical computing
- plotly - Interactive visualizations
- seaborn - Statistical plotting
- matplotlib - Base plotting
- scipy - Statistical tests
- statsmodels - Statistical modeling

---

## Statistical Methods

### Extraction Comparison
- **Test**: Welch's t-test (unequal variances)
- **Null Hypothesis**: No difference in mean DNA concentration
- **Significance**: α = 0.05

### PCR Success by Quality
- **Test**: One-way ANOVA
- **Groups**: Excellent, Good, Fair, Poor
- **Null Hypothesis**: No difference in mean PCR yield
- **Significance**: α = 0.05

### Quality Zones
- **Excellent**: 260/280 = 1.8-2.0 AND 260/230 ≥ 2.0
- **Good**: One ratio in ideal range
- **Fair**: Borderline contamination
- **Poor**: Both ratios outside ideal range

---

## Interpreting Results

### Extraction Methods
Compare yields between magnetic bead and column methods:
- **p < 0.05**: Significant difference in extraction efficiency
- Look for consistency across sample types

### DNA Quality
Check if samples fall in green quality zone:
- **260/280 < 1.8**: Protein contamination
- **260/280 > 2.0**: RNA or degraded DNA
- **260/230 < 2.0**: Organic contamination

### PCR Success
Examine correlation between quality grade and PCR yield:
- **p < 0.05**: DNA quality significantly affects PCR
- Higher quality DNA should produce higher PCR yields

### Sequencing Success
Identify factors contributing to QC failures:
- Compare QC pass rates across students
- Check if quality/PCR metrics predict sequencing success

---

## Example Questions to Explore

1. Does extraction method affect DNA concentration?
2. Which students had the best DNA quality?
3. Does DNA quality predict PCR success?
4. What factors caused sequencing failures?
5. Which sample types worked best?
6. How do teams compare (Spin vs Magnet)?

---

## Troubleshooting

**Module not found errors:**
```bash
# Ensure you're in the correct directory
pwd  # Should show path to dna-barcoding-analysis
```

**Missing data files:**
```bash
# Check data directory exists
ls -la data/class_data/
# Should show 9 CSV files + gels/ directory
```

**Import errors (missing packages):**
```bash
# Run inside Docker container
# All required packages are pre-installed
```

---

## Customization

### Adding New Analyses

Edit `analyze_lab_data.py` to add custom analyses:

```python
def my_custom_analysis(self):
    """Your custom analysis."""
    data = self.data['dataset_name']

    # Your analysis code here

    # Create visualization
    fig = px.scatter(...)
    fig.write_html(self.output_dir / "plots" / "my_plot.html")
```

### Modifying Plot Styles

Change colors, themes, or layouts:

```python
# Change color scheme
color_discrete_map={
    'Excellent': '#your_color',
    'Good': '#your_color',
    ...
}

# Change template
fig.update_layout(template='plotly_dark')  # or seaborn, simple, ggplot2
```

---

## Contact

Questions about the analysis module?
- Check the tracking system: `tracking/log.md`
- Review the code comments in `analyze_lab_data.py`

---

**Generated by**: DNA Barcoding Analysis Pipeline
**Last Updated**: 2025-11-26
