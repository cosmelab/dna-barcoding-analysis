#!/usr/bin/env python3
"""
Extraction Comparison Plots - Lab Data Analysis Module

Creates two comparison plots:
1. HMW DNA extraction: Compare preservation methods (Frozen, Ethanol, Silica)
2. Column extraction: Compare yields for Aedes vs Culex species

Uses Plotly for interactive visualizations with professional Dracula-lite styling.
Includes statistical tests (ANOVA/t-test) with p-value annotations.
"""

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
import numpy as np
from pathlib import Path
import json
from data_loader import LabData


# Clean white style with Dracula accents
COLORS = {
    'purple': '#bd93f9',
    'cyan': '#8be9fd',
    'green': '#50fa7b',
    'orange': '#ffb86c',
    'red': '#ff5555',
    'background': 'white',
    'text': '#2d3436',
    'grid': '#e0e0e0'
}

# Color mapping for preservation methods
PRESERVATION_COLORS = {
    'Frozen': COLORS['cyan'],
    'Ethanol': COLORS['purple'],
    'Silica': COLORS['orange']
}

# Color mapping for species
SPECIES_COLORS = {
    'Aedes albopictus': COLORS['green'],
    'Culex spp.': COLORS['red']
}


def perform_anova(data, groups, values):
    """
    Perform one-way ANOVA to test differences between groups.

    Args:
        data: DataFrame with the data
        groups: Column name for groups
        values: Column name for values

    Returns:
        tuple: (F-statistic, p-value, group_data)
    """
    # Get data for each group
    unique_groups = data[groups].unique()
    group_data = [data[data[groups] == g][values].dropna() for g in unique_groups]

    # Perform ANOVA
    f_stat, p_value = stats.f_oneway(*group_data)

    return f_stat, p_value, group_data


def perform_ttest(group1, group2):
    """
    Perform independent t-test between two groups.

    Args:
        group1: Array-like data for group 1
        group2: Array-like data for group 2

    Returns:
        tuple: (t-statistic, p-value)
    """
    t_stat, p_value = stats.ttest_ind(group1, group2, equal_var=False)
    return t_stat, p_value


def format_pvalue(p):
    """Format p-value for display."""
    if p < 0.001:
        return "p < 0.001"
    elif p < 0.01:
        return f"p = {p:.3f}"
    else:
        return f"p = {p:.2f}"


def create_hmw_extraction_plot(data):
    """
    Create comparison plot for HMW DNA extraction across preservation methods.
    """
    print("Creating HMW DNA extraction comparison plot...")

    df = data.pcr_summary.copy()
    f_stat, p_value, group_data = perform_anova(df, 'Sample_Type', 'Mean')

    fig = go.Figure()
    preservation_methods = ['Frozen', 'Ethanol', 'Silica']

    for method in preservation_methods:
        method_data = df[df['Sample_Type'] == method]
        values = method_data['Mean'].values
        students = method_data['Student_Code'].values

        # Calculate stats for box hover
        q1 = np.percentile(values, 25)
        median_val = np.percentile(values, 50)
        q3 = np.percentile(values, 75)
        mean_val = np.mean(values)
        min_val = np.min(values)
        max_val = np.max(values)

        # Box - no hover (stats shown via annotation if needed)
        fig.add_trace(go.Box(
            y=values,
            name=method,
            marker_color=PRESERVATION_COLORS[method],
            boxmean='sd',
            showlegend=True,
            boxpoints=False,
            hoverinfo='skip'
        ))

        # Separate scatter for points with student:value hover
        fig.add_trace(go.Scatter(
            x=[method] * len(values),
            y=values,
            mode='markers',
            marker=dict(
                size=10,
                color=PRESERVATION_COLORS[method],
                line=dict(width=1, color='black')
            ),
            showlegend=False,
            hovertemplate="%{text}<extra></extra>",
            text=[f"{s}: {v:.2f}" for s, v in zip(students, values)]
        ))

    # Statistical annotation
    max_y = df['Mean'].max()
    fig.add_annotation(
        x=1, y=max_y * 1.15,
        text=f"ANOVA: F={f_stat:.2f}, {format_pvalue(p_value)}",
        showarrow=False,
        font=dict(size=12, color=COLORS['text']),
        bgcolor='white',
        bordercolor=COLORS['grid'],
        borderwidth=1,
        borderpad=4
    )

    fig.update_layout(
        title=dict(
            text='HMW DNA Extraction: Preservation Method Comparison',
            font=dict(size=18, color=COLORS['text']),
            x=0.5, xanchor='center'
        ),
        xaxis=dict(title='Preservation Method', tickfont=dict(size=11, color=COLORS['text']), gridcolor=COLORS['grid']),
        yaxis=dict(title='DNA Concentration (ng/uL)', tickfont=dict(size=11, color=COLORS['text']), gridcolor=COLORS['grid']),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(color=COLORS['text'], size=12),
        showlegend=True,
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='center', x=0.5, bgcolor='rgba(255,255,255,0.9)'),
        hovermode='closest',
        height=450,
        margin=dict(l=60, r=40, t=80, b=60)
    )

    # Create explanation text
    explanation = f"""
<h3>HMW DNA Extraction: Preservation Method Comparison</h3>

<p><strong>Research Question:</strong> Does preservation method affect DNA extraction yield?</p>

<p><strong>Methods:</strong> Students extracted high molecular weight (HMW) DNA from mosquito
samples preserved in three different ways: frozen (-80°C), 95% ethanol, and silica gel.
DNA concentrations were measured using NanoDrop spectrophotometry and averaged across
triplicate measurements.</p>

<p><strong>Statistical Analysis:</strong> One-way ANOVA was performed to test for significant
differences between preservation methods (F = {f_stat:.2f}, {format_pvalue(p_value)}).</p>

<p><strong>Results:</strong>
{'<span style="color: ' + COLORS['green'] + ';">Significant differences were detected</span>' if p_value < 0.05 else '<span style="color: ' + COLORS['orange'] + ';">No significant differences were detected</span>'}
between preservation methods (alpha = 0.05).
{'This suggests that preservation method has a measurable impact on DNA extraction efficiency.' if p_value < 0.05 else 'This suggests that all three preservation methods yield comparable DNA extraction results.'}
</p>

<p><strong>Key Observations:</strong></p>
<ul>
<li>Mean concentration for Frozen samples: {df[df['Sample_Type']=='Frozen']['Mean'].mean():.2f} ± {df[df['Sample_Type']=='Frozen']['Mean'].std():.2f} ng/uL</li>
<li>Mean concentration for Ethanol samples: {df[df['Sample_Type']=='Ethanol']['Mean'].mean():.2f} ± {df[df['Sample_Type']=='Ethanol']['Mean'].std():.2f} ng/uL</li>
<li>Mean concentration for Silica samples: {df[df['Sample_Type']=='Silica']['Mean'].mean():.2f} ± {df[df['Sample_Type']=='Silica']['Mean'].std():.2f} ng/uL</li>
</ul>

<p><strong>Biological Interpretation:</strong> Different preservation methods can affect DNA
integrity and extractability. Freezing typically maintains DNA quality well but requires
consistent cold chain. Ethanol and silica are field-friendly options that don't require
refrigeration, making them practical for remote sampling locations.</p>
"""

    print(f"  Statistical test: F = {f_stat:.2f}, {format_pvalue(p_value)}")

    return fig, explanation


def create_column_extraction_plot(data):
    """
    Create comparison plot for column extraction yields between species.
    """
    print("Creating column extraction comparison plot...")

    df = data.column_extraction.copy()
    aedes_data = df[df['Species'] == 'Aedes albopictus']['Total_DNA_ng']
    culex_data = df[df['Species'] == 'Culex spp.']['Total_DNA_ng']
    t_stat, p_value = perform_ttest(aedes_data, culex_data)

    fig = go.Figure()
    species_list = ['Aedes albopictus', 'Culex spp.']

    for species in species_list:
        species_data = df[df['Species'] == species]
        values = species_data['Total_DNA_ng'].values
        students = species_data['Student_Code'].values

        # Calculate stats for box hover
        q1 = np.percentile(values, 25)
        median_val = np.percentile(values, 50)
        q3 = np.percentile(values, 75)
        mean_val = np.mean(values)
        min_val = np.min(values)
        max_val = np.max(values)

        # Box - no hover (Plotly limitation: can't customize box hover)
        fig.add_trace(go.Box(
            y=values,
            name=species,
            marker_color=SPECIES_COLORS[species],
            boxmean='sd',
            showlegend=True,
            boxpoints=False,
            hoverinfo='skip'
        ))

        # Separate scatter for points with student:value hover
        fig.add_trace(go.Scatter(
            x=[species] * len(values),
            y=values,
            mode='markers',
            marker=dict(
                size=10,
                color=SPECIES_COLORS[species],
                line=dict(width=1, color='black')
            ),
            showlegend=False,
            hovertemplate="%{text}<extra></extra>",
            text=[f"{s}: {v:.1f}" for s, v in zip(students, values)]
        ))

    # Statistical annotation
    max_y = df['Total_DNA_ng'].max()
    fig.add_annotation(
        x=0.5, y=max_y * 1.10,
        text=f"t-test: t={t_stat:.2f}, {format_pvalue(p_value)}",
        showarrow=False,
        font=dict(size=12, color=COLORS['text']),
        bgcolor='white',
        bordercolor=COLORS['grid'],
        borderwidth=1,
        borderpad=4
    )

    fig.update_layout(
        title=dict(
            text='Column Extraction: Species Comparison',
            font=dict(size=18, color=COLORS['text']),
            x=0.5, xanchor='center'
        ),
        xaxis=dict(title='Species', tickfont=dict(size=11, color=COLORS['text']), gridcolor=COLORS['grid']),
        yaxis=dict(title='Total DNA Yield (ng)', tickfont=dict(size=11, color=COLORS['text']), gridcolor=COLORS['grid']),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(color=COLORS['text'], size=12),
        showlegend=True,
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='center', x=0.5, bgcolor='rgba(255,255,255,0.9)'),
        hovermode='closest',
        height=450,
        margin=dict(l=60, r=40, t=80, b=60)
    )

    # Create explanation text
    explanation = f"""
<h3>Column Extraction: Species Comparison</h3>

<p><strong>Research Question:</strong> Does mosquito species affect DNA extraction yield
using column-based purification?</p>

<p><strong>Methods:</strong> Students extracted DNA from two mosquito species
(<em>Aedes albopictus</em> and <em>Culex</em> spp.) using the Zymo Research Quick-DNA
Miniprep Plus Kit. This column-based method uses silica membrane technology for DNA
binding and purification. DNA was eluted in 30 uL of buffer and quantified using
NanoDrop spectrophotometry. Total DNA yield was calculated as:
Concentration (ng/uL) × Elution Volume (30 uL).</p>

<p><strong>Statistical Analysis:</strong> Independent samples t-test was performed to
compare yields between species (t = {t_stat:.2f}, {format_pvalue(p_value)}).</p>

<p><strong>Results:</strong>
{'<span style="color: ' + COLORS['green'] + ';">Significant differences were detected</span>' if p_value < 0.05 else '<span style="color: ' + COLORS['orange'] + ';">No significant differences were detected</span>'}
between species (alpha = 0.05).
{'This suggests that species identity affects DNA extraction efficiency with the column-based method.' if p_value < 0.05 else 'This suggests that both species yield comparable DNA amounts with the column-based method.'}
</p>

<p><strong>Key Observations:</strong></p>
<ul>
<li><em>Aedes albopictus</em> mean yield: {aedes_data.mean():.1f} ± {aedes_data.std():.1f} ng (n={len(aedes_data)})</li>
<li><em>Culex</em> spp. mean yield: {culex_data.mean():.1f} ± {culex_data.std():.1f} ng (n={len(culex_data)})</li>
<li>Percent difference: {abs(aedes_data.mean() - culex_data.mean()) / ((aedes_data.mean() + culex_data.mean())/2) * 100:.1f}%</li>
</ul>

<p><strong>Biological Interpretation:</strong> Differences in DNA yield between species
could be due to several factors: (1) differences in body size and biomass, (2) variation
in cell density or nuclear content, (3) differences in exoskeleton composition affecting
lysis efficiency, or (4) presence of PCR inhibitors that affect quantification accuracy.
<em>Aedes albopictus</em> (Asian tiger mosquito) is typically larger than many
<em>Culex</em> species, which could contribute to higher DNA yields.</p>

<p><strong>Methodological Note:</strong> The column-based extraction method is optimized
for speed and convenience, making it ideal for high-throughput applications like DNA
barcoding. However, yields are typically lower than traditional phenol-chloroform
extraction methods.</p>
"""

    print(f"  Statistical test: t = {t_stat:.2f}, {format_pvalue(p_value)}")

    return fig, explanation


def save_figure_json(fig, output_path, name):
    """
    Save Plotly figure as JSON for later embedding.

    Args:
        fig: Plotly figure object
        output_path: Path object for output directory
        name: Name for the JSON file (without extension)
    """
    json_file = output_path / f"{name}.json"

    with open(json_file, 'w') as f:
        json.dump(fig.to_dict(), f, indent=2)

    print(f"  Saved: {json_file}")


def main():
    """Main function to create both extraction comparison plots."""
    print("\n" + "="*70)
    print("EXTRACTION COMPARISON PLOTS")
    print("="*70 + "\n")

    # Load data
    data = LabData().load_all()

    # Create output directory
    output_dir = Path(__file__).parent.parent.parent / "results" / "lab_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}\n")

    # Create HMW extraction plot
    hmw_fig, hmw_explanation = create_hmw_extraction_plot(data)

    # Save HMW plot (JSON only - HTML combined in report)
    save_figure_json(hmw_fig, output_dir, "hmw_extraction_comparison")

    # Create column extraction plot
    column_fig, column_explanation = create_column_extraction_plot(data)

    # Save column plot (JSON only)
    save_figure_json(column_fig, output_dir, "column_extraction_comparison")

    # Save explanation texts
    explanations = {
        'hmw_extraction': hmw_explanation,
        'column_extraction': column_explanation
    }

    explanation_file = output_dir / "extraction_explanations.json"
    with open(explanation_file, 'w') as f:
        json.dump(explanations, f, indent=2)
    print(f"  Saved: {explanation_file}\n")

    print("="*70)
    print("EXTRACTION DATA SAVED")
    print("="*70 + "\n")

    return {
        'hmw_fig': hmw_fig,
        'column_fig': column_fig,
        'hmw_explanation': hmw_explanation,
        'column_explanation': column_explanation
    }


if __name__ == "__main__":
    main()
