#!/usr/bin/env python3
"""
Sequencing Results Visualization - Lab Data Analysis Module

Creates interactive visualizations of DNA sequencing success rates and
quality control metrics for student lab submissions.

Generates:
- Bar chart showing QC pass rates by student
- Visual breakdown of sequencing pipeline success
- Educational annotations explaining QC metrics
"""

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
import json
from data_loader import LabData


# Clean white style with Dracula accents
COLORS = {
    'high_pass': '#50fa7b',      # Green - high success (>70%)
    'medium_pass': '#f1fa8c',    # Yellow - medium success (40-70%)
    'low_pass': '#ff5555',       # Red - low success (<40%)
    'total': '#bd93f9',          # Purple - total sequences
    'passed': '#50fa7b',         # Green - passed QC
    'failed': '#ff5555',         # Red - failed QC
    'consensus': '#8be9fd',      # Cyan - consensus generated
    'identified': '#ffb86c',     # Orange - species identified
    'background': 'white',       # White background
    'text': '#2d3436',           # Dark text
    'grid': '#e0e0e0'            # Light grid
}


def get_pass_rate_color(pass_rate):
    """
    Return color based on QC pass rate threshold.

    Args:
        pass_rate: Pass rate percentage (0-100)

    Returns:
        Color hex code
    """
    if pass_rate >= 70:
        return COLORS['high_pass']
    elif pass_rate >= 40:
        return COLORS['medium_pass']
    else:
        return COLORS['low_pass']


def create_pass_rate_chart(df):
    """
    Create bar chart showing QC pass rate by student.

    Args:
        df: DataFrame with sequencing data

    Returns:
        Plotly figure object
    """
    # Sort by pass rate descending
    df_sorted = df.sort_values('QC_Pass_Rate_%', ascending=True)

    # Assign colors based on pass rate
    colors = [get_pass_rate_color(rate) for rate in df_sorted['QC_Pass_Rate_%']]

    # Create horizontal bar chart
    fig = go.Figure()

    fig.add_trace(go.Bar(
        x=df_sorted['QC_Pass_Rate_%'],
        y=df_sorted['Student_Code'],
        orientation='h',
        marker=dict(
            color=colors,
            line=dict(color=COLORS['grid'], width=1)
        ),
        text=df_sorted['QC_Pass_Rate_%'].apply(lambda x: f'{x:.0f}%'),
        textposition='outside',
        hovertemplate='%{x:.0f}%<extra></extra>'
    ))

    fig.update_layout(
        title={
            'text': 'DNA Sequencing QC Pass Rate by Student',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'color': COLORS['text']}
        },
        xaxis=dict(
            title='QC Pass Rate (%)',
            range=[0, 105],
            gridcolor=COLORS['grid'],
            tickfont=dict(size=11, color=COLORS['text'])
        ),
        yaxis=dict(
            title='Student Code',
            gridcolor=COLORS['grid'],
            tickfont=dict(size=11, color=COLORS['text'])
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(color=COLORS['text'], size=12),
        height=450,
        margin=dict(l=60, r=40, t=80, b=60)
    )

    return fig


def create_sequencing_breakdown(df):
    """
    Create stacked bar chart showing sequencing pipeline breakdown.

    Shows:
    - Total sequences submitted
    - Passed QC count
    - Consensus generated count
    - Species identified count

    Args:
        df: DataFrame with sequencing data

    Returns:
        Plotly figure object
    """
    # Sort by total sequences descending
    df_sorted = df.sort_values('Total_Sequences', ascending=True)

    fig = go.Figure()

    # Total sequences (as background)
    fig.add_trace(go.Bar(
        y=df_sorted['Student_Code'],
        x=df_sorted['Total_Sequences'],
        name='Total Submitted',
        orientation='h',
        marker=dict(color=COLORS['total'], opacity=0.3),
        hovertemplate='%{x}<extra></extra>'
    ))

    # Passed QC
    fig.add_trace(go.Bar(
        y=df_sorted['Student_Code'],
        x=df_sorted['Passed_QC'],
        name='Passed QC',
        orientation='h',
        marker=dict(color=COLORS['passed']),
        hovertemplate='%{x}<extra></extra>'
    ))

    # Consensus generated
    fig.add_trace(go.Bar(
        y=df_sorted['Student_Code'],
        x=df_sorted['Consensus_Generated'],
        name='Consensus Generated',
        orientation='h',
        marker=dict(color=COLORS['consensus']),
        hovertemplate='%{x}<extra></extra>'
    ))

    # Species identified (count non-None entries)
    species_count = df_sorted['Species_Identified'].apply(
        lambda x: 0 if x == 'None' else len(str(x).split(';'))
    )

    fig.add_trace(go.Bar(
        y=df_sorted['Student_Code'],
        x=species_count,
        name='Species Identified',
        orientation='h',
        marker=dict(color=COLORS['identified']),
        hovertemplate='%{x}<extra></extra>'
    ))

    fig.update_layout(
        title={
            'text': 'Sequencing Pipeline Success by Stage',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'color': COLORS['text']}
        },
        xaxis=dict(
            title='Number of Sequences',
            gridcolor=COLORS['grid'],
            tickfont=dict(size=11, color=COLORS['text'])
        ),
        yaxis=dict(
            title='Student Code',
            gridcolor=COLORS['grid'],
            tickfont=dict(size=11, color=COLORS['text'])
        ),
        barmode='overlay',
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(color=COLORS['text'], size=12),
        height=450,
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=1.02,
            xanchor='center',
            x=0.5,
            bgcolor='rgba(255,255,255,0.9)'
        ),
        margin=dict(l=60, r=40, t=80, b=60)
    )

    return fig


def create_educational_summary(df):
    """
    Create text summary with educational context.

    Args:
        df: DataFrame with sequencing data

    Returns:
        HTML formatted string with educational content
    """
    total_sequences = df['Total_Sequences'].sum()
    total_passed = df['Passed_QC'].sum()
    total_failed = df['Failed_QC'].sum()
    total_consensus = df['Consensus_Generated'].sum()
    overall_pass_rate = (total_passed / total_sequences * 100) if total_sequences > 0 else 0

    # Count species identified
    species_identified = 0
    for species in df['Species_Identified']:
        if species != 'None':
            species_identified += len(str(species).split(';'))

    summary = f"""
    <div style="font-family: monospace; color: {COLORS['text']}; background-color: {COLORS['background']}; padding: 20px; border-radius: 5px;">
        <h3 style="color: {COLORS['high_pass']};">Understanding DNA Sequencing Quality Control</h3>

        <h4 style="color: {COLORS['text']};">Class Summary:</h4>
        <ul>
            <li><strong>Total Sequences Submitted:</strong> {total_sequences}</li>
            <li><strong>Passed QC:</strong> {total_passed} ({overall_pass_rate:.1f}%)</li>
            <li><strong>Failed QC:</strong> {total_failed}</li>
            <li><strong>Consensus Sequences Generated:</strong> {total_consensus}</li>
            <li><strong>Species Successfully Identified:</strong> {species_identified}</li>
        </ul>

        <h4 style="color: {COLORS['high_pass']};">What is Quality Control (QC)?</h4>
        <p>
        Quality Control in DNA sequencing evaluates the reliability of chromatogram data.
        Each sequencing reaction produces a chromatogram showing peaks for A, T, G, and C bases.
        QC software analyzes:
        </p>
        <ul>
            <li><strong>Peak height:</strong> Strong, clear peaks indicate good signal</li>
            <li><strong>Peak resolution:</strong> Well-separated peaks show clean reads</li>
            <li><strong>Baseline noise:</strong> Low background noise improves accuracy</li>
            <li><strong>Quality scores:</strong> Phred scores quantify base-calling confidence</li>
        </ul>

        <h4 style="color: {COLORS['low_pass']};">Why Do Sequences Fail QC?</h4>
        <ul>
            <li><strong>Poor DNA quality:</strong> Degraded or contaminated template DNA</li>
            <li><strong>Low concentration:</strong> Insufficient DNA for strong signal</li>
            <li><strong>Mixed templates:</strong> Multiple PCR products create overlapping peaks</li>
            <li><strong>Primer issues:</strong> Non-specific binding or primer dimers</li>
            <li><strong>Technical problems:</strong> Instrument errors or reaction failure</li>
        </ul>

        <h4 style="color: {COLORS['consensus']};">What is a Consensus Sequence?</h4>
        <p>
        A consensus sequence is generated by aligning forward and reverse sequencing reads
        of the same PCR product. The software:
        </p>
        <ol>
            <li>Aligns the complementary reads</li>
            <li>Compares base calls at each position</li>
            <li>Uses quality scores to resolve disagreements</li>
            <li>Produces a single, high-confidence sequence</li>
        </ol>
        <p>
        Consensus sequences are more accurate than individual reads because errors are
        unlikely to occur at the same position in both forward and reverse reads.
        </p>

        <h4 style="color: {COLORS['identified']};">Species Identification:</h4>
        <p>
        After generating consensus sequences, they are compared to reference databases
        (like BOLD or GenBank) using BLAST. Matches above a similarity threshold
        (typically 97-99% for COI) allow species-level identification.
        </p>
    </div>
    """

    return summary


def main():
    """
    Main function to generate sequencing results visualizations.
    """
    print("\n" + "="*70)
    print("DNA SEQUENCING RESULTS VISUALIZATION")
    print("="*70 + "\n")

    # Load data
    print("Loading sequencing data...")
    data = LabData()
    data.load_all()

    if data.sequencing is None:
        print("ERROR: Sequencing data not found!")
        return

    df = data.sequencing
    print(f"Loaded data for {len(df)} students\n")

    # Create output directory
    output_dir = Path(__file__).parent.parent.parent / "results" / "lab_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate visualizations
    print("Generating visualizations...")

    # 1. QC Pass Rate Chart
    print("  - QC pass rate chart")
    fig_pass_rate = create_pass_rate_chart(df)
    pass_rate_file = output_dir / "sequencing_pass_rate.json"
    fig_pass_rate.write_json(str(pass_rate_file))
    print(f"    Saved: {pass_rate_file}")

    # 2. Sequencing Breakdown Chart
    print("  - Pipeline breakdown chart")
    fig_breakdown = create_sequencing_breakdown(df)
    breakdown_file = output_dir / "sequencing_breakdown.json"
    fig_breakdown.write_json(str(breakdown_file))
    print(f"    Saved: {breakdown_file}")

    # 3. Educational Summary
    print("  - Educational summary")
    summary_html = create_educational_summary(df)
    summary_file = output_dir / "sequencing_summary.html"
    with open(summary_file, 'w') as f:
        f.write(summary_html)
    print(f"    Saved: {summary_file}")

    # Print summary statistics
    print("\n" + "-"*70)
    print("SUMMARY STATISTICS")
    print("-"*70)
    print(f"{'Student':<10} {'Total':<8} {'Passed':<8} {'Failed':<8} {'Rate':<10} {'Consensus':<10} {'Species'}")
    print("-"*70)

    for _, row in df.sort_values('QC_Pass_Rate_%', ascending=False).iterrows():
        species = row['Species_Identified']
        # Handle NaN, None, or 'None' string
        if pd.isna(species) or species == 'None':
            species_display = 'None'
        else:
            # Convert to string and truncate long species names
            species_str = str(species)
            species_display = species_str[:30] + '...' if len(species_str) > 30 else species_str

        print(f"{row['Student_Code']:<10} "
              f"{row['Total_Sequences']:<8} "
              f"{row['Passed_QC']:<8} "
              f"{row['Failed_QC']:<8} "
              f"{row['QC_Pass_Rate_%']:<10.0f} "
              f"{row['Consensus_Generated']:<10} "
              f"{species_display}")

    print("-"*70)
    print(f"\nTotal sequences: {df['Total_Sequences'].sum()}")
    print(f"Overall QC pass rate: {(df['Passed_QC'].sum() / df['Total_Sequences'].sum() * 100):.1f}%")
    print(f"Total consensus sequences: {df['Consensus_Generated'].sum()}")

    print("\n" + "="*70)
    print("VISUALIZATION COMPLETE")
    print("="*70 + "\n")
    print(f"Output directory: {output_dir}")
    print("\nFiles generated:")
    print(f"  - {pass_rate_file.name}")
    print(f"  - {breakdown_file.name}")
    print(f"  - {summary_file.name}")
    print()


if __name__ == "__main__":
    main()
