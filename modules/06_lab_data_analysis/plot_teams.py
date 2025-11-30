#!/usr/bin/env python3
"""
Team Comparison - Gamification Plot
Compares Spin vs Magnet teams across key lab metrics.

Shows side-by-side team performance with competitive scoreboard.
Uses student codes only (no names).
"""

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
from pathlib import Path
from data_loader import LabData


# Clean white style with Dracula accents
COLORS = {
    'Spin': '#bd93f9',    # Purple
    'Magnet': '#f1c40f',  # Golden
    'text': '#2d3436',    # Dark text
    'grid': '#e0e0e0'     # Light grid
}


def calculate_team_metrics(data):
    """Calculate team performance metrics."""

    metrics = {}

    # Metric 1: DNA Extraction Yield (average total DNA)
    # column_extraction already has Team column
    if data.column_extraction is not None:
        df = data.column_extraction.copy()
        team_yields = df.groupby('Team')['Total_DNA_ng'].mean()
        metrics['DNA Yield (ng)'] = {
            'Spin': team_yields.get('Spin', 0),
            'Magnet': team_yields.get('Magnet', 0),
            'higher_better': True
        }

    # Metric 2: DNA Quality Score (% samples with good 260/280 ratio)
    # nanodrop does NOT have Team column - need to merge
    if data.nanodrop is not None and data.teams is not None:
        df = data.nanodrop.copy()
        # Add team info
        df = df.merge(data.teams, on='Student_Code', how='left')

        # Filter for column extraction only (main DNA samples)
        df = df[df['Extraction_Method'] == 'Column']

        # Good quality: 260/280 ratio between 1.8 and 2.0
        df['quality_good'] = ((df['Ratio_260_280'] >= 1.8) &
                               (df['Ratio_260_280'] <= 2.0)).astype(int)

        team_quality = df.groupby('Team')['quality_good'].mean() * 100
        metrics['DNA Quality (%)'] = {
            'Spin': team_quality.get('Spin', 0),
            'Magnet': team_quality.get('Magnet', 0),
            'higher_better': True
        }

    # Metric 3: PCR Success Rate
    # gel_results already has Team column
    if data.gel_results is not None:
        df = data.gel_results.copy()
        # Exclude positive controls
        df = df[df['Sample_Type'] != 'control+']

        team_pcr = df.groupby('Team')['PCR_Success'].mean() * 100
        metrics['PCR Success (%)'] = {
            'Spin': team_pcr.get('Spin', 0),
            'Magnet': team_pcr.get('Magnet', 0),
            'higher_better': True
        }

    # Metric 4: Sequencing QC Pass Rate
    # sequencing does NOT have Team column - need to merge
    if data.sequencing is not None and data.teams is not None:
        df = data.sequencing.copy()
        # Add team info
        df = df.merge(data.teams, on='Student_Code', how='left')

        team_seq = df.groupby('Team')['QC_Pass_Rate_%'].mean()
        metrics['Sequencing QC (%)'] = {
            'Spin': team_seq.get('Spin', 0),
            'Magnet': team_seq.get('Magnet', 0),
            'higher_better': True
        }

    return metrics


def create_team_comparison(data):
    """Create team comparison visualization with scoreboard."""

    metrics = calculate_team_metrics(data)

    if not metrics:
        print("No metrics available for team comparison")
        return None

    # Calculate scoreboard
    spin_wins = 0
    magnet_wins = 0

    for metric_name, values in metrics.items():
        spin_val = values['Spin']
        magnet_val = values['Magnet']
        if spin_val > magnet_val:
            spin_wins += 1
        elif magnet_val > spin_val:
            magnet_wins += 1

    # Create simple bar chart (no subplots to avoid table issues)
    fig = go.Figure()

    # Prepare data for main comparison
    metric_names = list(metrics.keys())
    spin_values = [metrics[m]['Spin'] for m in metric_names]
    magnet_values = [metrics[m]['Magnet'] for m in metric_names]

    # Add bars for each team
    fig.add_trace(
        go.Bar(
            name='Spin Team',
            x=metric_names,
            y=spin_values,
            marker_color=COLORS['Spin'],
            text=[f"{v:.1f}" for v in spin_values],
            textposition='outside',
            hovertemplate='%{y:.1f}<extra></extra>'
        )
    )

    fig.add_trace(
        go.Bar(
            name='Magnet Team',
            x=metric_names,
            y=magnet_values,
            marker_color=COLORS['Magnet'],
            text=[f"{v:.1f}" for v in magnet_values],
            textposition='outside',
            hovertemplate='%{y:.1f}<extra></extra>'
        )
    )

    # Determine winner for title
    if spin_wins > magnet_wins:
        winner_text = f"Spin Team leads {spin_wins}-{magnet_wins}!"
    elif magnet_wins > spin_wins:
        winner_text = f"Magnet Team leads {magnet_wins}-{spin_wins}!"
    else:
        winner_text = f"Teams tied {spin_wins}-{magnet_wins}!"

    fig.update_layout(
        title={
            'text': f'Spin vs Magnet: Lab Performance Challenge — {winner_text}',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 16, 'color': COLORS['text']}
        },
        xaxis=dict(
            title='Metrics',
            gridcolor=COLORS['grid'],
            tickfont=dict(size=11, color=COLORS['text'])
        ),
        yaxis=dict(
            title='Value',
            gridcolor=COLORS['grid'],
            tickfont=dict(size=11, color=COLORS['text'])
        ),
        barmode='group',
        height=450,
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.15,
            xanchor="center",
            x=0.5,
            bgcolor='rgba(255,255,255,0.9)'
        ),
        font=dict(size=12, color=COLORS['text']),
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=60, r=40, t=60, b=80)
    )

    return fig


def save_plot(fig, output_dir):
    """Save plot as Plotly JSON."""
    output_path = Path(output_dir) / "team_comparison.json"

    with open(output_path, 'w') as f:
        json.dump(fig.to_dict(), f, indent=2)

    print(f"Saved: {output_path}")
    return output_path


def main():
    """Main execution function."""

    print("=" * 60)
    print("TEAM COMPARISON - Spin vs Magnet")
    print("=" * 60)
    print()

    # Load data
    data = LabData().load_all()

    # Verify team assignments
    print("Team Rosters:")
    spin_students = data.get_students_by_team('Spin')
    magnet_students = data.get_students_by_team('Magnet')
    print(f"  Spin Team: {', '.join(spin_students)}")
    print(f"  Magnet Team: {', '.join(magnet_students)}")
    print()

    # Create comparison plot
    print("Calculating team metrics...")
    fig = create_team_comparison(data)

    if fig is None:
        print("Could not create team comparison (insufficient data)")
        return

    # Save output
    output_dir = Path(__file__).parent.parent.parent / "results" / "lab_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    save_plot(fig, output_dir)

    print()
    print("Team comparison complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
