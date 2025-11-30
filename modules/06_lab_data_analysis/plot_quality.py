#!/usr/bin/env python3
"""
DNA Quality Visualization - Lab Data Analysis Module

Creates scatter plots showing NanoDrop quality metrics with reference zones.
Split by extraction method (Column vs HMW/Magnetic Beads) in side-by-side facets.
Uses student codes only (HV, JR, MA, etc.) - no names.
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
from pathlib import Path
from data_loader import LabData


# Colors - white background style
COLORS = {
    'Excellent': '#50fa7b',  # Green
    'Good': '#50fa7b',       # Green
    'Fair': '#ffb86c',       # Orange
    'Poor': '#ff5555',       # Red
    'background': 'white',
    'text': '#2d3436',
    'grid': '#e0e0e0'
}

# Display names for extraction methods
METHOD_NAMES = {
    'Column': 'Column (Spin)',
    'Magnetic_Beads': 'HMW (Magnetic Beads)'
}


def create_quality_plot(data):
    """Create DNA quality scatter plot with reference zones, split by extraction method."""

    # Merge nanodrop with quality assessment
    df = data.nanodrop.merge(
        data.quality_assessment[['Student_Code', 'Extraction_Method', 'Sample_Type', 'Quality_Grade']],
        on=['Student_Code', 'Extraction_Method', 'Sample_Type'],
        how='left'
    )

    # Fill missing quality grades
    df['Quality_Grade'] = df['Quality_Grade'].fillna('Unknown')

    # Get extraction methods
    methods = ['Column', 'Magnetic_Beads']

    # Create subplots - side by side
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=[METHOD_NAMES.get(m, m) for m in methods],
        horizontal_spacing=0.12
    )

    # Track which grades we've added to legend
    legend_added = set()

    for col_idx, method in enumerate(methods, 1):
        method_df = df[df['Extraction_Method'] == method]

        # Add reference zones (ideal zone)
        # Ideal zone (260/280: 1.8-2.0, 260/230: >2.0)
        fig.add_shape(
            type="rect",
            x0=1.8, x1=2.0, y0=2.0, y1=4.0,
            fillcolor='rgba(80, 250, 123, 0.15)',
            line=dict(width=0),
            layer='below',
            row=1, col=col_idx
        )

        # Add reference lines
        # Ideal 260/280 range
        fig.add_vline(x=1.8, line_dash="dash", line_color='rgba(80, 250, 123, 0.5)',
                      line_width=2, row=1, col=col_idx)
        fig.add_vline(x=2.0, line_dash="dash", line_color='rgba(80, 250, 123, 0.5)',
                      line_width=2, row=1, col=col_idx)

        # Ideal 260/230 threshold
        fig.add_hline(y=2.0, line_dash="dash", line_color='rgba(80, 250, 123, 0.5)',
                      line_width=2, row=1, col=col_idx)

        # Add ideal zone annotation
        fig.add_annotation(
            x=1.9, y=3.0,
            text="Ideal",
            showarrow=False,
            font=dict(size=9, color='#50fa7b'),
            bgcolor='rgba(255,255,255,0.8)',
            row=1, col=col_idx
        )

        # Plot points by quality grade
        for grade in ['Excellent', 'Good', 'Fair', 'Poor', 'Unknown']:
            grade_df = method_df[method_df['Quality_Grade'] == grade]
            if len(grade_df) == 0:
                continue

            color = COLORS.get(grade, '#cccccc')

            # Create hover text
            hover_text = [
                f"{row['Student_Code']}: {row['Sample_Type']}<br>" +
                f"260/280: {row['Ratio_260_280']:.2f}<br>" +
                f"260/230: {row['Ratio_260_230']:.2f}"
                for _, row in grade_df.iterrows()
            ]

            # Only show legend for first occurrence
            show_legend = grade not in legend_added
            if show_legend:
                legend_added.add(grade)

            fig.add_trace(go.Scatter(
                x=grade_df['Ratio_260_280'],
                y=grade_df['Ratio_260_230'],
                mode='markers',
                name=grade,
                marker=dict(
                    size=12,
                    color=color,
                    line=dict(width=1, color='black')
                ),
                hovertemplate='%{text}<extra></extra>',
                text=hover_text,
                showlegend=show_legend,
                legendgroup=grade
            ), row=1, col=col_idx)

    # Update layout - clean white style
    fig.update_layout(
        title=dict(
            text="DNA Quality Assessment by Extraction Method",
            font=dict(size=18, color=COLORS['text']),
            x=0.5,
            xanchor='center'
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(color=COLORS['text'], size=12),
        hovermode='closest',
        legend=dict(
            title=dict(text="Quality"),
            orientation='h',
            yanchor='bottom',
            y=1.08,
            xanchor='center',
            x=0.5,
            bgcolor='rgba(255,255,255,0.9)',
        ),
        height=450,
        margin=dict(l=60, r=40, t=100, b=60),
    )

    # Update axes for both subplots
    for col_idx in [1, 2]:
        fig.update_xaxes(
            title_text="260/280 Ratio",
            tickfont=dict(size=11, color=COLORS['text']),
            gridcolor=COLORS['grid'],
            range=[-10, 20],
            zeroline=True,
            zerolinecolor=COLORS['grid'],
            row=1, col=col_idx
        )
        fig.update_yaxes(
            title_text="260/230 Ratio" if col_idx == 1 else "",
            tickfont=dict(size=11, color=COLORS['text']),
            gridcolor=COLORS['grid'],
            range=[-10, 5],
            zeroline=True,
            zerolinecolor=COLORS['grid'],
            row=1, col=col_idx
        )

    # Style subplot titles
    for annotation in fig['layout']['annotations']:
        if annotation['text'] in [METHOD_NAMES.get(m, m) for m in methods]:
            annotation['font'] = dict(size=14, color=COLORS['text'])

    return fig


def main():
    """Main function to create and save quality visualization."""
    print("DNA Quality Visualization")
    print("=" * 50)

    # Load data
    data = LabData().load_all()

    # Create plot
    print("\nCreating quality scatter plot (split by extraction method)...")
    fig = create_quality_plot(data)

    # Save as Plotly JSON only (no standalone HTML - will be embedded in report)
    output_dir = Path(__file__).parent.parent.parent / "results" / "lab_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "quality_plot.json"
    with open(output_file, 'w') as f:
        json.dump(fig.to_dict(), f, indent=2)

    print(f"Saved: {output_file}")

    # Print summary statistics
    print("\n" + "=" * 50)
    print("Quality Summary by Extraction Method:")
    print("=" * 50)

    # Merge for analysis
    df = data.nanodrop.merge(
        data.quality_assessment[['Student_Code', 'Extraction_Method', 'Sample_Type', 'Quality_Grade']],
        on=['Student_Code', 'Extraction_Method', 'Sample_Type'],
        how='left'
    )

    for method in ['Column', 'Magnetic_Beads']:
        method_df = df[df['Extraction_Method'] == method]
        print(f"\n{METHOD_NAMES.get(method, method)}:")
        quality_counts = method_df['Quality_Grade'].value_counts()
        for grade, count in quality_counts.items():
            print(f"  {grade}: {count} samples")

        # Samples in ideal range
        ideal = method_df[
            (method_df['Ratio_260_280'] >= 1.8) &
            (method_df['Ratio_260_280'] <= 2.0) &
            (method_df['Ratio_260_230'] >= 2.0)
        ]
        print(f"  In ideal range: {len(ideal)} / {len(method_df)}")

    print("\nDone!")


if __name__ == "__main__":
    main()
