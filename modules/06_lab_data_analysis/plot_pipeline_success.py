#!/usr/bin/env python3
"""
Pipeline Success Correlation - Lab Data Analysis Module

Shows how DNA quality correlates with downstream success:
Quality → PCR Success → Sequencing QC → Species ID

Demonstrates that upstream quality determines downstream success.
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import json
from pathlib import Path
from data_loader import LabData


# Colors
COLORS = {
    'success': '#50fa7b',    # Green
    'partial': '#ffb86c',    # Orange
    'fail': '#ff5555',       # Red
    'neutral': '#8be9fd',    # Cyan
    'purple': '#bd93f9',
    'background': 'white',
    'text': '#2d3436',
    'grid': '#e0e0e0'
}

# Quality grade colors
QUALITY_COLORS = {
    'Excellent': '#50fa7b',
    'Good': '#50fa7b',
    'Fair': '#ffb86c',
    'Poor': '#ff5555'
}


def create_pipeline_funnel(data):
    """Create a funnel/Sankey showing pipeline success per student."""

    # Build student summary
    students = []

    # Get unique students
    all_students = set(data.quality_assessment['Student_Code'].unique())

    for student in sorted(all_students):
        # Quality: get predominant grade
        student_quality = data.quality_assessment[
            data.quality_assessment['Student_Code'] == student
        ]['Quality_Grade'].value_counts()

        if len(student_quality) > 0:
            # Calculate quality score (Poor=1, Fair=2, Good=3, Excellent=4)
            grade_scores = {'Poor': 1, 'Fair': 2, 'Good': 3, 'Excellent': 4}
            quality_data = data.quality_assessment[
                data.quality_assessment['Student_Code'] == student
            ]['Quality_Grade']
            avg_quality = quality_data.map(grade_scores).mean()
            best_quality = student_quality.index[0]
        else:
            avg_quality = 0
            best_quality = 'Unknown'

        # PCR success
        student_pcr = data.gel_results[
            data.gel_results['Student_Code'] == student
        ]
        pcr_total = len(student_pcr[student_pcr['Sample_Type'] != 'control+'])
        pcr_success = student_pcr[
            (student_pcr['PCR_Success'] == 1) &
            (student_pcr['Sample_Type'] != 'control+')
        ].shape[0]
        pcr_rate = pcr_success / pcr_total * 100 if pcr_total > 0 else 0

        # Sequencing
        student_seq = data.sequencing[
            data.sequencing['Student_Code'] == student
        ]
        if len(student_seq) > 0:
            seq_total = student_seq['Total_Sequences'].values[0]
            seq_passed = student_seq['Passed_QC'].values[0]
            consensus = student_seq['Consensus_Generated'].values[0]
            species = student_seq['Species_Identified'].values[0]
            has_species = pd.notna(species) and species != ''
        else:
            seq_total = 0
            seq_passed = 0
            consensus = 0
            has_species = False
            species = None

        students.append({
            'Student': student,
            'Quality_Score': avg_quality,
            'Best_Quality': best_quality,
            'PCR_Rate': pcr_rate,
            'PCR_Success': pcr_success,
            'PCR_Total': pcr_total,
            'Seq_Total': seq_total,
            'Seq_Passed': seq_passed,
            'Consensus': consensus,
            'Has_Species': has_species,
            'Species': species
        })

    df = pd.DataFrame(students)

    # Sort by pipeline success (species ID first, then consensus, then PCR)
    df['Success_Score'] = (
        df['Has_Species'].astype(int) * 100 +
        df['Consensus'] * 10 +
        df['PCR_Rate']
    )
    df = df.sort_values('Success_Score', ascending=True)

    # Create horizontal bar chart showing pipeline stages
    fig = make_subplots(
        rows=1, cols=4,
        subplot_titles=['DNA Quality', 'PCR Success', 'Seq QC Pass', 'Species ID'],
        horizontal_spacing=0.08,
        shared_yaxes=True
    )

    # Column 1: Quality Score (1-4 scale)
    quality_colors = [QUALITY_COLORS.get(q, '#cccccc') for q in df['Best_Quality']]
    fig.add_trace(go.Bar(
        y=df['Student'],
        x=df['Quality_Score'],
        orientation='h',
        marker_color=quality_colors,
        text=[f"{q}" for q in df['Best_Quality']],
        textposition='inside',
        hovertemplate='%{y}: %{text}<extra></extra>',
        showlegend=False
    ), row=1, col=1)

    # Column 2: PCR Success Rate
    pcr_colors = [
        COLORS['success'] if r >= 66 else COLORS['partial'] if r >= 33 else COLORS['fail']
        for r in df['PCR_Rate']
    ]
    fig.add_trace(go.Bar(
        y=df['Student'],
        x=df['PCR_Rate'],
        orientation='h',
        marker_color=pcr_colors,
        text=[f"{int(r)}%" for r in df['PCR_Rate']],
        textposition='inside',
        hovertemplate='%{y}: %{x:.0f}% (%{customdata})<extra></extra>',
        customdata=[f"{s}/{t}" for s, t in zip(df['PCR_Success'], df['PCR_Total'])],
        showlegend=False
    ), row=1, col=2)

    # Column 3: Sequencing QC Pass
    seq_rates = []
    for _, row in df.iterrows():
        if row['Seq_Total'] > 0:
            seq_rates.append(row['Seq_Passed'] / row['Seq_Total'] * 100)
        else:
            seq_rates.append(0)

    seq_colors = [
        COLORS['success'] if r >= 66 else COLORS['partial'] if r >= 33 else COLORS['fail']
        for r in seq_rates
    ]
    fig.add_trace(go.Bar(
        y=df['Student'],
        x=seq_rates,
        orientation='h',
        marker_color=seq_colors,
        text=[f"{int(r)}%" if r > 0 else "N/A" for r in seq_rates],
        textposition='inside',
        hovertemplate='%{y}: %{customdata}<extra></extra>',
        customdata=[f"{p}/{t}" if t > 0 else "Not sequenced"
                   for p, t in zip(df['Seq_Passed'], df['Seq_Total'])],
        showlegend=False
    ), row=1, col=3)

    # Column 4: Species Identified (binary)
    species_colors = [COLORS['success'] if s else COLORS['fail'] for s in df['Has_Species']]
    fig.add_trace(go.Bar(
        y=df['Student'],
        x=[100 if s else 0 for s in df['Has_Species']],
        orientation='h',
        marker_color=species_colors,
        text=['Yes' if s else 'No' for s in df['Has_Species']],
        textposition='inside',
        hovertemplate='%{y}: %{customdata}<extra></extra>',
        customdata=[s if pd.notna(s) else "No ID" for s in df['Species']],
        showlegend=False
    ), row=1, col=4)

    # Update layout
    fig.update_layout(
        title=dict(
            text="Pipeline Success: From DNA Quality to Species ID",
            font=dict(size=18, color=COLORS['text']),
            x=0.5,
            xanchor='center'
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(color=COLORS['text'], size=12),
        height=450,
        margin=dict(l=60, r=40, t=100, b=60),
        bargap=0.3
    )

    # Update x-axes
    fig.update_xaxes(range=[0, 4.5], title_text="Score (1-4)", row=1, col=1,
                     gridcolor=COLORS['grid'])
    fig.update_xaxes(range=[0, 105], title_text="%", row=1, col=2,
                     gridcolor=COLORS['grid'])
    fig.update_xaxes(range=[0, 105], title_text="%", row=1, col=3,
                     gridcolor=COLORS['grid'])
    fig.update_xaxes(range=[0, 105], title_text="", row=1, col=4,
                     gridcolor=COLORS['grid'], showticklabels=False)

    # Style subplot titles
    for annotation in fig['layout']['annotations']:
        annotation['font'] = dict(size=13, color=COLORS['text'])

    return fig, df


def create_quality_pcr_scatter(data):
    """Scatter plot: Quality metrics vs PCR success rate per student."""

    # Calculate average quality metrics per student
    quality_summary = data.nanodrop.groupby('Student_Code').agg({
        'Ratio_260_280': 'mean',
        'Ratio_260_230': 'mean'
    }).reset_index()

    # Calculate PCR success rate per student
    pcr_by_student = data.gel_results[
        data.gel_results['Sample_Type'] != 'control+'
    ].groupby('Student_Code').agg({
        'PCR_Success': ['sum', 'count']
    }).reset_index()
    pcr_by_student.columns = ['Student_Code', 'PCR_Success', 'PCR_Total']
    pcr_by_student['PCR_Rate'] = pcr_by_student['PCR_Success'] / pcr_by_student['PCR_Total'] * 100

    # Merge
    df = quality_summary.merge(pcr_by_student, on='Student_Code', how='left')
    df['PCR_Rate'] = df['PCR_Rate'].fillna(0)

    # Create scatter plot
    fig = go.Figure()

    # Color by PCR success rate
    colors = [
        COLORS['success'] if r >= 66 else COLORS['partial'] if r >= 33 else COLORS['fail']
        for r in df['PCR_Rate']
    ]

    fig.add_trace(go.Scatter(
        x=df['Ratio_260_280'],
        y=df['Ratio_260_230'],
        mode='markers+text',
        marker=dict(
            size=df['PCR_Rate'] / 5 + 10,  # Size by PCR rate
            color=colors,
            line=dict(width=1, color='black')
        ),
        text=df['Student_Code'],
        textposition='top center',
        hovertemplate=(
            '<b>%{text}</b><br>' +
            '260/280: %{x:.2f}<br>' +
            '260/230: %{y:.2f}<br>' +
            'PCR Success: %{customdata:.0f}%<extra></extra>'
        ),
        customdata=df['PCR_Rate']
    ))

    # Add ideal zone
    fig.add_shape(
        type="rect",
        x0=1.8, x1=2.0, y0=2.0, y1=4.0,
        fillcolor='rgba(80, 250, 123, 0.15)',
        line=dict(width=0),
        layer='below'
    )

    # Reference lines
    fig.add_vline(x=1.8, line_dash="dash", line_color='rgba(80, 250, 123, 0.5)', line_width=2)
    fig.add_vline(x=2.0, line_dash="dash", line_color='rgba(80, 250, 123, 0.5)', line_width=2)
    fig.add_hline(y=2.0, line_dash="dash", line_color='rgba(80, 250, 123, 0.5)', line_width=2)

    fig.add_annotation(
        x=1.9, y=3.0,
        text="Ideal Zone",
        showarrow=False,
        font=dict(size=10, color='#50fa7b'),
        bgcolor='rgba(255,255,255,0.8)'
    )

    fig.update_layout(
        title=dict(
            text="DNA Quality vs PCR Success by Student",
            subtitle=dict(text="Bubble size = PCR success rate | Color: green >66%, orange 33-66%, red <33%"),
            font=dict(size=16, color=COLORS['text']),
            x=0.5,
            xanchor='center'
        ),
        xaxis_title="Average 260/280 Ratio",
        yaxis_title="Average 260/230 Ratio",
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(color=COLORS['text'], size=12),
        height=450,
        margin=dict(l=60, r=40, t=100, b=60),
        xaxis=dict(gridcolor=COLORS['grid'], range=[0, 3]),
        yaxis=dict(gridcolor=COLORS['grid'], range=[-2, 4])
    )

    return fig, df


def main():
    """Main function to create pipeline success visualizations."""
    print("Pipeline Success Correlation")
    print("=" * 50)

    # Load data
    data = LabData().load_all()

    # Create plots
    print("\nCreating pipeline funnel chart...")
    funnel_fig, funnel_df = create_pipeline_funnel(data)

    print("Creating quality vs PCR scatter...")
    scatter_fig, scatter_df = create_quality_pcr_scatter(data)

    # Save plots
    output_dir = Path(__file__).parent.parent.parent / "results" / "lab_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    funnel_file = output_dir / "pipeline_success.json"
    with open(funnel_file, 'w') as f:
        f.write(funnel_fig.to_json())
    print(f"Saved: {funnel_file}")

    scatter_file = output_dir / "quality_vs_pcr.json"
    with open(scatter_file, 'w') as f:
        f.write(scatter_fig.to_json())
    print(f"Saved: {scatter_file}")

    # Print summary
    print("\n" + "=" * 50)
    print("Pipeline Success Summary:")
    print("=" * 50)

    # Students who got species ID
    success_students = funnel_df[funnel_df['Has_Species']]['Student'].tolist()
    print(f"\nStudents with successful species ID: {', '.join(success_students) if success_students else 'None'}")

    # Correlation insight
    print("\nKey Insight:")
    high_quality = funnel_df[funnel_df['Quality_Score'] >= 2.5]
    low_quality = funnel_df[funnel_df['Quality_Score'] < 2.5]

    if len(high_quality) > 0 and len(low_quality) > 0:
        high_pcr = high_quality['PCR_Rate'].mean()
        low_pcr = low_quality['PCR_Rate'].mean()
        print(f"  Higher quality samples (score >= 2.5): {high_pcr:.0f}% avg PCR success")
        print(f"  Lower quality samples (score < 2.5): {low_pcr:.0f}% avg PCR success")

    print("\nDone!")


if __name__ == "__main__":
    main()
