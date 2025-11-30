#!/usr/bin/env python3
"""
PCR Results and Gel Image Analysis

Creates visualizations of PCR success rates and embeds gel electrophoresis images.
Shows which DNA extraction methods yielded successful PCR amplification.

Educational context:
- Gel electrophoresis separates DNA by size using electric current
- A visible band = successful PCR amplification (DNA present)
- No band = failed amplification (no DNA or amplification issues)
- Band intensity indicates relative DNA concentration
"""

import pandas as pd
import plotly.graph_objects as go
import base64
from pathlib import Path
from data_loader import LabData


# Clean white style
COLORS = {
    'success': '#50fa7b',     # Green - successful PCR
    'failure': '#ff5555',     # Red - failed PCR
    'background': 'white',
    'text': '#2d3436',
    'grid': '#e0e0e0'
}

# Preservation method colors (consistent with extraction plots)
PRESERVATION_COLORS = {
    'Frozen': '#8be9fd',      # Cyan
    'Ethanol': '#bd93f9',     # Purple
    'Silica': '#ffb86c'       # Orange
}


def create_pcr_success_chart(gel_data):
    """
    Create heatmap grid showing PCR success by student and preservation method.

    Args:
        gel_data: DataFrame with PCR results (Student, PCR_Success columns)

    Returns:
        Plotly figure object
    """
    import numpy as np

    # Exclude controls
    student_data = gel_data[~gel_data['Sample_Type'].str.contains('control', case=False, na=False)]

    # Create pivot table: students x methods
    pivot = student_data.pivot_table(
        index='Student_Code',
        columns='Sample_Type',
        values='PCR_Success',
        aggfunc='first'
    )

    # Order methods consistently
    methods = ['Frozen', 'Ethanol', 'Silica']
    pivot = pivot.reindex(columns=[m for m in methods if m in pivot.columns])

    # Sort students alphabetically
    pivot = pivot.sort_index()

    # Create heatmap values (1=success, 0=fail)
    z_values = pivot.values

    # Create text annotations (✓ or ✗)
    text_values = np.where(z_values == 1, '✓', '✗')

    # Custom colorscale: red for fail, green for success
    colorscale = [[0, COLORS['failure']], [1, COLORS['success']]]

    fig = go.Figure(data=go.Heatmap(
        z=z_values,
        x=pivot.columns.tolist(),
        y=pivot.index.tolist(),
        colorscale=colorscale,
        showscale=False,
        text=text_values,
        texttemplate='%{text}',
        textfont=dict(size=16, color='white'),
        hovertemplate='%{y}: %{x}<br>%{text}<extra></extra>',
        xgap=3,
        ygap=3
    ))

    # Update layout
    fig.update_layout(
        title=dict(
            text='PCR Success by Student & Method',
            x=0.5,
            xanchor='center',
            font=dict(size=18, color=COLORS['text'])
        ),
        xaxis=dict(
            title='Preservation Method',
            tickfont=dict(size=12, color=COLORS['text']),
            side='bottom'
        ),
        yaxis=dict(
            title='Student',
            tickfont=dict(size=12, color=COLORS['text']),
            autorange='reversed'
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(color=COLORS['text'], size=12),
        height=400,
        margin=dict(l=60, r=40, t=80, b=60)
    )

    return fig


def create_success_by_method_chart(gel_data):
    """
    Create bar chart showing PCR success rate by extraction method.

    Args:
        gel_data: DataFrame with PCR results

    Returns:
        Plotly figure object
    """
    # Exclude controls
    sample_data = gel_data[~gel_data['Sample_Type'].str.contains('control', case=False, na=False)]

    # Calculate success rate per method
    method_stats = sample_data.groupby('Sample_Type').agg({
        'PCR_Success': ['sum', 'count']
    }).reset_index()
    method_stats.columns = ['Method', 'Successful', 'Total']
    method_stats['Success_Rate'] = (method_stats['Successful'] / method_stats['Total'] * 100)

    # Sort by success rate
    method_stats = method_stats.sort_values('Success_Rate', ascending=False)

    # Create bar chart with colors matching preservation methods
    fig = go.Figure()

    # Get colors for each method
    bar_colors = [PRESERVATION_COLORS.get(method, COLORS['success']) for method in method_stats['Method']]

    fig.add_trace(go.Bar(
        x=method_stats['Method'],
        y=method_stats['Success_Rate'],
        marker_color=bar_colors,
        text=[f"{rate:.1f}%" for rate in method_stats['Success_Rate']],
        textposition='outside',
        hovertemplate='%{y:.1f}%<extra></extra>'
    ))

    # Update layout - clean white style
    fig.update_layout(
        title=dict(
            text='PCR Success Rate by Extraction Method',
            x=0.5,
            xanchor='center',
            font=dict(size=18, color=COLORS['text'])
        ),
        xaxis=dict(
            title='Extraction Method',
            gridcolor=COLORS['grid'],
            tickfont=dict(size=11, color=COLORS['text']),
        ),
        yaxis=dict(
            title='Success Rate (%)',
            gridcolor=COLORS['grid'],
            tickfont=dict(size=11, color=COLORS['text']),
            range=[0, 105]
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(color=COLORS['text'], size=12),
        height=400,
        margin=dict(l=60, r=40, t=80, b=60)
    )

    return fig


def encode_image_to_base64(image_path):
    """
    Encode image file to base64 string for HTML embedding.

    Args:
        image_path: Path to image file

    Returns:
        Base64 encoded string (data URI format)
    """
    image_path = Path(image_path)

    if not image_path.exists():
        return None

    with open(image_path, 'rb') as f:
        image_data = f.read()

    # Encode to base64
    base64_data = base64.b64encode(image_data).decode('utf-8')

    # Return as data URI
    return f"data:image/png;base64,{base64_data}"


def get_gel_images_with_metadata(data_loader):
    """
    Get all gel images with their metadata and base64 encoding.

    Args:
        data_loader: LabData instance

    Returns:
        List of dictionaries with gel image info
    """
    gel_images = []
    gel_paths = data_loader.get_gel_images()

    # Define metadata for each gel
    gel_metadata = {
        'coi_gel1': {
            'title': 'COI PCR Gel 1: WL, KG, MA, HV',
            'students': ['WL', 'KG', 'MA', 'HV'],
            'description': 'First gel showing PCR results for four students. Each student has three samples (Frozen, Ethanol, Silica) plus positive control.'
        },
        'coi_gel2': {
            'title': 'COI PCR Gel 2: TW, BR, WA, JR',
            'students': ['TW', 'BR', 'WA', 'JR'],
            'description': 'Second gel showing PCR results for four students testing different preservation methods.'
        },
        'coi_gel3': {
            'title': 'COI PCR Gel 3: JA, MA',
            'students': ['JA', 'MA'],
            'description': 'Third gel showing PCR results for two students, including repeat samples from MA.'
        },
        'wolbachia': {
            'title': 'Wolbachia PCR Gel',
            'students': ['Various'],
            'description': 'Wolbachia screening gel testing for bacterial endosymbiont infection in samples.'
        }
    }

    for gel_key, gel_path in gel_paths.items():
        metadata = gel_metadata.get(gel_key, {})

        gel_info = {
            'key': gel_key,
            'title': metadata.get('title', gel_path.name),
            'students': metadata.get('students', []),
            'description': metadata.get('description', ''),
            'filename': gel_path.name,
            'path': str(gel_path),
            'base64': encode_image_to_base64(gel_path)
        }

        gel_images.append(gel_info)

    return gel_images


def create_educational_text():
    """
    Generate educational text explaining gel electrophoresis and PCR results.

    Returns:
        Dictionary with educational content sections
    """
    return {
        'gel_basics': {
            'title': 'Understanding Gel Electrophoresis',
            'content': """
Gel electrophoresis separates DNA fragments by size using an electric field:
- DNA is negatively charged and moves toward the positive electrode
- Smaller fragments move faster through the gel matrix
- DNA is visualized using fluorescent dye (GelRed or SYBR)
- UV light makes the dye glow, showing bands where DNA is present
            """
        },
        'band_interpretation': {
            'title': 'Reading the Gel',
            'content': """
What the bands tell us:
- <strong style="color: #50fa7b;">Visible band = Successful PCR</strong>
  - DNA was successfully amplified
  - Band appears at expected size (~710 bp for COI)
  - Brighter band = more DNA

- <strong style="color: #ff5555;">No band = Failed PCR</strong>
  - No amplification occurred
  - Possible causes: poor DNA quality, PCR inhibitors, low DNA concentration
  - Sample may need re-extraction or optimization
            """
        },
        'controls': {
            'title': 'Positive Controls',
            'content': """
Control lanes verify the PCR worked correctly:
- Positive control (control+): Known DNA template that should always amplify
- If control+ shows a band, the PCR reagents and conditions are working
- If control+ fails, something is wrong with the PCR setup itself
- Student samples are only valid when control+ is successful
            """
        },
        'success_factors': {
            'title': 'Factors Affecting PCR Success',
            'content': """
Why some samples succeed while others fail:
1. <strong>DNA Quality</strong>: Degraded DNA won't amplify well
2. <strong>DNA Concentration</strong>: Too little DNA = no amplification
3. <strong>Inhibitors</strong>: Ethanol, salts, or other contaminants can block PCR
4. <strong>Extraction Method</strong>: Different methods work better for different samples
5. <strong>Sample Preservation</strong>: Frozen, ethanol, or silica each affect DNA quality
            """
        }
    }


def main():
    """
    Main function to generate all PCR analysis outputs.
    """
    print("=" * 60)
    print("PCR Results and Gel Image Analysis")
    print("=" * 60)

    # Load data
    print("\nLoading data...")
    data = LabData()
    data.load_all()

    if data.gel_results is None:
        print("Error: Could not load gel results data")
        return

    # Create output directory
    output_dir = Path(__file__).parent.parent.parent / "results" / "lab_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create visualizations (JSON only - HTML will be combined in report)
    print("\nCreating visualizations...")

    # 1. PCR success by student
    print("  - PCR success rate by student...")
    fig_student = create_pcr_success_chart(data.gel_results)
    fig_student.write_json(output_dir / "pcr_success_by_student.json")

    # 2. PCR success by extraction method
    print("  - PCR success rate by method...")
    fig_method = create_success_by_method_chart(data.gel_results)
    fig_method.write_json(output_dir / "pcr_success_by_method.json")

    # 3. Process gel images
    print("\nProcessing gel images...")
    gel_images = get_gel_images_with_metadata(data)

    for gel in gel_images:
        print(f"  - {gel['filename']}: {len(gel['base64']) if gel['base64'] else 0} bytes encoded")

    # Save gel image data as JSON for embedding
    import json
    with open(output_dir / "gel_images.json", 'w') as f:
        json.dump(gel_images, f, indent=2)

    # 4. Generate educational content
    print("\nGenerating educational content...")
    edu_text = create_educational_text()
    with open(output_dir / "educational_text.json", 'w') as f:
        json.dump(edu_text, f, indent=2)

    # Summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)

    # Overall success rate
    sample_data = data.gel_results[~data.gel_results['Sample_Type'].str.contains('control', case=False, na=False)]
    total_pcr = len(sample_data)
    successful_pcr = sample_data['PCR_Success'].sum()
    success_rate = (successful_pcr / total_pcr * 100) if total_pcr > 0 else 0

    print(f"\nOverall PCR Results:")
    print(f"  Total reactions: {total_pcr}")
    print(f"  Successful: {successful_pcr}")
    print(f"  Failed: {total_pcr - successful_pcr}")
    print(f"  Success rate: {success_rate:.1f}%")

    # By extraction method
    print(f"\nBy Extraction Method:")
    method_stats = sample_data.groupby('Sample_Type')['PCR_Success'].agg(['sum', 'count'])
    for method, row in method_stats.iterrows():
        method_rate = (row['sum'] / row['count'] * 100) if row['count'] > 0 else 0
        print(f"  {method}: {row['sum']}/{row['count']} ({method_rate:.1f}%)")

    # By student
    print(f"\nBy Student:")
    student_stats = sample_data.groupby('Student_Code')['PCR_Success'].agg(['sum', 'count'])
    for student, row in student_stats.iterrows():
        student_rate = (row['sum'] / row['count'] * 100) if row['count'] > 0 else 0
        print(f"  {student}: {row['sum']}/{row['count']} ({student_rate:.1f}%)")

    print(f"\n" + "=" * 60)
    print(f"Output saved to: {output_dir}")
    print("  - pcr_success_by_student.json / .html")
    print("  - pcr_success_by_method.json / .html")
    print("  - gel_images.json (base64 encoded)")
    print("  - educational_text.json")
    print("=" * 60)


if __name__ == "__main__":
    main()
