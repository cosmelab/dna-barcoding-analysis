#!/usr/bin/env python3
"""
Generate personalized student report for 5-minute presentation.
Usage: python generate_student_report.py <STUDENT_CODE>
"""

import sys
from pathlib import Path
from datetime import datetime
from data_loader import LabData
import numpy as np


# Colors matching the main report
COLORS = {
    'purple': '#bd93f9',
    'cyan': '#8be9fd',
    'green': '#50fa7b',
    'orange': '#ffb86c',
    'red': '#ff5555',
    'gold': '#f1c40f',
    'text': '#2d3436',
}


def get_student_data(data, student_code):
    """Extract all data for a specific student."""
    student = {}

    # Basic info from teams
    if data.teams is not None:
        team_row = data.teams[data.teams['Student_Code'] == student_code]
        if len(team_row) > 0:
            student['team'] = team_row['Team'].values[0]
        else:
            student['team'] = 'Unknown'

    # Column extraction data
    if data.column_extraction is not None:
        col_data = data.column_extraction[data.column_extraction['Student_Code'] == student_code]
        if len(col_data) > 0:
            student['column_extraction'] = col_data.to_dict('records')
            student['avg_dna_yield'] = col_data['Total_DNA_ng'].mean()
        else:
            student['column_extraction'] = []
            student['avg_dna_yield'] = 0

    # NanoDrop quality data
    if data.nanodrop is not None:
        nano_data = data.nanodrop[data.nanodrop['Student_Code'] == student_code]
        student['nanodrop'] = nano_data.to_dict('records')

    # Quality assessment
    if data.quality_assessment is not None:
        qa_data = data.quality_assessment[data.quality_assessment['Student_Code'] == student_code]
        student['quality'] = qa_data.to_dict('records')
        # Count quality grades
        if len(qa_data) > 0:
            student['quality_counts'] = qa_data['Quality_Grade'].value_counts().to_dict()
        else:
            student['quality_counts'] = {}

    # PCR/Gel results
    if data.gel_results is not None:
        gel_data = data.gel_results[data.gel_results['Student_Code'] == student_code]
        # Exclude controls
        gel_data = gel_data[~gel_data['Sample_Type'].str.contains('control', case=False, na=False)]
        student['pcr_results'] = gel_data.to_dict('records')
        student['pcr_success'] = gel_data['PCR_Success'].sum()
        student['pcr_total'] = len(gel_data)
        student['pcr_rate'] = (student['pcr_success'] / student['pcr_total'] * 100) if student['pcr_total'] > 0 else 0

    # Sequencing results
    if data.sequencing is not None:
        seq_data = data.sequencing[data.sequencing['Student_Code'] == student_code]
        student['sequencing'] = seq_data.to_dict('records')
        if len(seq_data) > 0:
            student['seq_passed'] = int(seq_data['Passed_QC'].values[0]) if 'Passed_QC' in seq_data.columns else 0
            student['seq_total'] = int(seq_data['Total_Sequences'].values[0]) if 'Total_Sequences' in seq_data.columns else 0
            student['species_id'] = seq_data['Species_Identified'].values[0] if 'Species_Identified' in seq_data.columns else None
        else:
            student['seq_passed'] = 0
            student['seq_total'] = 0
            student['species_id'] = None

    return student


def get_class_averages(data):
    """Calculate class averages for comparison."""
    averages = {}

    # DNA yield average
    if data.column_extraction is not None:
        averages['dna_yield'] = data.column_extraction['Total_DNA_ng'].mean()

    # PCR success rate
    if data.gel_results is not None:
        gel = data.gel_results[~data.gel_results['Sample_Type'].str.contains('control', case=False, na=False)]
        averages['pcr_rate'] = (gel['PCR_Success'].sum() / len(gel) * 100) if len(gel) > 0 else 0

    return averages


def get_team_standing(data, team):
    """Get team performance metrics."""
    if data.column_extraction is None:
        return {}

    team_data = data.column_extraction[data.column_extraction['Team'] == team]
    other_team = 'Magnet' if team == 'Spin' else 'Spin'
    other_data = data.column_extraction[data.column_extraction['Team'] == other_team]

    return {
        'team': team,
        'team_yield': team_data['Total_DNA_ng'].mean() if len(team_data) > 0 else 0,
        'other_team': other_team,
        'other_yield': other_data['Total_DNA_ng'].mean() if len(other_data) > 0 else 0,
    }


def generate_html_report(student_code, student, class_avg, team_standing):
    """Generate the HTML report."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Determine performance indicators
    yield_vs_class = "above" if student.get('avg_dna_yield', 0) > class_avg.get('dna_yield', 0) else "below"
    pcr_vs_class = "above" if student.get('pcr_rate', 0) > class_avg.get('pcr_rate', 0) else "below"

    # Team color
    team_color = COLORS['purple'] if student.get('team') == 'Spin' else COLORS['gold']

    # PCR results table rows
    pcr_rows = ""
    for pcr in student.get('pcr_results', []):
        status = "✓" if pcr.get('PCR_Success') else "✗"
        color = COLORS['green'] if pcr.get('PCR_Success') else COLORS['red']
        pcr_rows += f"""
        <tr>
            <td>{pcr.get('Sample_Type', 'N/A')}</td>
            <td style="color: {color}; font-weight: bold;">{status}</td>
        </tr>
        """

    # Quality breakdown
    quality_html = ""
    for grade, count in student.get('quality_counts', {}).items():
        color = COLORS['green'] if grade in ['Good', 'Excellent'] else (COLORS['orange'] if grade == 'Fair' else COLORS['red'])
        quality_html += f'<span style="color: {color};">{grade}: {count}</span> &nbsp; '

    # Extraction results table
    extraction_rows = ""
    for ext in student.get('column_extraction', []):
        extraction_rows += f"""
        <tr>
            <td>{ext.get('Species', 'N/A')}</td>
            <td>{ext.get('Concentration_ng_uL', 0):.1f}</td>
            <td>{ext.get('Total_DNA_ng', 0):.1f}</td>
        </tr>
        """

    # Sequencing results
    seq_html = ""
    seq_passed = student.get('seq_passed', 0)
    seq_total = student.get('seq_total', 0)
    species_id = student.get('species_id')

    if seq_total > 0:
        pass_rate = (seq_passed / seq_total * 100) if seq_total > 0 else 0
        status_color = COLORS['green'] if pass_rate >= 50 else COLORS['red']
        seq_html = f"""
        <div class="metric-grid">
            <div class="metric">
                <div class="value">{seq_passed}/{seq_total}</div>
                <div class="label">Sequences Passed QC</div>
            </div>
            <div class="metric">
                <div class="value" style="color: {status_color};">{pass_rate:.0f}%</div>
                <div class="label">QC Pass Rate</div>
            </div>
        </div>
        """
        if species_id and str(species_id) != 'nan':
            seq_html += f"""
            <div class="highlight">
                <strong>Species Identified:</strong> <em>{species_id}</em>
            </div>
            """
        else:
            seq_html += """
            <div class="highlight" style="border-left-color: #ffb86c;">
                <strong>Note:</strong> No consensus sequences were generated for species identification.
            </div>
            """
    else:
        seq_html = "<p>No sequencing data available for this student.</p>"

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Lab Report - {student_code}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
            line-height: 1.6;
            color: {COLORS['text']};
            background: linear-gradient(135deg, #f5f0ff 0%, #f0fcff 50%, #fff0f7 100%);
            background-attachment: fixed;
            padding: 2rem;
        }}
        .container {{ max-width: 900px; margin: 0 auto; }}
        header {{
            background: #2d2d44;
            color: white;
            padding: 2rem;
            border-radius: 12px;
            text-align: center;
            margin-bottom: 2rem;
        }}
        header h1 {{ font-size: 2.5rem; margin-bottom: 0.5rem; }}
        .team-badge {{
            display: inline-block;
            background: {team_color};
            color: #2d2d44;
            padding: 0.3rem 1rem;
            border-radius: 20px;
            font-weight: bold;
            margin-top: 0.5rem;
        }}
        .card {{
            background: white;
            border-radius: 12px;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }}
        .card h2 {{
            color: #2d2d44;
            border-bottom: 3px solid {COLORS['purple']};
            padding-bottom: 0.5rem;
            margin-bottom: 1rem;
        }}
        .metric-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1rem;
            margin: 1rem 0;
        }}
        .metric {{
            background: #f8f9fa;
            padding: 1rem;
            border-radius: 8px;
            text-align: center;
        }}
        .metric .value {{
            font-size: 2rem;
            font-weight: bold;
            color: {COLORS['purple']};
        }}
        .metric .label {{ color: #636e72; font-size: 0.9rem; }}
        .comparison {{
            font-size: 0.85rem;
            margin-top: 0.3rem;
        }}
        .above {{ color: {COLORS['green']}; }}
        .below {{ color: {COLORS['orange']}; }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 1rem 0;
        }}
        th, td {{
            padding: 0.75rem;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        th {{ background: #f8f9fa; font-weight: 600; }}
        .highlight {{
            background: linear-gradient(90deg, {COLORS['purple']}22, {COLORS['cyan']}22);
            border-left: 4px solid {COLORS['purple']};
            padding: 1rem;
            border-radius: 0 8px 8px 0;
            margin: 1rem 0;
        }}
        .team-comparison {{
            display: flex;
            justify-content: space-around;
            margin: 1rem 0;
        }}
        .team-box {{
            text-align: center;
            padding: 1rem;
            border-radius: 8px;
            min-width: 150px;
        }}
        footer {{
            text-align: center;
            color: #636e72;
            font-size: 0.9rem;
            margin-top: 2rem;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>Lab Report: {student_code}</h1>
            <p>ENTM 201L DNA Barcoding Lab</p>
            <div class="team-badge">Team {student.get('team', 'Unknown')}</div>
        </header>

        <!-- Summary Metrics -->
        <div class="card">
            <h2>Your Results at a Glance</h2>
            <div class="metric-grid">
                <div class="metric">
                    <div class="value">{student.get('avg_dna_yield', 0):.0f}</div>
                    <div class="label">Avg DNA Yield (ng)</div>
                    <div class="comparison {yield_vs_class}">
                        {yield_vs_class.capitalize()} class avg ({class_avg.get('dna_yield', 0):.0f} ng)
                    </div>
                </div>
                <div class="metric">
                    <div class="value">{student.get('pcr_rate', 0):.0f}%</div>
                    <div class="label">PCR Success Rate</div>
                    <div class="comparison {pcr_vs_class}">
                        {pcr_vs_class.capitalize()} class avg ({class_avg.get('pcr_rate', 0):.0f}%)
                    </div>
                </div>
                <div class="metric">
                    <div class="value">{student.get('pcr_success', 0)}/{student.get('pcr_total', 0)}</div>
                    <div class="label">Successful PCRs</div>
                </div>
                <div class="metric">
                    <div class="value">{student.get('seq_passed', 0)}/{student.get('seq_total', 0)}</div>
                    <div class="label">Sequences Passed QC</div>
                </div>
            </div>
        </div>

        <!-- DNA Extraction -->
        <div class="card">
            <h2>1. DNA Extraction Results</h2>
            <p>Column-based extraction using Zymo Quick-DNA Miniprep Plus Kit.</p>

            <table>
                <thead>
                    <tr>
                        <th>Species</th>
                        <th>Concentration (ng/µL)</th>
                        <th>Total DNA (ng)</th>
                    </tr>
                </thead>
                <tbody>
                    {extraction_rows if extraction_rows else "<tr><td colspan='3'>No extraction data available</td></tr>"}
                </tbody>
            </table>

            <div class="highlight">
                <strong>DNA Quality Assessment:</strong> {quality_html if quality_html else "No quality data available"}
            </div>
        </div>

        <!-- PCR Results -->
        <div class="card">
            <h2>2. PCR Amplification</h2>
            <p>COI gene amplification (~710 bp) for species identification.</p>

            <table>
                <thead>
                    <tr>
                        <th>Preservation Method</th>
                        <th>Result</th>
                    </tr>
                </thead>
                <tbody>
                    {pcr_rows if pcr_rows else "<tr><td colspan='2'>No PCR data available</td></tr>"}
                </tbody>
            </table>

            <div class="highlight">
                <strong>Your PCR Success Rate:</strong> {student.get('pcr_rate', 0):.0f}%
                ({student.get('pcr_success', 0)} of {student.get('pcr_total', 0)} reactions)
            </div>
        </div>

        <!-- Sequencing Results -->
        <div class="card">
            <h2>3. Sequencing & Species ID</h2>
            <p>Sanger sequencing of successful PCR products.</p>
            {seq_html}
        </div>

        <!-- Team Competition -->
        <div class="card">
            <h2>4. Team Challenge: Spin vs Magnet</h2>
            <div class="team-comparison">
                <div class="team-box" style="background: {COLORS['purple']}33;">
                    <div style="font-size: 1.5rem; font-weight: bold;">Team Spin</div>
                    <div style="font-size: 2rem; color: {COLORS['purple']};">{team_standing.get('team_yield' if team_standing.get('team') == 'Spin' else 'other_yield', 0):.0f} ng</div>
                    <div>Avg DNA Yield</div>
                </div>
                <div style="display: flex; align-items: center; font-size: 2rem; color: #636e72;">vs</div>
                <div class="team-box" style="background: {COLORS['gold']}33;">
                    <div style="font-size: 1.5rem; font-weight: bold;">Team Magnet</div>
                    <div style="font-size: 2rem; color: {COLORS['gold']};">{team_standing.get('other_yield' if team_standing.get('team') == 'Spin' else 'team_yield', 0):.0f} ng</div>
                    <div>Avg DNA Yield</div>
                </div>
            </div>
            <div class="highlight">
                <strong>Your contribution:</strong> You are on Team {student.get('team', 'Unknown')}
                with an average yield of {student.get('avg_dna_yield', 0):.0f} ng.
            </div>
        </div>

        <!-- Presentation Tips -->
        <div class="card">
            <h2>5. Presentation Talking Points</h2>
            <ul style="margin-left: 1.5rem;">
                <li><strong>Introduction:</strong> Briefly explain DNA barcoding and the COI gene</li>
                <li><strong>Your Results:</strong> Highlight your extraction yields and PCR success</li>
                <li><strong>Challenges:</strong> Discuss any failed samples and possible reasons</li>
                <li><strong>Species ID:</strong> What species did you identify? Any surprises?</li>
                <li><strong>Team Performance:</strong> How did your team compare?</li>
            </ul>
        </div>

        <footer>
            <p>Generated on {timestamp} | ENTM 201L DNA Barcoding Lab</p>
        </footer>
    </div>
</body>
</html>
'''
    return html


def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_student_report.py <STUDENT_CODE>")
        print("Available codes: BR, HV, JA, JM, JR, KG, MA, TW, WA, WL")
        sys.exit(1)

    student_code = sys.argv[1].upper()
    valid_codes = ['BR', 'HV', 'JA', 'JM', 'JR', 'KG', 'MA', 'TW', 'WA', 'WL']

    if student_code not in valid_codes:
        print(f"Error: '{student_code}' is not a valid student code.")
        print(f"Valid codes: {', '.join(valid_codes)}")
        sys.exit(1)

    print(f"Generating report for student: {student_code}")

    # Load data
    data = LabData().load_all()

    # Get student-specific data
    student = get_student_data(data, student_code)
    class_avg = get_class_averages(data)
    team_standing = get_team_standing(data, student.get('team', 'Spin'))

    # Generate HTML
    html = generate_html_report(student_code, student, class_avg, team_standing)

    # Save report
    output_dir = Path(__file__).parent.parent.parent / "results" / "student_reports"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{student_code}_report.html"
    with open(output_file, 'w') as f:
        f.write(html)

    print(f"Report saved: {output_file}")


if __name__ == "__main__":
    main()
