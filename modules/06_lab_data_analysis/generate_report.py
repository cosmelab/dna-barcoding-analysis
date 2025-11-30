#!/usr/bin/env python3
"""
Generate Lab Analysis Report - ONE HTML file with ALL plots inline.
Linear story: Extraction → Quality → PCR → Sequencing → Teams
"""

import json
from pathlib import Path
from datetime import datetime
import plotly.io as pio


def generate_report(output_dir):
    """Generate single HTML report with all plots embedded."""

    output_dir = Path(output_dir)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Load all plots and convert to HTML divs
    plots_html = {}
    plot_files = [
        'hmw_extraction_comparison',
        'column_extraction_comparison',
        'quality_plot',
        'pcr_success_by_student',
        'pcr_success_by_method',
        'sequencing_pass_rate',
        'sequencing_breakdown',
        'team_comparison',
        'pipeline_success',
        'quality_vs_pcr'
    ]

    print("Loading plots...")
    for name in plot_files:
        json_file = output_dir / f"{name}.json"
        if json_file.exists():
            with open(json_file) as f:
                fig_dict = json.load(f)
            # Get full HTML with embedded Plotly.js for first plot, without for rest
            if not plots_html:  # First plot - include plotly.js
                plots_html[name] = pio.to_html(fig_dict, include_plotlyjs='cdn', full_html=False)
            else:
                plots_html[name] = pio.to_html(fig_dict, include_plotlyjs=False, full_html=False)
            print(f"  Loaded: {name}")

    # Load gel images
    gel_html = ""
    gel_file = output_dir / "gel_images.json"
    if gel_file.exists():
        with open(gel_file) as f:
            gels = json.load(f)
        for gel in gels:
            if gel.get('base64'):
                gel_html += f'''
                <div class="gel-card">
                    <img src="{gel['base64']}" alt="{gel['title']}">
                    <div class="gel-info">
                        <h4>{gel['title']}</h4>
                        <p>{gel['description']}</p>
                    </div>
                </div>
                '''

    # Build the report
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Lab Data Analysis - DNA Barcoding</title>
    <style>
        :root {{
            --purple: #bd93f9;
            --green: #50fa7b;
            --orange: #ffb86c;
            --red: #ff5555;
            --cyan: #8be9fd;
            --bg-dark: #2d2d44;
        }}
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
            line-height: 1.6;
            color: #2d3436;
            background: linear-gradient(135deg, #f5f0ff 0%, #f0fcff 50%, #fff0f7 100%);
            background-attachment: fixed;
        }}
        header {{
            background: var(--bg-dark);
            color: white;
            padding: 3rem 2rem;
            text-align: center;
        }}
        header h1 {{ font-size: 2.5rem; margin-bottom: 0.5rem; }}
        header p {{ color: var(--purple); font-size: 1.1rem; }}
        main {{
            max-width: 1000px;
            margin: 0 auto;
            padding: 2rem;
        }}
        .section {{
            background: white;
            border-radius: 12px;
            padding: 2rem;
            margin-bottom: 2rem;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }}
        .section h2 {{
            color: var(--bg-dark);
            font-size: 1.8rem;
            border-bottom: 3px solid var(--purple);
            padding-bottom: 0.5rem;
            margin-bottom: 1rem;
        }}
        .section h3 {{
            color: #636e72;
            margin: 2rem 0 1rem;
        }}
        .info-box {{
            background: #f5f0ff;
            border-left: 4px solid var(--purple);
            padding: 1rem 1.5rem;
            margin: 1.5rem 0;
            border-radius: 0 8px 8px 0;
        }}
        .info-box.good {{
            background: #f0fff5;
            border-left-color: var(--green);
        }}
        .info-box.warning {{
            background: #fff8f0;
            border-left-color: var(--orange);
        }}
        .info-box.bad {{
            background: #fff0f0;
            border-left-color: var(--red);
        }}
        .plot-container {{
            margin: 1.5rem 0;
            border-radius: 8px;
            overflow: hidden;
        }}
        .gel-card {{
            background: #f8f9fa;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            margin: 1.5rem 0;
        }}
        .gel-card img {{ width: 100%; max-width: 800px; display: block; margin: 0 auto; }}
        .gel-info {{ padding: 1rem; text-align: center; }}
        .gel-info h4 {{ color: var(--bg-dark); margin-bottom: 0.5rem; }}
        .gel-info p {{ color: #636e72; font-size: 0.9rem; }}
        ul {{ margin: 1rem 0 1rem 2rem; }}
        li {{ margin-bottom: 0.5rem; }}
        .step-number {{
            display: inline-block;
            width: 32px;
            height: 32px;
            background: var(--purple);
            color: white;
            border-radius: 50%;
            text-align: center;
            line-height: 32px;
            font-weight: bold;
            margin-right: 0.5rem;
        }}
    </style>
</head>
<body>
    <header>
        <h1>Lab Data Analysis Report</h1>
        <p>ENTM 201L DNA Barcoding | Generated {timestamp}</p>
    </header>

    <main>
        <!-- STEP 1: DNA EXTRACTION -->
        <div class="section">
            <h2><span class="step-number">1</span> DNA Extraction</h2>

            <div class="info-box">
                <strong>What is DNA Extraction?</strong><br>
                DNA extraction isolates DNA from cells for downstream analysis. We used two methods:
                <ul>
                    <li><strong>HMW (Magnetic Beads):</strong> High molecular weight extraction using magnetic beads to bind DNA</li>
                    <li><strong>Column (Zymo Kit):</strong> Silica membrane columns for rapid purification</li>
                </ul>
            </div>

            <h3>HMW Extraction: Does Preservation Method Matter?</h3>
            <p>We tested three preservation methods: frozen (-80°C), 95% ethanol, and silica gel desiccant.</p>
            <div class="plot-container">
                {plots_html.get('hmw_extraction_comparison', '<p>Plot not available</p>')}
            </div>
            <div class="info-box">
                <strong>Result:</strong> No significant difference between preservation methods (ANOVA p = 0.86).
                All three methods work equally well for HMW extraction.
            </div>

            <h3>Column Extraction: Aedes vs Culex</h3>
            <p>We compared DNA yields between two mosquito genera using the Zymo Quick-DNA kit.</p>
            <div class="plot-container">
                {plots_html.get('column_extraction_comparison', '<p>Plot not available</p>')}
            </div>
            <div class="info-box good">
                <strong>Result:</strong> <em>Aedes albopictus</em> yields significantly more DNA than <em>Culex</em> spp.
                (t-test p = 0.002). This may be due to larger body size.
            </div>
        </div>

        <!-- STEP 2: DNA QUALITY -->
        <div class="section">
            <h2><span class="step-number">2</span> DNA Quality Assessment</h2>

            <div class="info-box">
                <strong>Understanding NanoDrop Ratios</strong><br>
                <ul>
                    <li><strong>260/280 ratio:</strong> Measures protein contamination. Ideal: 1.8-2.0</li>
                    <li><strong>260/230 ratio:</strong> Measures organic contamination (salts, EDTA, phenol). Ideal: >2.0</li>
                </ul>
                The green box on the plot shows the ideal quality zone where both ratios are optimal.
            </div>

            <h3>Quality by Extraction Method: Column vs HMW</h3>
            <p>We compared NanoDrop quality metrics between our two extraction methods.</p>
            <div class="plot-container">
                {plots_html.get('quality_plot', '<p>Plot not available</p>')}
            </div>

            <div class="info-box bad">
                <strong>Result:</strong> 78% of samples rated "Poor" quality across both methods.
                Most samples show severe organic contamination (260/230 < 2.0).
                <br><br>
                <strong>Column (Spin):</strong> 17 Poor, 3 Fair (0 in ideal range)<br>
                <strong>HMW (Magnetic Beads):</strong> 24 Poor, 5 Fair, 3 Good (0 in ideal range)
            </div>
        </div>

        <!-- STEP 3: PCR -->
        <div class="section">
            <h2><span class="step-number">3</span> PCR Amplification</h2>

            <div class="info-box">
                <strong>What is PCR?</strong><br>
                PCR (Polymerase Chain Reaction) amplifies specific DNA regions. We targeted the COI gene (~710 bp)
                for species identification. A visible band on the gel = successful amplification.
            </div>

            <h3>PCR Success by Student</h3>
            <div class="plot-container">
                {plots_html.get('pcr_success_by_student', '<p>Plot not available</p>')}
            </div>

            <h3>PCR Success by Extraction Method</h3>
            <div class="plot-container">
                {plots_html.get('pcr_success_by_method', '<p>Plot not available</p>')}
            </div>

            <div class="info-box warning">
                <strong>Result:</strong> 60% overall PCR success rate (18/30 reactions).
                HV and JR achieved 100% success. KG had 0% success.
                Ethanol preservation had the highest success rate (70%).
            </div>

            <h3>Gel Electrophoresis Images</h3>
            <p>Bright bands at ~710 bp indicate successful COI amplification. No band = failed PCR.</p>
            {gel_html}
        </div>

        <!-- STEP 4: SEQUENCING -->
        <div class="section">
            <h2><span class="step-number">4</span> DNA Sequencing</h2>

            <div class="info-box">
                <strong>Sequencing Pipeline</strong><br>
                Successful PCR products were sent for Sanger sequencing. Quality control (QC) filters out
                low-quality reads based on Phred scores. High-quality forward and reverse reads are combined
                into consensus sequences for species identification.
            </div>

            <h3>QC Pass Rate by Student</h3>
            <div class="plot-container">
                {plots_html.get('sequencing_pass_rate', '<p>Plot not available</p>')}
            </div>

            <h3>Sequencing Pipeline Breakdown</h3>
            <div class="plot-container">
                {plots_html.get('sequencing_breakdown', '<p>Plot not available</p>')}
            </div>

            <div class="info-box bad">
                <strong>Result:</strong> Only 40% of sequences passed QC (12/30).
                Just 4 consensus sequences were generated, identifying <em>Aedes albopictus</em> and <em>Culex</em> spp.
                Poor DNA quality from extraction likely contributed to low sequencing success.
            </div>
        </div>

        <!-- STEP 5: TEAM CHALLENGE -->
        <div class="section">
            <h2><span class="step-number">5</span> Team Challenge: Spin vs Magnet</h2>

            <div class="info-box">
                <strong>The Competition</strong><br>
                Students were divided into two teams to compare overall lab performance:
                <ul>
                    <li><strong>Team Spin:</strong> JR, HV, TW, JM, WL</li>
                    <li><strong>Team Magnet:</strong> MA, BR, WA, KG, JA</li>
                </ul>
            </div>

            <div class="plot-container">
                {plots_html.get('team_comparison', '<p>Plot not available</p>')}
            </div>
        </div>

        <!-- STEP 6: PIPELINE CORRELATION -->
        <div class="section">
            <h2><span class="step-number">6</span> Quality → Success Correlation</h2>

            <div class="info-box">
                <strong>Does DNA Quality Predict Success?</strong><br>
                We tracked each student's journey through the pipeline to see if upstream DNA quality
                predicts downstream success. The answer: <em>absolutely yes</em>.
            </div>

            <h3>Pipeline Success by Student</h3>
            <p>This chart shows each student's progress through all four stages: DNA quality, PCR success,
               sequencing QC, and species identification.</p>
            <div class="plot-container">
                {plots_html.get('pipeline_success', '<p>Plot not available</p>')}
            </div>

            <h3>DNA Quality vs PCR Success</h3>
            <p>Students with better NanoDrop ratios (closer to the ideal zone) had higher PCR success rates.
               Bubble size indicates PCR success rate.</p>
            <div class="plot-container">
                {plots_html.get('quality_vs_pcr', '<p>Plot not available</p>')}
            </div>

            <div class="info-box good">
                <strong>Key Finding:</strong> Only 3 students (HV, JM, WL) achieved successful species identification.
                All three had relatively better DNA quality or higher PCR success rates.
                The pipeline bottleneck is clearly at the extraction/quality step.
            </div>
        </div>

        <!-- SUMMARY -->
        <div class="section">
            <h2>Summary</h2>
            <ul>
                <li><strong>Extraction:</strong> Both HMW and column methods worked. Preservation method didn't significantly affect HMW yields.</li>
                <li><strong>Quality:</strong> Major issue - 78% poor quality due to organic contamination. This is the bottleneck.</li>
                <li><strong>PCR:</strong> 60% success rate. Some students (HV, JR) achieved 100%.</li>
                <li><strong>Sequencing:</strong> Low success (40% QC pass) due to upstream quality issues.</li>
            </ul>

            <div class="info-box">
                <strong>Key Lesson:</strong> DNA quality at extraction determines downstream success.
                Even with successful PCR bands, contaminated DNA often fails sequencing QC.
            </div>
        </div>

        <!-- RECOMMENDATIONS -->
        <div class="section">
            <h2>Recommendations: Improving DNA Quality</h2>

            <p>Since 78% of samples showed poor quality (mainly low 260/230 ratios indicating organic contamination),
               here are concrete steps to improve results in future experiments:</p>

            <div class="info-box good">
                <strong>Isopropanol/NaOAc Cleanup Protocol</strong><br>
                A simple ethanol precipitation can remove organic contaminants and dramatically improve 260/230 ratios.
                <br><br>
                <strong>Quick Overview:</strong>
                <ol>
                    <li>Add 0.1x volume 3M NaOAc (pH 5.2)</li>
                    <li>Add 0.7x volume isopropanol</li>
                    <li>Incubate 10-20 min at RT</li>
                    <li>Spin at low speed (3000-5000 x g) to pellet</li>
                    <li>Wash 2x with 70% ethanol</li>
                    <li>Air dry and resuspend in Tris buffer</li>
                </ol>
                <br>
                See the full protocol: <code>docs/hmw_cleanup_protocol.md</code>
            </div>

            <h3>Other Improvements</h3>
            <ul>
                <li><strong>Fresh reagents:</strong> Use fresh 70% ethanol (old ethanol absorbs water)</li>
                <li><strong>Extra washes:</strong> Add an extra wash step during column extraction</li>
                <li><strong>Elution optimization:</strong> Pre-warm elution buffer to 65°C</li>
                <li><strong>Sample handling:</strong> Minimize freeze-thaw cycles</li>
                <li><strong>Timing:</strong> Process samples promptly after collection</li>
            </ul>

            <div class="info-box warning">
                <strong>Trade-off:</strong> Cleanup protocols typically lose 10-20% of DNA.
                Only use if you have sufficient starting material or if contamination is severe.
            </div>
        </div>
    </main>
</body>
</html>
'''

    # Write report
    report_file = output_dir / "lab_report.html"
    with open(report_file, 'w') as f:
        f.write(html)

    print(f"\nReport saved: {report_file}")
    return report_file


def main():
    print("=" * 60)
    print("GENERATING LAB ANALYSIS REPORT")
    print("=" * 60)

    script_dir = Path(__file__).parent
    output_dir = script_dir.parent.parent / "results" / "lab_analysis"

    report_path = generate_report(output_dir)

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)
    print(f"\nOpen: {report_path}")


if __name__ == "__main__":
    main()
