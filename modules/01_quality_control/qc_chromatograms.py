#!/usr/bin/env python3
"""
Quality Control for Sanger Chromatograms
Analyzes .ab1 files and generates QC report with chromatogram visualization
"""

import os
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqIO import AbiIO
import pandas as pd
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import base64
from io import BytesIO

def plot_chromatogram(ab1_file, start_base=50, num_bases=150):
    """Generate chromatogram plot with sequence overlay

    Shows middle region by default (skipping poor quality start/end)
    """
    try:
        with open(ab1_file, 'rb') as f:
            record = AbiIO.AbiIterator(f).__next__()

        # Get trace data and base calls
        channels = {'DATA9': 'G', 'DATA10': 'A', 'DATA11': 'T', 'DATA12': 'C'}
        colors = {'G': 'black', 'A': 'green', 'T': 'red', 'C': 'blue'}

        # Get base calls and positions
        seq_record = SeqIO.read(ab1_file, 'abi')
        sequence = str(seq_record.seq)
        quality = seq_record.letter_annotations.get('phred_quality', [])

        # Get peak positions (where each base was called)
        peak_positions = record.annotations['abif_raw'].get('PLOC1', [])

        # Determine region to show
        end_base = min(start_base + num_bases, len(sequence))
        if start_base >= len(sequence):
            start_base = max(0, len(sequence) - num_bases)
            end_base = len(sequence)

        # Get trace region
        if peak_positions and len(peak_positions) > end_base:
            trace_start = peak_positions[start_base] if start_base < len(peak_positions) else 0
            trace_end = peak_positions[end_base-1] if end_base-1 < len(peak_positions) else len(record.annotations['abif_raw']['DATA9'])
        else:
            # Approximate if no peak positions
            trace_start = start_base * 10
            trace_end = end_base * 10

        fig, ax = plt.subplots(figsize=(18, 5))

        # Plot traces
        for channel, base in channels.items():
            if channel in record.annotations['abif_raw']:
                trace = record.annotations['abif_raw'][channel][trace_start:trace_end]
                ax.plot(trace, color=colors[base], label=base, linewidth=0.8, alpha=0.7)

        # Add base calls on top
        if peak_positions and len(peak_positions) > end_base:
            for i in range(start_base, end_base):
                if i < len(sequence) and i < len(peak_positions):
                    pos = peak_positions[i] - trace_start
                    base = sequence[i]

                    # Color code by quality
                    qual = quality[i] if i < len(quality) else 0
                    if qual >= 30:
                        qual_color = 'darkgreen'
                    elif qual >= 20:
                        qual_color = 'orange'
                    else:
                        qual_color = 'red'

                    ax.text(pos, ax.get_ylim()[1] * 0.95, base,
                           ha='center', va='top', fontsize=8,
                           color=qual_color, fontweight='bold')

        ax.set_xlabel('Trace Position', fontsize=10)
        ax.set_ylabel('Signal Intensity', fontsize=10)
        ax.set_title(f'{ab1_file.name} - Bases {start_base}-{end_base} (green=Q≥30, orange=Q≥20, red=Q<20)',
                    fontsize=11, fontweight='bold')
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.2)

        # Convert to base64
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=120, bbox_inches='tight')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        plt.close()

        return img_base64
    except Exception as e:
        print(f"  Warning: Could not generate chromatogram plot: {e}")
        return None

def analyze_chromatogram(ab1_file, output_dir=None):
    """Analyze a single .ab1 chromatogram file"""
    try:
        record = SeqIO.read(ab1_file, "abi")

        # Extract quality scores
        quality_scores = record.letter_annotations.get("phred_quality", [])

        # Calculate QC metrics
        seq_length = len(record.seq)
        avg_quality = sum(quality_scores) / len(quality_scores) if quality_scores else 0
        high_quality_bases = sum(1 for q in quality_scores if q >= 20)
        low_quality_bases = sum(1 for q in quality_scores if q < 20)

        # Determine pass/fail
        qc_pass = avg_quality >= 20 and high_quality_bases / seq_length >= 0.8

        result = {
            "file": ab1_file.name,
            "length": seq_length,
            "avg_quality": round(avg_quality, 2),
            "high_quality_bases": high_quality_bases,
            "low_quality_bases": low_quality_bases,
            "percent_high_quality": round(100 * high_quality_bases / seq_length, 2),
            "qc_status": "PASS" if qc_pass else "FAIL",
            "sequence": str(record.seq)
        }

        # Generate chromatogram plot if output_dir provided
        if output_dir:
            result["chromatogram_img"] = plot_chromatogram(ab1_file)

        return result
    except Exception as e:
        return {
            "file": ab1_file.name,
            "error": str(e),
            "qc_status": "ERROR"
        }

def generate_html_report(results, output_file):
    """Generate improved HTML report with collapsible chromatograms"""
    from datetime import datetime

    html = """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DNA Barcoding QC Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }}
        h1 {{ color: #333; }}
        .summary {{ background-color: white; padding: 20px; margin-bottom: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        table {{ border-collapse: collapse; width: 100%; background-color: white; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; position: sticky; top: 0; }}
        .pass {{ background-color: #d4edda; }}
        .fail {{ background-color: #f8d7da; }}
        .error {{ background-color: #fff3cd; }}
        .chromatogram-row {{ background-color: #fafafa; }}
        .chromatogram-container {{
            padding: 20px;
            overflow-x: auto;
            max-width: 100%;
        }}
        .chromatogram-container img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
            background-color: white;
        }}
        details {{ margin: 10px 0; }}
        summary {{
            cursor: pointer;
            padding: 10px;
            background-color: #e8f5e9;
            border-radius: 4px;
            font-weight: bold;
        }}
        summary:hover {{ background-color: #c8e6c9; }}
        .sequence-display {{
            font-family: 'Courier New', monospace;
            background-color: #f0f0f0;
            padding: 10px;
            margin: 10px 0;
            border-radius: 4px;
            overflow-x: auto;
            white-space: pre;
            font-size: 11px;
        }}
    </style>
</head>
<body>
    <div class="summary">
        <h1>DNA Barcoding Quality Control Report</h1>
        <p><strong>Generated:</strong> {date}</p>
        <p><strong>Note:</strong> Chromatograms show middle region (bases 50-200) where quality is typically highest.
        Click "Show Chromatogram" to view sequence traces.</p>
    </div>
    <table>
        <tr>
            <th>File</th>
            <th>Length (bp)</th>
            <th>Avg Quality</th>
            <th>High Quality %</th>
            <th>Status</th>
        </tr>
""".format(date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    for result in results:
        if "error" in result:
            row_class = "error"
            row = f"""        <tr class="{row_class}">
            <td>{result['file']}</td>
            <td colspan="3">ERROR: {result['error']}</td>
            <td>{result['qc_status']}</td>
        </tr>
"""
        else:
            row_class = "pass" if result['qc_status'] == "PASS" else "fail"
            row = f"""        <tr class="{row_class}">
            <td>{result['file']}</td>
            <td>{result['length']}</td>
            <td>{result['avg_quality']}</td>
            <td>{result['percent_high_quality']}%</td>
            <td>{result['qc_status']}</td>
        </tr>
"""
            html += row

            # Add chromatogram visualization in collapsible section
            if "chromatogram_img" in result and result["chromatogram_img"]:
                # Show first 60 bases of sequence as preview
                seq_preview = result.get('sequence', '')[:60] + '...' if len(result.get('sequence', '')) > 60 else result.get('sequence', '')

                html += f"""        <tr class="chromatogram-row">
            <td colspan="5">
                <details>
                    <summary>📊 Show Chromatogram & Sequence</summary>
                    <div class="chromatogram-container">
                        <h4>Chromatogram View (Bases 50-200)</h4>
                        <p><em>Base colors indicate quality: Green (Q≥30), Orange (Q≥20), Red (Q&lt;20)</em></p>
                        <img src="data:image/png;base64,{result['chromatogram_img']}" alt="Chromatogram">
                        <h4>Full Sequence ({result['length']} bp)</h4>
                        <div class="sequence-display">{result.get('sequence', 'N/A')}</div>
                    </div>
                </details>
            </td>
        </tr>
"""

    html += """    </table>
</body>
</html>"""

    with open(output_file, 'w') as f:
        f.write(html)

def main():
    """Main QC function"""
    if len(sys.argv) < 2:
        print("Usage: python qc_chromatograms.py <input_directory> [output_directory]")
        print("Example: python qc_chromatograms.py data/my_sequences/ results/")
        sys.exit(1)

    input_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all .ab1 files
    ab1_files = list(input_dir.glob("*.ab1"))

    if not ab1_files:
        print(f"No .ab1 files found in {input_dir}")
        sys.exit(1)

    print(f"Found {len(ab1_files)} chromatogram files")

    # Analyze each file
    results = []
    passed_sequences = []

    for ab1_file in ab1_files:
        print(f"Analyzing {ab1_file.name}...")
        result = analyze_chromatogram(ab1_file, output_dir=output_dir)
        results.append(result)

        # Save passed sequences to FASTA
        if result['qc_status'] == "PASS" and 'sequence' in result:
            passed_sequences.append({
                'id': ab1_file.stem,
                'seq': result['sequence']
            })

    # Generate outputs
    print("\nGenerating reports...")

    # HTML report
    html_file = output_dir / "qc_report.html"
    generate_html_report(results, html_file)
    print(f"HTML report: {html_file}")

    # JSON report
    json_file = output_dir / "qc_results.json"
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"JSON report: {json_file}")

    # FASTA file for passed sequences
    if passed_sequences:
        fasta_file = output_dir / "passed_sequences.fasta"
        with open(fasta_file, 'w') as f:
            for seq in passed_sequences:
                f.write(f">{seq['id']}\n{seq['seq']}\n")
        print(f"Passed sequences: {fasta_file}")

    # Summary
    passed = sum(1 for r in results if r['qc_status'] == "PASS")
    failed = sum(1 for r in results if r['qc_status'] == "FAIL")
    errors = sum(1 for r in results if r['qc_status'] == "ERROR")

    print(f"\nSummary:")
    print(f"  Total: {len(results)}")
    print(f"  Passed: {passed}")
    print(f"  Failed: {failed}")
    print(f"  Errors: {errors}")

if __name__ == "__main__":
    main()
