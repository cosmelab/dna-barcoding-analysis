#!/usr/bin/env python3
"""
Quality Control for Sanger Chromatograms
Analyzes .ab1 files and generates QC report
"""

import os
import sys
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import json

def analyze_chromatogram(ab1_file):
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

        return {
            "file": ab1_file.name,
            "length": seq_length,
            "avg_quality": round(avg_quality, 2),
            "high_quality_bases": high_quality_bases,
            "low_quality_bases": low_quality_bases,
            "percent_high_quality": round(100 * high_quality_bases / seq_length, 2),
            "qc_status": "PASS" if qc_pass else "FAIL",
            "sequence": str(record.seq)
        }
    except Exception as e:
        return {
            "file": ab1_file.name,
            "error": str(e),
            "qc_status": "ERROR"
        }

def generate_html_report(results, output_file):
    """Generate simple HTML report"""
    html = """<!DOCTYPE html>
<html>
<head>
    <title>DNA Barcoding QC Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #333; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #4CAF50; color: white; }
        .pass { background-color: #d4edda; }
        .fail { background-color: #f8d7da; }
        .error { background-color: #fff3cd; }
    </style>
</head>
<body>
    <h1>Quality Control Report</h1>
    <p>Generated: {date}</p>
    <table>
        <tr>
            <th>File</th>
            <th>Length (bp)</th>
            <th>Avg Quality</th>
            <th>High Quality %</th>
            <th>Status</th>
        </tr>
"""

    from datetime import datetime
    html = html.format(date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

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
    output_dir.mkdir(exist_ok=True)

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
        result = analyze_chromatogram(ab1_file)
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
