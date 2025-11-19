#!/usr/bin/env python3
"""
Sequence Alignment with MAFFT
Aligns sequences and generates visualization
"""

import os
import sys
import subprocess
from pathlib import Path
from Bio import SeqIO, AlignIO
import json

def run_mafft(input_fasta, output_aligned):
    """Run MAFFT alignment"""
    print(f"Running MAFFT alignment...")

    try:
        with open(output_aligned, 'w') as outfile:
            result = subprocess.run(
                ['mafft', '--auto', str(input_fasta)],
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        print(f"✓ Alignment complete: {output_aligned}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ MAFFT failed: {e.stderr}")
        return False
    except FileNotFoundError:
        print("✗ MAFFT not found. Please install MAFFT.")
        return False

def generate_alignment_stats(alignment_file):
    """Generate alignment statistics"""
    alignment = AlignIO.read(alignment_file, "fasta")

    stats = {
        "num_sequences": len(alignment),
        "alignment_length": alignment.get_alignment_length(),
        "sequences": []
    }

    for record in alignment:
        seq_stats = {
            "id": record.id,
            "length": len(str(record.seq).replace('-', '')),
            "gaps": str(record.seq).count('-'),
            "percent_gaps": round(100 * str(record.seq).count('-') / len(record.seq), 2)
        }
        stats["sequences"].append(seq_stats)

    return stats

def generate_html_report(stats, alignment_file, output_file):
    """Generate simple HTML report"""
    from datetime import datetime

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Sequence Alignment Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #333; }}
        table {{ border-collapse: collapse; width: 100%; margin-top: 20px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
        .summary {{ background-color: #f0f0f0; padding: 15px; margin: 20px 0; border-radius: 5px; }}
    </style>
</head>
<body>
    <h1>Sequence Alignment Report</h1>
    <p>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>

    <div class="summary">
        <h2>Alignment Summary</h2>
        <p><strong>Number of sequences:</strong> {stats['num_sequences']}</p>
        <p><strong>Alignment length:</strong> {stats['alignment_length']} bp</p>
        <p><strong>Alignment file:</strong> {Path(alignment_file).name}</p>
    </div>

    <h2>Sequence Details</h2>
    <table>
        <tr>
            <th>Sequence ID</th>
            <th>Original Length (bp)</th>
            <th>Gaps</th>
            <th>Gap %</th>
        </tr>
"""

    for seq in stats['sequences']:
        html += f"""        <tr>
            <td>{seq['id']}</td>
            <td>{seq['length']}</td>
            <td>{seq['gaps']}</td>
            <td>{seq['percent_gaps']}%</td>
        </tr>
"""

    html += """    </table>
</body>
</html>"""

    with open(output_file, 'w') as f:
        f.write(html)

def main():
    """Main alignment function"""
    if len(sys.argv) < 2:
        print("Usage: python align_sequences.py <input_fasta> [output_directory]")
        print("Example: python align_sequences.py results/passed_sequences.fasta results/")
        sys.exit(1)

    input_fasta = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results")
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_fasta.exists():
        print(f"Input file not found: {input_fasta}")
        sys.exit(1)

    # Check number of sequences
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    if len(sequences) < 2:
        print(f"Need at least 2 sequences for alignment. Found: {len(sequences)}")
        sys.exit(1)

    print(f"Found {len(sequences)} sequences to align")

    # Run MAFFT
    aligned_file = output_dir / "aligned_sequences.fasta"
    if not run_mafft(input_fasta, aligned_file):
        sys.exit(1)

    # Generate statistics
    print("\nGenerating alignment statistics...")
    stats = generate_alignment_stats(aligned_file)

    # Save JSON stats
    json_file = output_dir / "alignment_stats.json"
    with open(json_file, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"Statistics saved: {json_file}")

    # Generate HTML report
    html_file = output_dir / "alignment_report.html"
    generate_html_report(stats, aligned_file, html_file)
    print(f"HTML report: {html_file}")

    print(f"\nAlignment Summary:")
    print(f"  Sequences aligned: {stats['num_sequences']}")
    print(f"  Alignment length: {stats['alignment_length']} bp")
    print(f"\nOutput files:")
    print(f"  - {aligned_file}")
    print(f"  - {html_file}")

if __name__ == "__main__":
    main()
