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

def generate_visual_alignment(alignment_file, max_seqs=20, chunk_size=80):
    """Generate visual alignment with conserved/variable positions"""
    alignment = AlignIO.read(alignment_file, "fasta")

    # Limit display to max_seqs sequences
    display_seqs = alignment[:max_seqs] if len(alignment) > max_seqs else alignment
    aln_length = alignment.get_alignment_length()

    # Calculate conservation for each position
    conservation = []
    for i in range(aln_length):
        column = alignment[:, i]
        # Position is conserved if all non-gap characters are the same
        bases = [b for b in column if b != '-']
        if bases:
            most_common = max(set(bases), key=bases.count)
            conservation_rate = bases.count(most_common) / len(bases)
            conservation.append(conservation_rate)
        else:
            conservation.append(0)

    # Generate HTML alignment view in chunks
    html_chunks = []
    for chunk_start in range(0, aln_length, chunk_size):
        chunk_end = min(chunk_start + chunk_size, aln_length)

        chunk_html = f'<div class="alignment-chunk"><div class="position-header">Position {chunk_start+1}-{chunk_end}</div>'

        for record in display_seqs:
            seq_chunk = str(record.seq[chunk_start:chunk_end])

            # Format sequence with conservation styling
            formatted_seq = '<span class="sequence">'
            for i, base in enumerate(seq_chunk):
                pos = chunk_start + i
                cons = conservation[pos]

                # Conserved (>80%) = UPPERCASE, Variable (<80%) = lowercase
                if cons >= 0.8:
                    display_base = base.upper()
                    css_class = 'conserved'
                else:
                    display_base = base.lower()
                    css_class = 'variable'

                # Color code bases
                if base == '-':
                    color_class = 'gap'
                elif base.upper() == 'A':
                    color_class = 'base-a'
                elif base.upper() == 'T':
                    color_class = 'base-t'
                elif base.upper() == 'G':
                    color_class = 'base-g'
                elif base.upper() == 'C':
                    color_class = 'base-c'
                else:
                    color_class = 'base-n'

                formatted_seq += f'<span class="{css_class} {color_class}">{display_base}</span>'

            formatted_seq += '</span>'

            # Truncate long IDs
            seq_id = record.id[:30] + '...' if len(record.id) > 30 else record.id
            chunk_html += f'<div class="alignment-row"><span class="seq-id">{seq_id}</span>{formatted_seq}</div>'

        chunk_html += '</div>'
        html_chunks.append(chunk_html)

    return '\n'.join(html_chunks)

def generate_html_report(stats, alignment_file, output_file):
    """Generate comprehensive HTML report with visual alignment"""
    from datetime import datetime

    visual_alignment = generate_visual_alignment(alignment_file)

    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Sequence Alignment Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; max-width: 1400px; }}
        h1 {{ color: #333; }}
        h2 {{ color: #555; margin-top: 30px; }}
        table {{ border-collapse: collapse; width: 100%; margin-top: 20px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
        .summary {{ background-color: #f0f0f0; padding: 15px; margin: 20px 0; border-radius: 5px; }}
        .legend {{ background-color: #fff3cd; padding: 15px; margin: 20px 0; border-radius: 5px; border: 1px solid #ffc107; }}
        .legend-item {{ display: inline-block; margin-right: 20px; }}

        /* Alignment visualization */
        .alignment-chunk {{ margin: 20px 0; padding: 10px; background: #f8f9fa; border-radius: 5px; overflow-x: auto; }}
        .position-header {{ font-weight: bold; color: #666; margin-bottom: 10px; }}
        .alignment-row {{ font-family: 'Courier New', monospace; margin: 2px 0; white-space: nowrap; }}
        .seq-id {{ display: inline-block; width: 200px; padding-right: 10px; font-weight: bold; color: #333; }}
        .sequence {{ letter-spacing: 1px; }}

        /* Conservation styling */
        .conserved {{ font-weight: bold; font-size: 1.1em; }}
        .variable {{ font-size: 0.9em; opacity: 0.8; }}

        /* Base colors */
        .base-a {{ color: #00C851; }}
        .base-t {{ color: #ff4444; }}
        .base-g {{ color: #ffbb33; }}
        .base-c {{ color: #33b5e5; }}
        .base-n {{ color: #999; }}
        .gap {{ color: #ddd; }}
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

    <div class="legend">
        <h3>Legend</h3>
        <div class="legend-item"><strong>UPPERCASE</strong> = Conserved position (≥80% identity)</div>
        <div class="legend-item"><strong>lowercase</strong> = Variable position (&lt;80% identity)</div>
        <br>
        <div class="legend-item"><span class="base-a">A</span> = Adenine</div>
        <div class="legend-item"><span class="base-t">T</span> = Thymine</div>
        <div class="legend-item"><span class="base-g">G</span> = Guanine</div>
        <div class="legend-item"><span class="base-c">C</span> = Cytosine</div>
        <div class="legend-item"><span class="gap">-</span> = Gap</div>
    </div>

    <h2>Visual Alignment</h2>
    <p><em>Scroll right to see full alignment. Conserved positions (UPPERCASE) vs variable positions (lowercase).</em></p>
    {visual_alignment}

    <h2>Sequence Statistics</h2>
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
