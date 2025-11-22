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

def generate_visual_alignment(alignment_file, max_seqs=20, chunk_size=100):
    """Generate heatmap-style visual alignment with grayscale conservation boxes"""
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

    # Generate HTML alignment view in chunks with heatmap boxes
    html_chunks = []
    for chunk_start in range(0, aln_length, chunk_size):
        chunk_end = min(chunk_start + chunk_size, aln_length)

        chunk_html = f'<div class="alignment-chunk"><div class="position-header">Position {chunk_start+1}-{chunk_end}</div>'

        for record in display_seqs:
            seq_chunk = str(record.seq[chunk_start:chunk_end])

            # Format sequence with heatmap boxes
            formatted_seq = '<div class="sequence">'
            for i, base in enumerate(seq_chunk):
                pos = chunk_start + i
                cons = conservation[pos]

                # Determine conservation class based on rate
                if base == '-':
                    box_class = 'gap-box'
                elif cons >= 0.95:
                    box_class = 'cons-100'
                elif cons >= 0.90:
                    box_class = 'cons-90'
                elif cons >= 0.80:
                    box_class = 'cons-80'
                elif cons >= 0.70:
                    box_class = 'cons-70'
                elif cons >= 0.60:
                    box_class = 'cons-60'
                elif cons >= 0.50:
                    box_class = 'cons-50'
                elif cons >= 0.40:
                    box_class = 'cons-40'
                elif cons >= 0.30:
                    box_class = 'cons-30'
                elif cons >= 0.20:
                    box_class = 'cons-20'
                else:
                    box_class = 'cons-10'

                formatted_seq += f'<span class="base-box {box_class}" title="{base.upper()} (conservation: {cons*100:.0f}%)">{base.upper()}</span>'

            formatted_seq += '</div>'

            # Truncate long IDs
            seq_id = record.id[:28] + '...' if len(record.id) > 28 else record.id
            chunk_html += f'<div class="alignment-row"><span class="seq-id">{seq_id}</span>{formatted_seq}</div>'

        chunk_html += '</div>'
        html_chunks.append(chunk_html)

    return '\n'.join(html_chunks)

def generate_html_report(stats, alignment_file, output_file):
    """Generate comprehensive HTML report with visual alignment using new design system"""
    from datetime import datetime
    from pathlib import Path

    visual_alignment = generate_visual_alignment(alignment_file)

    # Calculate metrics
    avg_gaps = sum(s['percent_gaps'] for s in stats['sequences']) / len(stats['sequences']) if stats['sequences'] else 0

    # Read CSS files and embed them
    project_root = Path(__file__).parent.parent.parent
    base_css = (project_root / "tracking/styles/base.css").read_text()
    components_css = (project_root / "tracking/styles/components.css").read_text()
    reports_css = (project_root / "tracking/styles/reports.css").read_text()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Sequence Alignment Report - DNA Barcoding</title>

    <!-- Embedded CSS for reliable loading -->
    <style>
{base_css}

{components_css}

{reports_css}

        /* Alignment-specific styles - Heatmap visualization */
        .alignment-chunk {{
            margin: 20px 0;
            padding: 15px;
            background: white;
            border-radius: var(--radius-md);
            overflow-x: auto;
            box-shadow: var(--shadow-md);
        }}

        .position-header {{
            font-weight: 600;
            color: var(--text-secondary);
            margin-bottom: 10px;
            font-family: var(--font-mono);
            font-size: 0.9rem;
        }}

        .alignment-row {{
            display: flex;
            align-items: center;
            margin: 1px 0;
            white-space: nowrap;
        }}

        .seq-id {{
            display: inline-block;
            width: 200px;
            padding-right: 15px;
            font-weight: 600;
            color: var(--text-primary);
            font-family: var(--font-primary);
            font-size: 0.85rem;
            flex-shrink: 0;
        }}

        .sequence {{
            display: flex;
            gap: 1px;
            font-family: var(--font-mono);
            font-size: 0.75rem;
        }}

        /* Heatmap box styling - darker = more conserved */
        .base-box {{
            width: 14px;
            height: 18px;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            border-radius: 2px;
            font-weight: 600;
            color: white;
            text-shadow: 0 1px 1px rgba(0, 0, 0, 0.3);
            transition: transform 0.1s ease;
        }}

        .base-box:hover {{
            transform: scale(1.5);
            z-index: 10;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.3);
        }}

        /* Conservation levels - grayscale heatmap */
        /* Highly conserved (90-100%) - dark gray/black */
        .cons-100 {{ background: #1a1a1a; }}
        .cons-90 {{ background: #2d2d2d; }}

        /* Moderately conserved (70-89%) - medium gray */
        .cons-80 {{ background: #4a4a4a; }}
        .cons-70 {{ background: #6b6b6b; }}

        /* Low conservation (50-69%) - light gray */
        .cons-60 {{ background: #8c8c8c; }}
        .cons-50 {{ background: #adadad; }}

        /* Variable (<50%) - very light gray */
        .cons-40 {{ background: #c9c9c9; color: #555; text-shadow: none; }}
        .cons-30 {{ background: #dedede; color: #555; text-shadow: none; }}
        .cons-20 {{ background: #efefef; color: #555; text-shadow: none; }}
        .cons-10 {{ background: #f5f5f5; color: #777; text-shadow: none; }}

        /* Gaps */
        .gap-box {{
            background: #ffd4d4;
            color: #999;
            text-shadow: none;
            border: 1px solid #ffb3b3;
        }}

        /* Legend */
        .conservation-legend {{
            margin: 20px 0;
            padding: 15px;
            background: var(--bg-purple-light);
            border-radius: var(--radius-md);
            border-left: 4px solid var(--purple);
        }}

        .legend-scale {{
            display: flex;
            gap: 2px;
            margin: 10px 0;
            align-items: center;
        }}

        .legend-box {{
            width: 40px;
            height: 25px;
            display: flex;
            align-items: center;
            justify-content: center;
            border-radius: 3px;
            font-size: 0.7rem;
            font-weight: 600;
        }}

        .legend-label {{
            font-size: 0.85rem;
            color: var(--text-secondary);
            margin: 0 10px;
        }}
    </style>
</head>
<body>
    <!-- Report Header -->
    <header class="report-header">
        <h1>🧬 Sequence Alignment Report</h1>
        <div class="progress-badge">Step 3 of 5</div>
        <div class="report-date">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</div>
    </header>

    <!-- Summary Dashboard -->
    <section class="summary-dashboard">
        <div class="metric-card metric-primary">
            <div class="metric-value">{stats['num_sequences']}</div>
            <div class="metric-label">Sequences Aligned</div>
            <div class="metric-icon">🧬</div>
        </div>

        <div class="metric-card metric-success">
            <div class="metric-value">{stats['alignment_length']}</div>
            <div class="metric-label">Alignment Length (bp)</div>
            <div class="metric-icon">📏</div>
        </div>

        <div class="metric-card metric-primary">
            <div class="metric-value">{avg_gaps:.1f}%</div>
            <div class="metric-label">Avg. Gap %</div>
            <div class="metric-icon">📊</div>
        </div>

        <div class="metric-card metric-success">
            <div class="metric-value">MAFFT</div>
            <div class="metric-label">Alignment Tool</div>
            <div class="metric-icon">⚙️</div>
        </div>
    </section>

    <!-- Main Content -->
    <main class="report-content">
        <div class="content-section">
            <h2>Multiple Sequence Alignment</h2>
            <p>Generated using MAFFT --auto algorithm</p>

            <div class="info-box info-tip">
                <strong>💡 Reading the Heatmap Alignment:</strong>
                <ul>
                    <li><strong>Dark boxes (black/dark gray):</strong> Highly conserved positions (90-100% identical across sequences)</li>
                    <li><strong>Medium gray boxes:</strong> Moderately conserved (50-90% identical)</li>
                    <li><strong>Light gray boxes:</strong> Variable positions (&lt;50% identical)</li>
                    <li><strong>Pink boxes:</strong> Gaps (-) representing insertions/deletions</li>
                    <li><strong>Hover over any box:</strong> See the base letter and conservation %</li>
                </ul>
            </div>

            <div class="conservation-legend">
                <h4 style="margin-top: 0;">Conservation Scale</h4>
                <div class="legend-scale">
                    <span class="legend-label">High</span>
                    <span class="legend-box cons-100" style="color: white;">100%</span>
                    <span class="legend-box cons-90" style="color: white;">90%</span>
                    <span class="legend-box cons-80" style="color: white;">80%</span>
                    <span class="legend-box cons-70" style="color: white;">70%</span>
                    <span class="legend-box cons-60" style="color: white;">60%</span>
                    <span class="legend-box cons-50" style="color: white;">50%</span>
                    <span class="legend-box cons-40" style="color: #555;">40%</span>
                    <span class="legend-box cons-30" style="color: #555;">30%</span>
                    <span class="legend-box cons-20" style="color: #555;">20%</span>
                    <span class="legend-box cons-10" style="color: #777;">10%</span>
                    <span class="legend-label">Low</span>
                    <span class="legend-box gap-box" style="color: #999; border: 1px solid #ffb3b3;">Gap</span>
                </div>
                <p style="font-size: 0.85rem; margin-top: 10px; margin-bottom: 0;">
                    <strong>Interpretation:</strong> Darker colors indicate higher sequence similarity at that position.
                    DNA barcoding regions (like COI) should show mostly dark boxes (high conservation within species)
                    with some lighter boxes showing intraspecific variation.
                </p>
            </div>

            <h3>Visual Alignment</h3>
            <p><em>Scroll horizontally to view the full alignment. Hover over boxes to see base and conservation %. Alignment shown in 100bp chunks.</em></p>

            {visual_alignment}

            <h3 style="margin-top: 2rem;">Sequence Statistics</h3>
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Sequence ID</th>
                        <th>Original Length (bp)</th>
                        <th>Gaps Added</th>
                        <th>Gap %</th>
                    </tr>
                </thead>
                <tbody>
"""

    for seq in stats['sequences']:
        row_class = "row-pass" if seq['percent_gaps'] < 5 else "row-warning" if seq['percent_gaps'] < 10 else "row-fail"
        html += f"""                    <tr class="{row_class}">
                        <td><code>{seq['id']}</code></td>
                        <td>{seq['length']}</td>
                        <td>{seq['gaps']}</td>
                        <td><span class="badge badge-{"success" if seq['percent_gaps'] < 5 else "warning" if seq['percent_gaps'] < 10 else "fail"}">{seq['percent_gaps']}%</span></td>
                    </tr>
"""

    html += """                </tbody>
            </table>
        </div>
    </main>

    <!-- Footer with Help -->
    <footer class="report-footer">
        <div class="help-section">
            <h3>Understanding Sequence Alignment</h3>

            <h4>What is multiple sequence alignment?</h4>
            <p>Multiple sequence alignment (MSA) arranges DNA sequences to identify regions of similarity. These similarities may indicate functional, structural, or evolutionary relationships between sequences.</p>

            <h4>Why align sequences?</h4>
            <ul>
                <li><strong>Phylogenetic analysis:</strong> Alignments are required for building evolutionary trees</li>
                <li><strong>Identify conserved regions:</strong> Find functionally important DNA regions shared across species</li>
                <li><strong>Compare sequences:</strong> See how your samples relate to reference sequences</li>
            </ul>

            <h4>What is MAFFT?</h4>
            <p>MAFFT is a state-of-the-art multiple sequence alignment program. The --auto option automatically selects the most appropriate algorithm based on dataset size and similarity.</p>

            <h4>What do gaps mean?</h4>
            <p>Gaps (-) represent insertions or deletions (indels) that occurred during evolution. MAFFT inserts gaps to maximize alignment quality. High gap percentages (&gt;10%) may indicate:</p>
            <ul>
                <li>Sequence quality issues</li>
                <li>Distant evolutionary relationships</li>
                <li>Sequencing errors or contamination</li>
            </ul>

            <h4>Conserved vs Variable Positions</h4>
            <p>Conserved positions (UPPERCASE) show ≥80% identity across all sequences, indicating functionally important or slowly-evolving regions. Variable positions (lowercase) may represent neutral mutations or rapidly-evolving regions.</p>
        </div>
    </footer>
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
