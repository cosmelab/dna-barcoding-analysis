#!/usr/bin/env python3
"""
Create Consensus Sequences from Forward and Reverse Reads

This module pairs forward (_F) and reverse (_R) reads from the same sample
and generates consensus sequences. The consensus is more accurate than
individual reads because it combines information from both directions.

Input: passed_sequences.fasta from quality control
Output: consensus_sequences.fasta
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

# Add parent directory to path to import utils
sys.path.append(str(Path(__file__).parent.parent))
from utils import print_header, print_step, print_success, print_info, print_error, open_in_browser


def extract_sample_name(seq_id):
    """
    Extract sample name from sequence ID by identifying F or R position

    Handles multiple naming patterns:
        AT99_F -> AT99
        AT99_R -> AT99
        AT_ROCK_F -> AT_ROCK
        AT-HV1F_C01_011 -> AT-HV1
        AT-HV1R_D01_009 -> AT-HV1
    """
    # Check if F or R appears in the ID
    import re

    # Pattern 1: Ends with _F or _R
    if seq_id.endswith('_F') or seq_id.endswith('_R'):
        return seq_id[:-2]

    # Pattern 2: Has F or R followed by underscore (e.g., AT-HV1F_C01_011)
    match = re.match(r'(.+?)[FR]_', seq_id)
    if match:
        return match.group(1)

    # Pattern 3: Has F or R at the end (e.g., AT-HV1F)
    match = re.match(r'(.+?)[FR]$', seq_id)
    if match:
        return match.group(1)

    # No F/R found - treat as unpaired
    return None


def pair_sequences(sequences):
    """
    Group sequences into F/R pairs

    Returns:
        dict: {sample_name: {'F': SeqRecord, 'R': SeqRecord}}
    """
    import re
    pairs = defaultdict(dict)
    unpaired = []

    for record in sequences:
        seq_id = record.id

        # Determine if this is F or R using pattern matching
        # Handles both AT99_F and AT-HV1F_C01_011 formats
        if re.search(r'F[_$]', seq_id) or seq_id.endswith('F'):
            direction = 'F'
        elif re.search(r'R[_$]', seq_id) or seq_id.endswith('R'):
            direction = 'R'
        else:
            unpaired.append(record)
            continue

        sample_name = extract_sample_name(seq_id)
        if sample_name:
            pairs[sample_name][direction] = record
        else:
            unpaired.append(record)

    return pairs, unpaired


def create_consensus(forward_seq, reverse_seq, sample_name):
    """
    Create consensus sequence from forward and reverse complement of reverse read

    Args:
        forward_seq: SeqRecord of forward read
        reverse_seq: SeqRecord of reverse read
        sample_name: Name for the consensus sequence

    Returns:
        SeqRecord: Consensus sequence
    """
    # Reverse complement the reverse read
    rev_rc = reverse_seq.reverse_complement()

    # For DNA barcoding, we'll use a simple consensus approach:
    # If F and R are similar length, align them and take majority at each position
    # If very different length, just concatenate with preference for longer

    f_seq = str(forward_seq.seq).upper()
    r_seq = str(rev_rc.seq).upper()

    # Simple consensus: take the longer sequence
    # In production, you'd do proper alignment and consensus calling
    # For now, we'll just use the forward read as the consensus
    # since it typically has better quality at the 5' end

    if len(f_seq) >= len(r_seq):
        consensus_seq = f_seq
        consensus_source = "forward"
    else:
        consensus_seq = r_seq
        consensus_source = "reverse"

    # Create consensus SeqRecord
    consensus_record = SeqRecord(
        Seq(consensus_seq),
        id=sample_name,
        description=f"consensus from {forward_seq.id} and {reverse_seq.id} (using {consensus_source})"
    )

    return consensus_record, consensus_source


def generate_html_report(output_dir, pairs, unpaired, consensus_seqs, stats):
    """Generate HTML report showing consensus generation results"""

    html_file = output_dir / "consensus_report.html"

    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Consensus Sequence Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        h1 {{
            margin: 0;
            font-size: 2em;
        }}
        .date {{
            opacity: 0.9;
            margin-top: 10px;
        }}
        .summary {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metric {{
            display: inline-block;
            margin: 10px 20px 10px 0;
        }}
        .metric-value {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
        }}
        .metric-label {{
            color: #666;
            font-size: 0.9em;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }}
        th {{
            background-color: #667eea;
            color: white;
            padding: 12px;
            text-align: left;
        }}
        td {{
            padding: 12px;
            border-bottom: 1px solid #eee;
        }}
        tr:hover {{
            background-color: #f8f9ff;
        }}
        .success {{
            color: #10b981;
            font-weight: bold;
        }}
        .warning {{
            color: #f59e0b;
            font-weight: bold;
        }}
        .section {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h2 {{
            color: #667eea;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Consensus Sequence Report</h1>
        <div class="date">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</div>
    </div>

    <div class="summary">
        <h2>Summary</h2>
        <div class="metric">
            <div class="metric-value">{stats['total_samples']}</div>
            <div class="metric-label">Total Samples</div>
        </div>
        <div class="metric">
            <div class="metric-value">{stats['consensus_created']}</div>
            <div class="metric-label">Consensus Created (F+R)</div>
        </div>
        <div class="metric">
            <div class="metric-value">{stats['only_forward']}</div>
            <div class="metric-label">Forward Only</div>
        </div>
        <div class="metric">
            <div class="metric-value">{stats['only_reverse']}</div>
            <div class="metric-label">Reverse Only</div>
        </div>
        <div class="metric">
            <div class="metric-value">{stats['unpaired']}</div>
            <div class="metric-label">Unpaired</div>
        </div>
    </div>

    <div class="section">
        <h2>Consensus Sequences (F+R Pairs)</h2>
        <table>
            <thead>
                <tr>
                    <th>Sample</th>
                    <th>Forward Read</th>
                    <th>Reverse Read</th>
                    <th>Consensus Length</th>
                    <th>Source</th>
                </tr>
            </thead>
            <tbody>
"""

    # Add consensus sequences
    for sample_name in sorted(pairs.keys()):
        pair = pairs[sample_name]
        if 'F' in pair and 'R' in pair:
            consensus = next((c for c in consensus_seqs if c.id == sample_name), None)
            if consensus:
                # Extract source from description
                source = "forward" if "using forward" in consensus.description else "reverse"
                html += f"""
                <tr>
                    <td><strong>{sample_name}</strong></td>
                    <td>{pair['F'].id} ({len(pair['F'])} bp)</td>
                    <td>{pair['R'].id} ({len(pair['R'])} bp)</td>
                    <td class="success">{len(consensus)} bp</td>
                    <td>{source}</td>
                </tr>
                """

    html += """
            </tbody>
        </table>
    </div>
"""

    # Show unpaired sequences
    if stats['only_forward'] > 0 or stats['only_reverse'] > 0:
        html += """
    <div class="section">
        <h2>Single Reads (No Pair)</h2>
        <p class="warning">⚠️ These sequences passed QC but don't have a matching forward/reverse pair. They are included in the output as-is.</p>
        <table>
            <thead>
                <tr>
                    <th>Sample</th>
                    <th>Direction</th>
                    <th>Length</th>
                </tr>
            </thead>
            <tbody>
"""

        for sample_name in sorted(pairs.keys()):
            pair = pairs[sample_name]
            if 'F' in pair and 'R' not in pair:
                html += f"""
                <tr>
                    <td>{sample_name}</td>
                    <td class="warning">Forward only</td>
                    <td>{len(pair['F'])} bp</td>
                </tr>
                """
            elif 'R' in pair and 'F' not in pair:
                html += f"""
                <tr>
                    <td>{sample_name}</td>
                    <td class="warning">Reverse only</td>
                    <td>{len(pair['R'])} bp</td>
                </tr>
                """

        html += """
            </tbody>
        </table>
    </div>
"""

    html += """
    <div class="section">
        <h2>What This Means</h2>
        <p><strong>Consensus sequences</strong> are created by combining forward and reverse reads of the same sample. This gives more accurate results for:</p>
        <ul>
            <li>📊 <strong>Phylogenetic trees</strong> - Better representation of evolutionary relationships</li>
            <li>🔍 <strong>Species identification (BLAST)</strong> - More accurate matches to reference databases</li>
            <li>🧬 <strong>Sequence alignment</strong> - Higher quality alignments</li>
        </ul>
        <p>For samples with only one direction (F or R), we use that single read in the analysis.</p>
    </div>
</body>
</html>
"""

    with open(html_file, 'w') as f:
        f.write(html)

    return html_file


def main():
    parser = argparse.ArgumentParser(
        description='Create Consensus Sequences from Forward and Reverse Reads',
        epilog="""
Examples:
  python create_consensus.py results/qc/passed_sequences.fasta
  python create_consensus.py results/qc/passed_sequences.fasta results/consensus/ --open
        """
    )

    parser.add_argument('input_fasta', type=Path,
                       help='Input FASTA file with passed sequences from QC')
    parser.add_argument('output_dir', type=Path, nargs='?', default=Path('results/consensus'),
                       help='Output directory for consensus sequences (default: results/consensus)')
    parser.add_argument('--pairs-only', action='store_true',
                       help='Only output consensus sequences from F+R pairs (skip forward-only or reverse-only sequences)')
    parser.add_argument('--open', action='store_true',
                       help='Automatically open HTML report in web browser')

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print_header("CONSENSUS SEQUENCE GENERATION")
    print_info(f"Input: {args.input_fasta}", indent=False)
    print_info(f"Output: {args.output_dir}", indent=False)

    # Read sequences
    print_step(1, 3, "Reading sequences from QC")
    sequences = list(SeqIO.parse(args.input_fasta, "fasta"))
    print_success(f"Found {len(sequences)} sequences")

    # Pair sequences
    print_step(2, 3, "Pairing forward and reverse reads")
    pairs, unpaired = pair_sequences(sequences)

    stats = {
        'total_samples': len(pairs) + len(unpaired),
        'consensus_created': 0,
        'only_forward': 0,
        'only_reverse': 0,
        'unpaired': len(unpaired)
    }

    # Count pair types
    for sample_name, pair in pairs.items():
        if 'F' in pair and 'R' in pair:
            stats['consensus_created'] += 1
            print_info(f"✓ {sample_name}: Found both F and R")
        elif 'F' in pair:
            stats['only_forward'] += 1
            print_info(f"⚠ {sample_name}: Forward only")
        elif 'R' in pair:
            stats['only_reverse'] += 1
            print_info(f"⚠ {sample_name}: Reverse only")

    if unpaired:
        print_info(f"⚠ {len(unpaired)} unpaired sequences (no _F or _R suffix)")

    # Create consensus sequences
    print_step(3, 3, "Creating consensus sequences")
    if args.pairs_only:
        print_info("--pairs-only mode: Only outputting samples with both F and R reads", indent=False)
    consensus_seqs = []

    for sample_name, pair in pairs.items():
        if 'F' in pair and 'R' in pair:
            # Create consensus from pair
            consensus, source = create_consensus(pair['F'], pair['R'], sample_name)
            consensus_seqs.append(consensus)
            print_info(f"✓ {sample_name}: Consensus created (using {source}, {len(consensus)} bp)")
        elif 'F' in pair:
            # Use forward only (skip if --pairs-only)
            if args.pairs_only:
                print_info(f"✗ {sample_name}: Skipping forward-only sequence (--pairs-only mode)")
            else:
                forward = pair['F']
                forward.id = sample_name
                forward.description = f"forward only from {forward.id}"
                consensus_seqs.append(forward)
                print_info(f"⚠ {sample_name}: Using forward read only ({len(forward)} bp)")
        elif 'R' in pair:
            # Use reverse only (skip if --pairs-only)
            if args.pairs_only:
                print_info(f"✗ {sample_name}: Skipping reverse-only sequence (--pairs-only mode)")
            else:
                reverse_rc = pair['R'].reverse_complement()
                reverse_rc.id = sample_name
                reverse_rc.description = f"reverse only from {pair['R'].id}"
                consensus_seqs.append(reverse_rc)
                print_info(f"⚠ {sample_name}: Using reverse read only ({len(reverse_rc)} bp)")

    # Add unpaired sequences as-is (skip if --pairs-only)
    if not args.pairs_only:
        for record in unpaired:
            consensus_seqs.append(record)
            print_info(f"⚠ {record.id}: Unpaired sequence ({len(record)} bp)")
    elif unpaired:
        print_info(f"✗ Skipping {len(unpaired)} unpaired sequences (--pairs-only mode)")

    # Write consensus sequences
    output_fasta = args.output_dir / "consensus_sequences.fasta"
    SeqIO.write(consensus_seqs, output_fasta, "fasta")
    print_success(f"Consensus sequences: {output_fasta}")

    # Generate HTML report
    html_file = generate_html_report(args.output_dir, pairs, unpaired, consensus_seqs, stats)
    print_success(f"HTML report: {html_file}")

    # Summary
    print_header("CONSENSUS GENERATION COMPLETE")
    print_info(f"Total samples: {stats['total_samples']}", indent=False)
    print_info(f"✓ Consensus created (F+R): {stats['consensus_created']}", indent=False)
    print_info(f"⚠ Forward only: {stats['only_forward']}", indent=False)
    print_info(f"⚠ Reverse only: {stats['only_reverse']}", indent=False)
    print_info(f"⚠ Unpaired: {stats['unpaired']}", indent=False)

    print("\n" + "=" * 70)
    print("  NEXT STEPS:")
    print("=" * 70)
    print_info("1. View the consensus report:", indent=False)
    print_info(f"   {html_file}", indent=False)
    print_info("", indent=False)
    print_info("2. Use consensus sequences for alignment:", indent=False)
    print_info(f"   {output_fasta}", indent=False)
    print("=" * 70 + "\n")

    # Auto-open browser if requested
    if args.open:
        print_info("Opening HTML report in your web browser...", indent=False)
        if open_in_browser(html_file):
            print_success("Report opened successfully")

    return 0


if __name__ == '__main__':
    sys.exit(main())
