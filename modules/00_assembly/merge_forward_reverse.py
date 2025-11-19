#!/usr/bin/env python3
"""
Merge Forward and Reverse Sanger Sequences
Pairs F/R sequences and generates consensus
"""

import os
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
import re

def pair_sequences(fasta_file):
    """Pair forward and reverse sequences by sample name"""
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    pairs = {}
    unpaired = []

    for record in sequences:
        # Extract sample name (everything before _F or _R)
        match = re.match(r'(.+?)_([FR])(?:_|$)', record.id)
        if match:
            sample_name = match.group(1)
            direction = match.group(2)

            if sample_name not in pairs:
                pairs[sample_name] = {}

            pairs[sample_name][direction] = record
        else:
            unpaired.append(record.id)

    return pairs, unpaired

def merge_sequences(forward_seq, reverse_seq, min_overlap=20):
    """Merge forward and reverse sequences

    For now, uses simple concatenation after reverse-complementing R sequence.
    In production, would use proper overlap detection and consensus calling.
    """
    # Reverse complement the reverse sequence
    rev_rc = reverse_seq.reverse_complement()

    # Simple merge: take forward sequence + reverse complement
    # In production, would detect overlap and create consensus
    merged_seq = str(forward_seq.seq) + str(rev_rc.seq)

    return merged_seq

def generate_consensus(pairs):
    """Generate consensus sequences from F/R pairs"""
    consensus_records = []
    stats = {
        "total_pairs": len(pairs),
        "successful_merges": 0,
        "forward_only": 0,
        "reverse_only": 0,
        "details": []
    }

    for sample_name, seqs in pairs.items():
        detail = {"sample": sample_name}

        if 'F' in seqs and 'R' in seqs:
            # Have both F and R
            merged = merge_sequences(seqs['F'], seqs['R'])
            consensus = SeqRecord(
                Seq(merged),
                id=sample_name,
                description=f"consensus from {seqs['F'].id} and {seqs['R'].id}"
            )
            consensus_records.append(consensus)
            stats["successful_merges"] += 1
            detail["status"] = "merged"
            detail["length"] = len(merged)

        elif 'F' in seqs:
            # Only forward
            consensus_records.append(SeqRecord(
                seqs['F'].seq,
                id=sample_name,
                description=f"forward only from {seqs['F'].id}"
            ))
            stats["forward_only"] += 1
            detail["status"] = "forward_only"
            detail["length"] = len(seqs['F'].seq)

        elif 'R' in seqs:
            # Only reverse - reverse complement it
            rev_rc = seqs['R'].reverse_complement()
            consensus_records.append(SeqRecord(
                rev_rc.seq,
                id=sample_name,
                description=f"reverse only (RC) from {seqs['R'].id}"
            ))
            stats["reverse_only"] += 1
            detail["status"] = "reverse_only"
            detail["length"] = len(rev_rc.seq)

        stats["details"].append(detail)

    return consensus_records, stats

def generate_html_report(stats, output_file):
    """Generate simple HTML report"""
    from datetime import datetime

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Sequence Assembly Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #333; }}
        .summary {{ background-color: #f0f0f0; padding: 15px; margin: 20px 0; border-radius: 5px; }}
        table {{ border-collapse: collapse; width: 100%; margin-top: 20px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
        .merged {{ background-color: #d4edda; }}
        .forward-only {{ background-color: #fff3cd; }}
        .reverse-only {{ background-color: #fff3cd; }}
    </style>
</head>
<body>
    <h1>Sequence Assembly Report</h1>
    <p>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>

    <div class="summary">
        <h2>Assembly Summary</h2>
        <p><strong>Total pairs:</strong> {stats['total_pairs']}</p>
        <p><strong>Successfully merged:</strong> {stats['successful_merges']}</p>
        <p><strong>Forward only:</strong> {stats['forward_only']}</p>
        <p><strong>Reverse only:</strong> {stats['reverse_only']}</p>
    </div>

    <h2>Sample Details</h2>
    <table>
        <tr>
            <th>Sample</th>
            <th>Status</th>
            <th>Consensus Length (bp)</th>
        </tr>
"""

    for detail in stats['details']:
        status_class = detail['status'].replace('_', '-')
        html += f"""        <tr class="{status_class}">
            <td>{detail['sample']}</td>
            <td>{detail['status'].replace('_', ' ').title()}</td>
            <td>{detail['length']}</td>
        </tr>
"""

    html += """    </table>
</body>
</html>"""

    with open(output_file, 'w') as f:
        f.write(html)

def main():
    """Main assembly function"""
    if len(sys.argv) < 2:
        print("Usage: python merge_forward_reverse.py <input_fasta> [output_directory]")
        print("Example: python merge_forward_reverse.py results/passed_sequences.fasta results/")
        print("\nInput FASTA should contain sequences with _F and _R suffixes")
        sys.exit(1)

    input_fasta = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results")
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_fasta.exists():
        print(f"Input file not found: {input_fasta}")
        sys.exit(1)

    print(f"Reading sequences from: {input_fasta}")

    # Pair sequences
    pairs, unpaired = pair_sequences(input_fasta)

    print(f"\nFound {len(pairs)} samples:")
    for sample, seqs in pairs.items():
        dirs = []
        if 'F' in seqs:
            dirs.append('F')
        if 'R' in seqs:
            dirs.append('R')
        print(f"  {sample}: {'+'.join(dirs)}")

    if unpaired:
        print(f"\nUnpaired sequences (skipped): {', '.join(unpaired)}")

    # Generate consensus
    print("\nGenerating consensus sequences...")
    consensus_records, stats = generate_consensus(pairs)

    # Save consensus FASTA
    consensus_file = output_dir / "consensus_sequences.fasta"
    SeqIO.write(consensus_records, consensus_file, "fasta")
    print(f"\nConsensus sequences: {consensus_file}")

    # Save JSON stats
    json_file = output_dir / "assembly_stats.json"
    with open(json_file, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"Statistics: {json_file}")

    # Generate HTML report
    html_file = output_dir / "assembly_report.html"
    generate_html_report(stats, html_file)
    print(f"HTML report: {html_file}")

    print(f"\nAssembly Summary:")
    print(f"  Total samples: {stats['total_pairs']}")
    print(f"  Merged (F+R): {stats['successful_merges']}")
    print(f"  Forward only: {stats['forward_only']}")
    print(f"  Reverse only: {stats['reverse_only']}")

if __name__ == "__main__":
    main()
