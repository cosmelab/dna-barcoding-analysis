#!/usr/bin/env python3
"""
Species Identification using BLAST
Queries sequences against NCBI GenBank database
"""

import os
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import time
import json

def blast_sequence(sequence, seq_id):
    """BLAST a single sequence against GenBank"""
    print(f"  BLASTing {seq_id}...", end=" ")

    try:
        # Run BLAST (this takes time - it's querying NCBI)
        result_handle = NCBIWWW.qblast(
            "blastn",
            "nt",
            sequence,
            hitlist_size=5,
            megablast=True
        )

        # Parse results
        blast_record = NCBIXML.read(result_handle)

        hits = []
        for alignment in blast_record.alignments[:5]:
            for hsp in alignment.hsps:
                hits.append({
                    "accession": alignment.accession,
                    "definition": alignment.title,
                    "length": alignment.length,
                    "e_value": hsp.expect,
                    "identity": hsp.identities,
                    "alignment_length": hsp.align_length,
                    "percent_identity": round(100 * hsp.identities / hsp.align_length, 2)
                })
                break  # Only get first HSP per alignment

        print(f"Done ({len(hits)} hits)")
        return {
            "sequence_id": seq_id,
            "status": "SUCCESS",
            "top_hit": hits[0] if hits else None,
            "all_hits": hits
        }

    except Exception as e:
        print(f"Error: {str(e)}")
        return {
            "sequence_id": seq_id,
            "status": "ERROR",
            "error": str(e)
        }

def generate_html_report(results, output_file):
    """Generate simple HTML report with BLAST results"""
    html = """<!DOCTYPE html>
<html>
<head>
    <title>Species Identification Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #333; }}
        h2 {{ color: #666; margin-top: 30px; }}
        table {{ border-collapse: collapse; width: 100%; margin-top: 10px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
        .success {{ background-color: #d4edda; }}
        .error {{ background-color: #f8d7da; }}
        .match-high {{ background-color: #d4edda; }}
        .match-medium {{ background-color: #fff3cd; }}
        .match-low {{ background-color: #f8d7da; }}
    </style>
</head>
<body>
    <h1>Species Identification Report</h1>
    <p>BLAST search against NCBI GenBank database</p>
    <p>Generated: {date}</p>
"""

    from datetime import datetime
    html = html.format(date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    for result in results:
        seq_id = result['sequence_id']
        html += f"\n    <h2>Sequence: {seq_id}</h2>\n"

        if result['status'] == "ERROR":
            html += f'    <p class="error">ERROR: {result["error"]}</p>\n'
            continue

        if not result['top_hit']:
            html += '    <p class="error">No BLAST hits found</p>\n'
            continue

        # Top hit summary
        top = result['top_hit']
        match_class = "match-high" if top['percent_identity'] >= 97 else "match-medium" if top['percent_identity'] >= 90 else "match-low"

        html += f'    <p class="{match_class}"><strong>Top Match:</strong> {top["definition"]} ({top["percent_identity"]}% identity)</p>\n'

        # All hits table
        html += """    <table>
        <tr>
            <th>Rank</th>
            <th>Species</th>
            <th>Accession</th>
            <th>Identity %</th>
            <th>E-value</th>
        </tr>
"""

        for i, hit in enumerate(result['all_hits'], 1):
            match_class = "match-high" if hit['percent_identity'] >= 97 else "match-medium" if hit['percent_identity'] >= 90 else "match-low"

            # Extract species name (first two words after first '>')
            species_name = " ".join(hit['definition'].split()[1:3])

            html += f"""        <tr class="{match_class}">
            <td>{i}</td>
            <td>{species_name}</td>
            <td>{hit['accession']}</td>
            <td>{hit['percent_identity']}%</td>
            <td>{hit['e_value']:.2e}</td>
        </tr>
"""

        html += "    </table>\n"

    html += """</body>
</html>"""

    with open(output_file, 'w') as f:
        f.write(html)

def main():
    """Main identification function"""
    if len(sys.argv) < 2:
        print("Usage: python identify_species.py <input_fasta> [output_directory]")
        print("Example: python identify_species.py results/passed_sequences.fasta results/")
        sys.exit(1)

    input_fasta = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results")
    output_dir.mkdir(exist_ok=True)

    if not input_fasta.exists():
        print(f"Input file not found: {input_fasta}")
        sys.exit(1)

    # Read sequences
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    print(f"Found {len(sequences)} sequences to identify")

    # BLAST each sequence
    results = []
    for i, record in enumerate(sequences, 1):
        print(f"\n[{i}/{len(sequences)}] Processing {record.id}")
        result = blast_sequence(str(record.seq), record.id)
        results.append(result)

        # Be nice to NCBI - wait between queries
        if i < len(sequences):
            print("  Waiting 3 seconds before next query...")
            time.sleep(3)

    # Generate outputs
    print("\n\nGenerating reports...")

    # HTML report
    html_file = output_dir / "identification_report.html"
    generate_html_report(results, html_file)
    print(f"HTML report: {html_file}")

    # JSON report
    json_file = output_dir / "identification_results.json"
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"JSON report: {json_file}")

    # Summary
    success = sum(1 for r in results if r['status'] == "SUCCESS" and r['top_hit'])
    errors = sum(1 for r in results if r['status'] == "ERROR")

    print(f"\nSummary:")
    print(f"  Total sequences: {len(results)}")
    print(f"  Successfully identified: {success}")
    print(f"  Errors: {errors}")

    if success > 0:
        print(f"\nTop matches:")
        for result in results:
            if result['status'] == "SUCCESS" and result['top_hit']:
                top = result['top_hit']
                species_name = " ".join(top['definition'].split()[1:3])
                print(f"  {result['sequence_id']}: {species_name} ({top['percent_identity']}%)")

if __name__ == "__main__":
    main()
