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

def extract_species_name(definition):
    """Extract species name from BLAST definition"""
    # Definition format: "gi|xxx|gb|YYY| Genus species ..."
    parts = definition.split()
    # Find genus and species (usually after accession info)
    for i, part in enumerate(parts):
        if not part.startswith(('gi|', 'gb|', 'ref|')) and len(part) > 2:
            # Found genus, get genus + species
            if i + 1 < len(parts):
                return f"{parts[i]} {parts[i+1]}"
            return parts[i]
    return " ".join(parts[1:3])  # Fallback

def get_common_name(species):
    """Get common name for known species"""
    common_names = {
        "Aedes albopictus": "Asian tiger mosquito",
        "Aedes aegypti": "Yellow fever mosquito",
        "Culex pipiens": "Northern house mosquito",
        "Culex quinquefasciatus": "Southern house mosquito",
        "Culex tarsalis": "Western encephalitis mosquito",
        "Anopheles gambiae": "African malaria mosquito",
        "Anopheles stephensi": "Asian malaria mosquito",
        "Culex erraticus": "Swamp mosquito",
    }
    for sci_name, common in common_names.items():
        if sci_name.lower() in species.lower():
            return common
    return ""

def generate_html_report(results, output_file):
    """Generate HTML report with BLAST results using new design system"""
    from datetime import datetime
    from pathlib import Path

    # Calculate summary statistics
    total = len(results)
    success = sum(1 for r in results if r['status'] == "SUCCESS" and r['top_hit'])
    errors = sum(1 for r in results if r['status'] == "ERROR")

    # Get unique species
    species_found = set()
    for result in results:
        if result['status'] == "SUCCESS" and result['top_hit']:
            species = extract_species_name(result['top_hit']['definition'])
            species_found.add(species)

    # Average identity
    identities = [r['top_hit']['percent_identity'] for r in results
                  if r['status'] == "SUCCESS" and r['top_hit']]
    avg_identity = sum(identities) / len(identities) if identities else 0

    # Read CSS files and embed them
    project_root = Path(__file__).parent.parent.parent
    base_css = (project_root / "modules/styles/base.css").read_text()
    components_css = (project_root / "modules/styles/components.css").read_text()
    reports_css = (project_root / "modules/styles/reports.css").read_text()

    # Build HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Species Identification Report - DNA Barcoding</title>

    <!-- Embedded CSS for reliable loading -->
    <style>
{base_css}

{components_css}

{reports_css}
    </style>
</head>
<body>
    <!-- Report Header -->
    <header class="report-header">
        <h1>🔬 Species Identification Report</h1>
        <div class="progress-badge">Step 5 of 5</div>
        <div class="report-date">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</div>
    </header>

    <!-- Summary Dashboard -->
    <section class="summary-dashboard">
        <div class="metric-card metric-success">
            <div class="metric-value">{success}</div>
            <div class="metric-label">Successfully Identified</div>
            <div class="metric-icon">✓</div>
        </div>

        <div class="metric-card metric-primary">
            <div class="metric-value">{len(species_found)}</div>
            <div class="metric-label">Unique Species</div>
            <div class="metric-icon">🦟</div>
        </div>

        <div class="metric-card metric-primary">
            <div class="metric-value">{avg_identity:.1f}%</div>
            <div class="metric-label">Avg. Identity</div>
            <div class="metric-icon">🎯</div>
        </div>

        <div class="metric-card metric-{"fail" if errors > 0 else "success"}">
            <div class="metric-value">{errors}</div>
            <div class="metric-label">Errors</div>
            <div class="metric-icon">{"⚠" if errors > 0 else "✓"}</div>
        </div>
    </section>

    <!-- Main Content -->
    <main class="report-content">
        <div class="content-section">
            <h2>Identification Results</h2>
            <p>BLAST search against NCBI GenBank nucleotide database (nt)</p>

            <div class="info-box info-tip">
                <strong>💡 How to interpret results:</strong>
                <ul>
                    <li><strong>≥97% identity:</strong> Same species (high confidence)</li>
                    <li><strong>90-97% identity:</strong> Related species (moderate confidence)</li>
                    <li><strong>&lt;90% identity:</strong> Distant match (low confidence)</li>
                </ul>
            </div>
"""

    # Generate results for each sequence
    for result in results:
        seq_id = result['sequence_id']

        if result['status'] == "ERROR":
            html += f"""
            <div class="collapsible-section">
                <button class="collapsible-toggle">
                    <span class="toggle-icon">▶</span>
                    <strong>{seq_id}</strong>
                    <span class="badge badge-fail">ERROR</span>
                </button>
                <div class="collapsible-content">
                    <div class="info-box info-error">
                        <strong>⚠️ BLAST Error:</strong> {result['error']}
                    </div>
                </div>
            </div>
"""
            continue

        if not result['top_hit']:
            html += f"""
            <div class="collapsible-section">
                <button class="collapsible-toggle">
                    <span class="toggle-icon">▶</span>
                    <strong>{seq_id}</strong>
                    <span class="badge badge-warning">NO HITS</span>
                </button>
                <div class="collapsible-content">
                    <div class="info-box info-warning">
                        <strong>⚠️ No BLAST hits found</strong>
                        This sequence may be too short or of poor quality.
                    </div>
                </div>
            </div>
"""
            continue

        # Successful identification
        top = result['top_hit']
        species = extract_species_name(top['definition'])
        common = get_common_name(species)
        identity = top['percent_identity']

        # Determine badge color
        badge_class = "badge-success" if identity >= 97 else "badge-warning" if identity >= 90 else "badge-fail"

        html += f"""
            <div class="collapsible-section expanded">
                <button class="collapsible-toggle">
                    <span class="toggle-icon">▶</span>
                    <strong>{seq_id}</strong>
                    <span class="badge {badge_class}">{identity}% match</span>
                    <em style="margin-left: auto; color: var(--text-secondary);">{species}</em>
                </button>
                <div class="collapsible-content">
                    <!-- Top Match Summary -->
                    <div class="info-box info-success">
                        <strong>🎯 Top Match:</strong>
                        <div style="margin-top: 0.5rem;">
                            <div style="font-size: 1.1rem; font-weight: 600; font-style: italic;">
                                {species}
                            </div>
                            {f'<div style="color: var(--text-secondary); margin-top: 0.25rem;">{common}</div>' if common else ''}
                            <div style="margin-top: 0.5rem;">
                                <span class="badge {badge_class}">{identity}% identity</span>
                                <span class="badge badge-secondary">E-value: {top['e_value']:.2e}</span>
                                <span class="badge badge-secondary">Accession: {top['accession']}</span>
                            </div>
                        </div>
                    </div>

                    <!-- All BLAST Hits Table -->
                    <h4 style="margin-top: 1.5rem;">All BLAST Hits (Top 5)</h4>
                    <table class="data-table">
                        <thead>
                            <tr>
                                <th>Rank</th>
                                <th>Species</th>
                                <th>Accession</th>
                                <th>Identity %</th>
                                <th>E-value</th>
                            </tr>
                        </thead>
                        <tbody>
"""

        # Add all hits
        for i, hit in enumerate(result['all_hits'], 1):
            hit_species = extract_species_name(hit['definition'])
            hit_identity = hit['percent_identity']
            row_class = "row-pass" if hit_identity >= 97 else "row-warning" if hit_identity >= 90 else "row-fail"
            badge_class = "badge-success" if hit_identity >= 97 else "badge-warning" if hit_identity >= 90 else "badge-fail"

            html += f"""
                            <tr class="{row_class}">
                                <td>{i}</td>
                                <td><em>{hit_species}</em></td>
                                <td><code>{hit['accession']}</code></td>
                                <td><span class="badge {badge_class}">{hit_identity}%</span></td>
                                <td>{hit['e_value']:.2e}</td>
                            </tr>
"""

        html += """
                        </tbody>
                    </table>
                </div>
            </div>
"""

    html += """
        </div>
    </main>

    <!-- Footer with Help -->
    <footer class="report-footer">
        <div class="help-section">
            <h3>Understanding Your Results</h3>

            <h4>What is BLAST?</h4>
            <p>BLAST (Basic Local Alignment Search Tool) compares your DNA sequence against millions of sequences in GenBank to find the closest matches.</p>

            <h4>What does % identity mean?</h4>
            <ul>
                <li><strong>99-100%:</strong> Almost certainly the same species</li>
                <li><strong>97-99%:</strong> Same species (standard DNA barcoding threshold)</li>
                <li><strong>95-97%:</strong> Likely same species, possibly subspecies variation</li>
                <li><strong>90-95%:</strong> Closely related species</li>
                <li><strong>&lt;90%:</strong> Different species, possibly same genus</li>
            </ul>

            <h4>Why do multiple hits have identical % identity?</h4>
            <div class="info-box info-success">
                <strong>✓ This is completely normal!</strong>
                <p style="margin-top: 0.5rem;">You might see several BLAST hits with the exact same % identity (e.g., 99.24%). This happens because:</p>
                <ul style="margin-top: 0.5rem; margin-bottom: 0;">
                    <li><strong>Same species, different database entries:</strong> GenBank contains multiple entries for the same species from different studies, locations, or researchers</li>
                    <li><strong>COI is conserved:</strong> The COI barcode region is highly conserved within species, so individuals have nearly identical sequences</li>
                    <li><strong>Different voucher specimens:</strong> Each entry represents a different specimen (e.g., RMNH.INS.1271356, Cx_pip_2) but from the same species</li>
                    <li><strong>Multiple accessions:</strong> The same organism's genome may be submitted multiple times with different accession numbers</li>
                </ul>
                <p style="margin-top: 0.5rem; margin-bottom: 0;"><strong>What to do:</strong> Look at the species names - if they're the same across the top hits, you have a confident identification! The multiple entries actually increase your confidence that the match is correct.</p>
            </div>

            <h4>What if my top match isn't a mosquito?</h4>
            <p>This could indicate:</p>
            <ul>
                <li>Contamination in your sample</li>
                <li>Poor sequence quality leading to incorrect matches</li>
                <li>A novel or undescribed species</li>
            </ul>
            <p>Review your chromatograms and quality scores if results seem unexpected.</p>
        </div>
    </footer>

    <!-- JavaScript for collapsible sections -->
    <script>
        document.querySelectorAll('.collapsible-toggle').forEach(button => {
            button.addEventListener('click', () => {
                const section = button.parentElement;
                section.classList.toggle('expanded');
            });
        });
    </script>
</body>
</html>
"""

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
