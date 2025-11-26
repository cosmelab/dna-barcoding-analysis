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
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.Align import PairwiseAligner
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


def trim_sequence_ends(seq_str, min_quality_window=20):
    """
    Trim N bases and low-quality regions from sequence ends

    Removes:
    - N bases from start and end
    - Ambiguous bases (non-ACGT) from ends

    This dramatically speeds up alignment and improves consensus quality

    Args:
        seq_str: DNA sequence string
        min_quality_window: minimum window size to keep (default 20bp)

    Returns:
        str: Trimmed sequence
    """
    seq = seq_str.upper()

    # Trim from start - remove N's and non-ACGT bases
    start = 0
    while start < len(seq) and seq[start] not in 'ACGT':
        start += 1

    # Trim from end - remove N's and non-ACGT bases
    end = len(seq)
    while end > start and seq[end-1] not in 'ACGT':
        end -= 1

    trimmed = seq[start:end]

    # Keep at least min_quality_window bp, otherwise return empty
    if len(trimmed) < min_quality_window:
        return ""

    return trimmed


def create_consensus(forward_seq, reverse_seq, sample_name):
    """
    Create consensus sequence from forward and reverse complement using proper alignment

    For COI barcoding (~700bp amplicon):
    - F and R reads (~650-700bp each) should overlap almost completely
    - Use pairwise alignment to find best overlap
    - Merge sequences in overlap region

    Args:
        forward_seq: SeqRecord of forward read
        reverse_seq: SeqRecord of reverse read
        sample_name: Name for the consensus sequence

    Returns:
        SeqRecord: Consensus sequence with description
    """
    # Reverse complement the reverse read
    rev_rc = reverse_seq.reverse_complement()

    # Trim low-quality ends BEFORE alignment (much faster!)
    f_seq = trim_sequence_ends(str(forward_seq.seq))
    r_seq = trim_sequence_ends(str(rev_rc.seq))

    # Check if sequences are long enough after trimming
    if len(f_seq) < 100 or len(r_seq) < 100:
        # Sequences too short after trimming - use original
        f_seq = str(forward_seq.seq).upper()
        r_seq = str(rev_rc.seq).upper()

    # For sequences with too many N's, use only central high-quality region
    # This prevents pairwise2 from hanging on low-complexity regions
    n_count_f = f_seq.count('N')
    n_count_r = r_seq.count('N')

    if n_count_f > len(f_seq) * 0.05 or n_count_r > len(r_seq) * 0.05:
        # Too many N's - use only central 400bp region for alignment
        f_mid = len(f_seq) // 2
        r_mid = len(r_seq) // 2
        f_align = f_seq[max(0, f_mid-200):f_mid+200]
        r_align = r_seq[max(0, r_mid-200):r_mid+200]
    else:
        f_align = f_seq
        r_align = r_seq

    # Use Bio.pairwise2 (deprecated but faster than PairwiseAligner for this case)
    # For COI barcoding, F and R should overlap almost completely
    # globalms: global alignment, match=2, mismatch=-1, gap_open=-2, gap_extend=-1
    # Returns list of alignments, we take the first (best-scoring) one
    alignments = pairwise2.align.globalms(f_align, r_align, 2, -1, -2, -1, one_alignment_only=True)

    if alignments:
        alignment = alignments[0]
        aligned_f = alignment.seqA
        aligned_r = alignment.seqB
        score = alignment.score

        # Calculate identity from alignment
        matches = sum(1 for a, b in zip(aligned_f, aligned_r) if a == b and a != '-')
        total = sum(1 for a, b in zip(aligned_f, aligned_r) if a != '-' or b != '-')
        identity = (matches / total * 100) if total > 0 else 0

        # Create consensus from alignment
        consensus_seq = ""
        for a, b in zip(aligned_f, aligned_r):
            if a == b:
                # Match - use the base
                consensus_seq += a if a != '-' else ''
            elif a == '-':
                # Gap in forward - use reverse
                consensus_seq += b
            elif b == '-':
                # Gap in reverse - use forward
                consensus_seq += a
            else:
                # Mismatch - use forward (arbitrary, could use quality scores)
                consensus_seq += a

        consensus_source = f"merged ({identity:.1f}% identity)"
        strategy = f"pairwise alignment with {identity:.1f}% identity, {matches} matches"

    else:
        # Fallback: if alignment fails, use longer sequence
        if len(f_seq) >= len(r_seq):
            consensus_seq = f_seq
            consensus_source = "forward"
            strategy = "alignment failed - using forward"
        else:
            consensus_seq = r_seq
            consensus_source = "reverse"
            strategy = "alignment failed - using reverse"

    # Create consensus SeqRecord
    consensus_record = SeqRecord(
        Seq(consensus_seq),
        id=sample_name,
        description=f"consensus from {forward_seq.id} and {reverse_seq.id} (strategy: {strategy})"
    )

    return consensus_record, consensus_source


def generate_html_report(output_dir, pairs, unpaired, consensus_seqs, stats):
    """Generate HTML report showing consensus generation results using new design system"""

    html_file = output_dir / "consensus_report.html"

    # Read CSS files and embed them
    import sys
    from pathlib import Path

    project_root = Path(__file__).parent.parent.parent
    base_css = (project_root / "modules/styles/base.css").read_text()
    components_css = (project_root / "modules/styles/components.css").read_text()
    reports_css = (project_root / "modules/styles/reports.css").read_text()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Consensus Sequences Report - DNA Barcoding</title>

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
        <h1>🧬 Consensus Sequences Report</h1>
        <div class="progress-badge">Step 2 of 5</div>
        <div class="report-date">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</div>
    </header>

    <!-- Summary Dashboard -->
    <section class="summary-dashboard">
        <div class="metric-card metric-success">
            <div class="metric-value">{stats['consensus_created']}</div>
            <div class="metric-label">Consensus Created (F+R)</div>
            <div class="metric-icon">✓</div>
        </div>

        <div class="metric-card metric-primary">
            <div class="metric-value">{stats['total_samples']}</div>
            <div class="metric-label">Total Samples</div>
            <div class="metric-icon">📊</div>
        </div>

        <div class="metric-card metric-{"warning" if stats['only_forward'] + stats['only_reverse'] > 0 else "success"}">
            <div class="metric-value">{stats['only_forward'] + stats['only_reverse']}</div>
            <div class="metric-label">Single Reads (F or R only)</div>
            <div class="metric-icon">{"⚠" if stats['only_forward'] + stats['only_reverse'] > 0 else "✓"}</div>
        </div>

        <div class="metric-card metric-{"warning" if stats['unpaired'] > 0 else "success"}">
            <div class="metric-value">{stats['unpaired']}</div>
            <div class="metric-label">Failed to Pair</div>
            <div class="metric-icon">{"⚠" if stats['unpaired'] > 0 else "✓"}</div>
        </div>
    </section>

    <!-- Main Content -->
    <main class="report-content">
        <div class="content-section">
            <!-- Instructions Section -->
            <div class="info-box info-tip" style="margin-bottom: 2rem;">
                <h3 style="margin-top: 0;">📖 How to Read This Report</h3>

                <h4 style="margin-top: 1rem;">Understanding Consensus Sequences</h4>
                <p>Sanger sequencing reads DNA from one direction at a time. By sequencing both <strong>forward (F)</strong> and <strong>reverse (R)</strong> strands, we can create more accurate consensus sequences.</p>

                <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin: 1rem 0;">
                    <div>
                        <h4 style="margin: 0.5rem 0;">What you'll see below:</h4>
                        <ul style="margin: 0.5rem 0;">
                            <li><strong>Forward Read:</strong> Original 5' → 3' sequence</li>
                            <li><strong>Reverse Read (Original):</strong> Raw reverse strand as sequenced</li>
                            <li><strong>Reverse Read (Rev Comp):</strong> Reverse read converted to match forward orientation</li>
                            <li><strong>Consensus Sequence:</strong> The selected sequence used for analysis (longer/better quality)</li>
                        </ul>
                    </div>
                    <div>
                        <h4 style="margin: 0.5rem 0;">Consensus Strategy:</h4>
                        <ul style="margin: 0.5rem 0;">
                            <li>Forward and reverse reads are compared</li>
                            <li>The <strong>longer</strong> sequence is selected as consensus</li>
                            <li>This provides better coverage of the target region</li>
                            <li>All sequences below are shown in <strong>full length</strong> for your review</li>
                        </ul>
                    </div>
                </div>

                <p style="margin-top: 1rem;"><strong>💡 Tip:</strong> Click each sample name below to expand and view all four sequences side-by-side. Scroll horizontally to view the complete sequences.</p>
            </div>

            <h2>Consensus Sequences</h2>
            <p class="text-secondary"><strong>👇 Click the arrows below to expand</strong> and view detailed sequence comparison and alignment</p>
"""

    # Show consensus sequences with actual sequence display
    if stats['consensus_created'] > 0:
        first_sample = True
        for sample_name in sorted(pairs.keys()):
            pair = pairs[sample_name]
            if 'F' in pair and 'R' in pair:
                consensus = next((c for c in consensus_seqs if c.id == sample_name), None)
                if consensus:
                    # Extract source from description
                    source = "forward" if "using forward" in consensus.description else "reverse"

                    # Get sequences
                    def format_seq_blocks(seq, block_size=80):
                        """Format sequence in blocks of specified size"""
                        return '\n'.join(seq[i:i+block_size] for i in range(0, len(seq), block_size))

                    f_seq_raw = str(pair['F'].seq).upper()
                    r_seq_raw = str(pair['R'].seq).upper()
                    r_rc_seq_raw = str(pair['R'].reverse_complement().seq).upper()
                    cons_seq_raw = str(consensus.seq).upper()

                    # Create alignment visualization using proper pairwise alignment
                    def create_alignment_html(seq1, seq2, seq1_name, seq2_name):
                        """Create HTML visualization using BioPython pairwise alignment"""

                        # Use pairwise2 (faster for HTML visualization)
                        alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -2, -1)

                        if not alignments:
                            return "<p>Alignment failed</p>"

                        alignment = alignments[0]
                        aligned_seq1 = alignment.seqA
                        aligned_seq2 = alignment.seqB

                        # Calculate statistics
                        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
                        mismatches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != b and a != '-' and b != '-')
                        gaps = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == '-' or b == '-')
                        total = len(aligned_seq1)
                        identity = (matches / (total - gaps) * 100) if (total - gaps) > 0 else 0

                        # Create alignment display in blocks of 80
                        alignment_html = ""
                        block_size = 80
                        for i in range(0, total, block_size):
                            block_seq1 = aligned_seq1[i:i+block_size]
                            block_seq2 = aligned_seq2[i:i+block_size]

                            # Create match line
                            match_line = ""
                            for a, b in zip(block_seq1, block_seq2):
                                if a == b and a != '-':
                                    match_line += '|'
                                elif a == '-' or b == '-':
                                    match_line += ' '
                                else:
                                    match_line += 'x'

                            alignment_html += f'''
                            <div style="margin-bottom: 1.5rem; line-height: 1.2;">
                                <div style="color: var(--cyan); font-weight: 600; margin: 0;">{seq1_name:12s} {i+1:5d}  {block_seq1}</div>
                                <div style="color: var(--text-secondary); margin: 0.1rem 0;">            {' '*5}  {match_line}</div>
                                <div style="color: var(--purple); font-weight: 600; margin: 0;">{seq2_name:12s} {i+1:5d}  {block_seq2}</div>
                            </div>'''

                        # Determine identity color and add explanation
                        identity_color = "var(--status-pass)" if identity > 70 else "var(--status-warn)" if identity > 40 else "var(--status-fail)"

                        explanation = ""
                        if identity < 70:
                            quality_msg = "⚠️ Warning" if identity > 40 else "❌ Poor Quality"
                            explanation = f'''
                            <div style="padding: 1rem; background: rgba(255, 184, 108, 0.1); border-left: 4px solid var(--status-warn); border-radius: var(--radius-md); margin-bottom: 1rem;">
                                <strong>{quality_msg}: {identity:.1f}% identity</strong>
                                <p style="margin: 0.5rem 0 0 0; font-size: 0.9rem;">Expected >95% identity for COI barcoding (~700bp amplicon with ~650-700bp reads). Low identity suggests:</p>
                                <ul style="margin: 0.5rem 0 0 1.5rem; font-size: 0.85rem;">
                                    <li>Sequencing quality issues</li>
                                    <li>Sample contamination or mixed samples</li>
                                    <li>Wrong sample pairing (F and R from different specimens)</li>
                                </ul>
                            </div>'''

                        stats_html = f'''
                        {explanation}
                        <div style="display: grid; grid-template-columns: repeat(4, 1fr); gap: 1rem; margin-bottom: 1rem; padding: 1rem; background: var(--bg-secondary); border-radius: var(--radius-md);">
                            <div style="text-align: center;">
                                <div style="font-size: 1.5rem; font-weight: 700; color: {identity_color};">{identity:.1f}%</div>
                                <div style="font-size: 0.75rem; color: var(--text-secondary); text-transform: uppercase;">Identity</div>
                            </div>
                            <div style="text-align: center;">
                                <div style="font-size: 1.5rem; font-weight: 700; color: var(--green);">{matches}</div>
                                <div style="font-size: 0.75rem; color: var(--text-secondary); text-transform: uppercase;">Matches</div>
                            </div>
                            <div style="text-align: center;">
                                <div style="font-size: 1.5rem; font-weight: 700; color: var(--status-fail);">{mismatches}</div>
                                <div style="font-size: 0.75rem; color: var(--text-secondary); text-transform: uppercase;">Mismatches</div>
                            </div>
                            <div style="text-align: center;">
                                <div style="font-size: 1.5rem; font-weight: 700; color: var(--text-secondary);">{gaps}</div>
                                <div style="font-size: 0.75rem; color: var(--text-secondary); text-transform: uppercase;">Gaps</div>
                            </div>
                        </div>'''

                        return stats_html + alignment_html

                    alignment_html = create_alignment_html(f_seq_raw, r_rc_seq_raw, "Forward", "Rev-Comp")

                    # Format sequences in blocks
                    f_seq = format_seq_blocks(f_seq_raw)
                    r_seq = format_seq_blocks(r_seq_raw)
                    r_rc_seq = format_seq_blocks(r_rc_seq_raw)
                    cons_seq = format_seq_blocks(cons_seq_raw)

                    # Create visual sequence comparison (expand first sample by default)
                    expanded_class = " expanded" if first_sample else ""
                    first_sample = False

                    html += f"""
            <div class="collapsible-section section-consensus{expanded_class}">
                <button class="collapsible-toggle" onclick="this.parentElement.classList.toggle('expanded')" style="cursor: pointer;">
                    <span class="toggle-icon">▶</span>
                    <strong>{sample_name}</strong>
                    <span class="badge badge-success">{len(consensus)} bp</span>
                    <span class="badge badge-primary">{source} used</span>
                </button>
                <div class="collapsible-content">
                    <table class="summary-table" style="margin-bottom: 1.5rem;">
                        <thead>
                            <tr>
                                <th>Forward Length</th>
                                <th>Reverse Length</th>
                                <th>Consensus Length</th>
                                <th>Source</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td>{len(pair['F'])} bp</td>
                                <td>{len(pair['R'])} bp</td>
                                <td><strong>{len(consensus)} bp</strong></td>
                                <td><span class="badge badge-primary">{source}</span></td>
                            </tr>
                        </tbody>
                    </table>

                    <h4 style="margin-top: 1rem;">Sequence Alignment</h4>
                    <p class="text-secondary" style="margin-bottom: 1rem;">Visual alignment showing overlap between Forward and Reverse Complement reads. <strong>|</strong> = match, <strong>x</strong> = mismatch</p>

                    <div class="sequence-box" style="background: var(--bg-purple-light); border-left: 3px solid var(--purple);">
                        <div class="sequence-header">Forward vs Reverse Complement Alignment</div>
                        <div class="sequence-content">{alignment_html}</div>
                    </div>

                    <h4 style="margin-top: 1.5rem;">Full Sequences</h4>

                    <div class="sequence-box" style="background: var(--bg-cyan-light); border-left: 3px solid var(--cyan);">
                        <div class="sequence-header">Forward Read ({pair['F'].id})</div>
                        <div class="sequence-content">{f_seq}</div>
                    </div>

                    <div class="sequence-box" style="background: var(--bg-pink-light); border-left: 3px solid var(--pink); margin-top: 1rem;">
                        <div class="sequence-header">Reverse Read - Original ({pair['R'].id})</div>
                        <div class="sequence-content" style="color: var(--text-secondary);">{r_seq}</div>
                    </div>

                    <div class="sequence-box" style="background: var(--bg-purple-light); border-left: 3px solid var(--purple); margin-top: 1rem;">
                        <div class="sequence-header">Reverse Read - Reverse Complement (for comparison)</div>
                        <div class="sequence-content">{r_rc_seq}</div>
                    </div>

                    <div class="sequence-box" style="background: var(--bg-green-light); border-left: 3px solid var(--green); margin-top: 1rem;">
                        <div class="sequence-header">✓ Consensus Sequence (Used for Analysis)</div>
                        <div class="sequence-content" style="color: var(--green); font-weight: 600;">{cons_seq}</div>
                        <div style="margin-top: 0.5rem; padding: 0 1rem; font-size: 0.85rem; color: var(--text-secondary);">
                            This sequence was selected as the consensus and will be used for downstream analysis
                        </div>
                    </div>
                </div>
            </div>
"""

        html += """
"""

    # Show single reads
    if stats['only_forward'] > 0 or stats['only_reverse'] > 0:
        html += """
            <h3 style="margin-top: 2rem;">Single Reads (No Pair)</h3>

            <div class="info-box info-warning">
                <strong>⚠ Single Direction Reads:</strong>
                These sequences passed QC but don't have a matching forward/reverse pair. They are included in the output as-is.
            </div>

            <table class="data-table">
                <thead>
                    <tr>
                        <th>Sample</th>
                        <th>Direction</th>
                        <th>Sequence ID</th>
                        <th>Length</th>
                    </tr>
                </thead>
                <tbody>
"""

        for sample_name in sorted(pairs.keys()):
            pair = pairs[sample_name]
            if 'F' in pair and 'R' not in pair:
                html += f"""
                    <tr class="row-warning">
                        <td><strong>{sample_name}</strong></td>
                        <td><span class="badge badge-warning">Forward only</span></td>
                        <td><code>{pair['F'].id}</code></td>
                        <td>{len(pair['F'])} bp</td>
                    </tr>
                """
            elif 'R' in pair and 'F' not in pair:
                html += f"""
                    <tr class="row-warning">
                        <td><strong>{sample_name}</strong></td>
                        <td><span class="badge badge-warning">Reverse only</span></td>
                        <td><code>{pair['R'].id}</code></td>
                        <td>{len(pair['R'])} bp</td>
                    </tr>
                """

        html += """
                </tbody>
            </table>
"""

    # Show unpaired sequences if any
    if stats['unpaired'] > 0:
        html += f"""
            <h3 style="margin-top: 2rem;">Unpaired Sequences</h3>

            <div class="info-box info-warning">
                <strong>⚠ Unpaired Sequences:</strong>
                {stats['unpaired']} sequence(s) could not be paired (missing _F or _R suffix). Included in output as-is.
            </div>
"""

    html += """
        </div>
    </main>

    <!-- Footer with Help -->
    <footer class="report-footer">
        <div class="help-section">
            <h3>Understanding Consensus Sequences</h3>

            <h4>Why create consensus sequences?</h4>
            <p>Sanger sequencing reads DNA from one direction at a time. By sequencing both forward (5' to 3') and reverse (3' to 5') strands, we can:</p>
            <ul>
                <li><strong>Improve accuracy:</strong> Combine information from both directions</li>
                <li><strong>Extend coverage:</strong> Forward and reverse reads may cover different regions</li>
                <li><strong>Validate base calls:</strong> Confirming bases from both strands reduces errors</li>
            </ul>

            <h4>How is the consensus created?</h4>
            <p>This pipeline uses a simple consensus strategy:</p>
            <ol>
                <li>Reverse complement the reverse (R) read so it matches forward (F) orientation</li>
                <li>Compare the lengths of F and R reads</li>
                <li>Use the longer read as the consensus (typically has better coverage)</li>
                <li>In cases where quality differs significantly, the higher quality read is selected</li>
            </ol>

            <h4>What about single reads?</h4>
            <p>Samples with only forward or reverse reads are still usable for:</p>
            <ul>
                <li>Phylogenetic analysis (included in tree construction)</li>
                <li>Species identification via BLAST</li>
                <li>Sequence alignment</li>
            </ul>
            <p>However, consensus sequences from F+R pairs are generally more reliable.</p>

            <h4>What if I have unpaired sequences?</h4>
            <p>Unpaired sequences lack the standard _F or _R naming convention. They are included in downstream analysis but should be reviewed to ensure proper sample naming.</p>
        </div>
    </footer>
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

    # Combine with reference sequences
    reference_file = Path("data/reference_sequences/socal_mosquitoes.fasta")
    if reference_file.exists():
        print_info("Adding reference sequences for phylogenetic analysis...")
        reference_seqs = list(SeqIO.parse(reference_file, "fasta"))
        combined_seqs = consensus_seqs + reference_seqs
        combined_fasta = args.output_dir / "combined_with_references.fasta"
        SeqIO.write(combined_seqs, combined_fasta, "fasta")
        print_success(f"Combined sequences: {combined_fasta} ({len(consensus_seqs)} samples + {len(reference_seqs)} references)")
    else:
        print_info(f"⚠ Reference file not found: {reference_file}")

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
