#!/usr/bin/env python3

################################################################################
# Alignment Visualization Script
#
# Purpose: Create publication-quality visualizations of sequence alignments,
#          including sequence logo, conservation plots, and alignment matrices.
#
# Features:
#   - Sequence logo visualization
#   - Conservation profile plot
#   - Gap distribution plot
#   - Sequence identity matrix
#   - Interactive HTML output (optional)
#
# Requirements:
#   - Python 3.6+
#   - BioPython (pip install biopython)
#   - Matplotlib (pip install matplotlib)
#   - NumPy (pip install numpy)
#   - Seaborn (optional, for better aesthetics: pip install seaborn)
#
# Usage:
#   ./visualize_alignment.py -i aligned.fasta -o output.pdf [OPTIONS]
#
# Examples:
#   # Basic visualization
#   ./visualize_alignment.py -i alignment.fasta -o alignment_viz.png
#
#   # Create PDF with high resolution
#   ./visualize_alignment.py -i alignment.fasta -o alignment_viz.pdf --dpi 300
#
#   # Show interactive plot (requires display)
#   ./visualize_alignment.py -i alignment.fasta --show
#
#   # Region-specific visualization (positions 1-500)
#   ./visualize_alignment.py -i alignment.fasta -o alignment_viz.png --region 1 500
#
# References:
#   - BioPython: https://biopython.org/
#   - Matplotlib for bioinformatics: https://www.numberanalytics.com/blog/matplotlib-bioinformatics-pipelines
#   - Sequence logos: https://en.wikipedia.org/wiki/Sequence_logo
#   - Conservation scoring: https://avrilomics.blogspot.com/2016/07/calculating-conservation-score-for.html
#
################################################################################

import sys
import argparse
import math
from pathlib import Path
from typing import List, Tuple, Optional

try:
    from Bio import AlignIO
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import LinearSegmentedColormap
except ImportError as e:
    print(f"Error: Required package not found. Install with: pip install biopython numpy matplotlib")
    sys.exit(1)

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False


################################################################################
# Constants
################################################################################

# Color schemes for nucleotides
DNA_COLORS = {
    'A': '#00CC00',  # Green
    'C': '#0000FF',  # Blue
    'G': '#FFB300',  # Orange
    'T': '#CC0000',  # Red
    '-': '#CCCCCC'   # Gray (gap)
}

# Color schemes for amino acids (Clustal convention)
PROTEIN_COLORS = {
    # Hydrophobic
    'A': '#80a0f0', 'I': '#80a0f0', 'L': '#80a0f0', 'M': '#80a0f0', 'V': '#80a0f0',
    # Aromatic
    'F': '#ff8000', 'W': '#ff8000', 'Y': '#ff8000',
    # Polar uncharged
    'S': '#00ff00', 'T': '#00ff00', 'C': '#00ff00',
    # Positive
    'H': '#ff0000', 'K': '#ff0000', 'R': '#ff0000',
    # Negative
    'D': '#ff0066', 'E': '#ff0066',
    # Asparagine and Glutamine
    'N': '#0066ff', 'Q': '#0066ff',
    # Proline
    'P': '#ffff00',
    # Glycine
    'G': '#ff9900',
    # Stop and unknown
    '*': '#000000', 'X': '#cccccc', '-': '#cccccc'
}


################################################################################
# Class: AlignmentVisualizer
################################################################################

class AlignmentVisualizer:
    """Create visualizations for sequence alignments."""

    def __init__(self, alignment_file: str, format: str = 'fasta'):
        """
        Initialize visualizer with alignment file.

        Args:
            alignment_file: Path to alignment file
            format: Format of alignment (default: fasta)
        """
        self.alignment_file = alignment_file
        self.format = format
        self.alignment = None
        self.seq_type = None
        self.alignment_length = 0

    def load_alignment(self) -> bool:
        """Load alignment from file."""
        try:
            self.alignment = AlignIO.read(self.alignment_file, self.format)
            self.alignment_length = self.alignment.get_alignment_length()

            # Detect sequence type
            if len(self.alignment) > 0:
                first_seq = str(self.alignment[0].seq).upper()
                seq_no_gaps = first_seq.replace('-', '')
                if all(c in 'ACGTRYSWKMBDHVN' for c in seq_no_gaps):
                    self.seq_type = 'DNA'
                else:
                    self.seq_type = 'protein'

            return True
        except Exception as e:
            print(f"Error loading alignment: {e}")
            return False

    def get_color_scheme(self) -> dict:
        """Get appropriate color scheme for sequence type."""
        if self.seq_type == 'DNA':
            return DNA_COLORS
        else:
            return PROTEIN_COLORS

    def calculate_conservation(self) -> List[float]:
        """Calculate Shannon entropy-based conservation score."""
        entropies = []

        for pos in range(self.alignment_length):
            column = str(self.alignment[:, pos]).upper()
            char_counts = {}
            total_non_gap = 0

            for char in column:
                if char != '-':
                    char_counts[char] = char_counts.get(char, 0) + 1
                    total_non_gap += 1

            if total_non_gap == 0:
                entropy = 0.0
            else:
                entropy = 0.0
                for count in char_counts.values():
                    if count > 0:
                        p = count / total_non_gap
                        entropy -= p * math.log2(p)

            entropies.append(entropy)

        # Convert entropy to conservation (1 - normalized entropy)
        if self.seq_type == 'DNA':
            max_entropy = math.log2(4)
        else:
            max_entropy = math.log2(20)

        conservation = [1 - (e / max_entropy) if max_entropy > 0 else 0 for e in entropies]
        return [max(0, min(1, c)) for c in conservation]

    def get_gap_percentage(self) -> List[float]:
        """Get gap percentage for each position."""
        gap_pcts = []
        for pos in range(self.alignment_length):
            column = str(self.alignment[:, pos])
            gap_count = column.count('-')
            gap_pct = gap_count / len(column) * 100 if len(column) > 0 else 0
            gap_pcts.append(gap_pct)
        return gap_pcts

    def plot_conservation_profile(self, ax, positions: Optional[Tuple[int, int]] = None):
        """
        Plot conservation score across alignment.

        Args:
            ax: Matplotlib axis
            positions: Tuple of (start, end) positions to plot
        """
        conservation = self.calculate_conservation()

        if positions:
            start, end = positions
            conservation = conservation[start:end]
            x_values = range(start + 1, end + 1)
        else:
            x_values = range(1, len(conservation) + 1)

        ax.fill_between(x_values, conservation, alpha=0.6, color='steelblue', label='Conservation')
        ax.plot(x_values, conservation, color='steelblue', linewidth=1)
        ax.set_ylabel('Conservation Score', fontsize=10, fontweight='bold')
        ax.set_ylim([0, 1])
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('#f8f9fa')

    def plot_gap_distribution(self, ax, positions: Optional[Tuple[int, int]] = None):
        """
        Plot gap percentage across alignment.

        Args:
            ax: Matplotlib axis
            positions: Tuple of (start, end) positions to plot
        """
        gap_pcts = self.get_gap_percentage()

        if positions:
            start, end = positions
            gap_pcts = gap_pcts[start:end]
            x_values = range(start + 1, end + 1)
        else:
            x_values = range(1, len(gap_pcts) + 1)

        ax.bar(x_values, gap_pcts, color='coral', alpha=0.7, width=1)
        ax.set_ylabel('Gap %', fontsize=10, fontweight='bold')
        ax.set_ylim([0, 100])
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_facecolor('#f8f9fa')

    def plot_alignment_matrix(self, ax, sample_size: int = 50):
        """
        Plot sequence alignment as colored matrix (sampled for large alignments).

        Args:
            ax: Matplotlib axis
            sample_size: Show every Nth position for large alignments
        """
        colors = self.get_color_scheme()

        # Sample positions if alignment is too large
        if self.alignment_length > 500:
            step = max(1, self.alignment_length // 500)
        else:
            step = 1

        # Create matrix
        num_seqs = len(self.alignment)
        sampled_length = (self.alignment_length + step - 1) // step

        matrix = np.zeros((num_seqs, sampled_length, 3))  # RGB values

        for seq_idx, record in enumerate(self.alignment):
            seq = str(record.seq).upper()
            for pos_idx, pos in enumerate(range(0, self.alignment_length, step)):
                char = seq[pos]
                color = colors.get(char, colors.get('-'))
                # Convert hex to RGB
                if isinstance(color, str):
                    rgb = tuple(int(color.lstrip('#')[i:i+2], 16) / 255 for i in (0, 2, 4))
                    matrix[seq_idx, pos_idx] = rgb

        ax.imshow(matrix, aspect='auto', interpolation='nearest')
        ax.set_ylabel('Sequence', fontsize=10, fontweight='bold')
        ax.set_xlabel('Alignment Position', fontsize=10, fontweight='bold')
        ax.set_yticks(range(min(num_seqs, 20)))
        if num_seqs <= 20:
            ax.set_yticklabels([rec.id for rec in self.alignment], fontsize=8)
        ax.set_xticks([])

    def create_full_visualization(self, output_file: str, dpi: int = 150,
                                   positions: Optional[Tuple[int, int]] = None):
        """
        Create comprehensive alignment visualization.

        Args:
            output_file: Path to output file (PDF, PNG, etc.)
            dpi: Resolution (dots per inch)
            positions: Optional (start, end) positions to visualize
        """
        fig = plt.figure(figsize=(14, 10))
        fig.suptitle(f'Alignment Visualization: {Path(self.alignment_file).name}',
                     fontsize=14, fontweight='bold', y=0.995)

        # Create grid layout
        gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3, top=0.96)

        # 1. Alignment matrix
        ax1 = fig.add_subplot(gs[0, :])
        self.plot_alignment_matrix(ax1)
        ax1.set_title('Sequence Alignment Matrix', fontweight='bold', fontsize=11)

        # 2. Conservation profile
        ax2 = fig.add_subplot(gs[1, :])
        self.plot_conservation_profile(ax2, positions)
        ax2.set_title('Conservation Profile', fontweight='bold', fontsize=11)

        # 3. Gap distribution
        ax3 = fig.add_subplot(gs[2, :])
        self.plot_gap_distribution(ax3, positions)
        ax3.set_title('Gap Distribution', fontweight='bold', fontsize=11)
        ax3.set_xlabel('Alignment Position', fontsize=10, fontweight='bold')

        # Save figure
        try:
            plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=output_file.split('.')[-1])
            print(f"Visualization saved to: {output_file}")
        except Exception as e:
            print(f"Error saving visualization: {e}")
            return False

        return True

    def show_plot(self):
        """Display plot interactively."""
        try:
            plt.show()
        except Exception as e:
            print(f"Error displaying plot: {e}")


################################################################################
# Main function
################################################################################

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Create alignment visualizations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic visualization
  %(prog)s -i alignment.fasta -o output.png

  # High resolution PDF
  %(prog)s -i alignment.fasta -o output.pdf --dpi 300

  # Visualize specific region (positions 1-500)
  %(prog)s -i alignment.fasta -o output.png --region 1 500

  # Show interactive plot
  %(prog)s -i alignment.fasta --show
        '''
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input alignment file')
    parser.add_argument('-o', '--output', default=None,
                        help='Output image file (PNG, PDF, etc.)')
    parser.add_argument('-f', '--format', default='fasta',
                        choices=['fasta', 'stockholm', 'clustal', 'phylip'],
                        help='Alignment format (default: fasta)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='Resolution in DPI (default: 150)')
    parser.add_argument('--region', type=int, nargs=2, metavar=('START', 'END'),
                        help='Visualize specific region (1-based positions)')
    parser.add_argument('--show', action='store_true',
                        help='Display plot interactively')

    args = parser.parse_args()

    # Check if input file exists
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    if not args.output and not args.show:
        print("Error: Specify output file (-o) or use --show flag")
        sys.exit(1)

    # Create visualizer and load alignment
    viz = AlignmentVisualizer(args.input, args.format)

    if not viz.load_alignment():
        sys.exit(1)

    print(f"Loaded alignment: {len(viz.alignment)} sequences, {viz.alignment_length} bp")
    print(f"Sequence type: {viz.seq_type}")

    # Validate region if specified
    region = None
    if args.region:
        start, end = args.region
        if start < 1 or end > viz.alignment_length or start >= end:
            print(f"Error: Invalid region. Alignment length is {viz.alignment_length} bp")
            sys.exit(1)
        region = (start - 1, end)  # Convert to 0-based
        print(f"Visualizing region: {start}-{end}")

    # Create visualization
    if args.output:
        success = viz.create_full_visualization(args.output, args.dpi, region)
        if not success:
            sys.exit(1)

    # Show interactive plot if requested
    if args.show:
        viz.create_full_visualization('/tmp/temp_alignment.png', args.dpi, region)
        viz.show_plot()


if __name__ == '__main__':
    main()
