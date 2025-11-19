#!/usr/bin/env python3

################################################################################
# Alignment Statistics Calculator
#
# Purpose: Calculate alignment quality metrics including gap content,
#          conservation scores, and sequence statistics.
#
# Features:
#   - Gap percentage per position and per sequence
#   - Conservation scores (identity and similarity)
#   - Shannon entropy for measuring sequence variation
#   - Alignment quality summary
#   - Position-by-position statistics
#
# Requirements:
#   - Python 3.6+
#   - BioPython (pip install biopython)
#   - NumPy (pip install numpy)
#
# Usage:
#   ./alignment_stats.py -i aligned.fasta [-o output.csv] [-f fasta|stockholm]
#
# Examples:
#   # Basic statistics on command line
#   ./alignment_stats.py -i alignment.fasta
#
#   # Save detailed statistics to CSV
#   ./alignment_stats.py -i alignment.fasta -o stats.csv
#
#   # Process Stockholm format alignment
#   ./alignment_stats.py -i alignment.sto -f stockholm
#
# References:
#   - BioPython AlignIO: https://biopython.org/wiki/AlignIO
#   - Shannon Entropy for conservation:
#     https://avrilomics.blogspot.com/2016/07/calculating-conservation-score-for.html
#   - Multiple Sequence Alignment Scoring:
#     https://www.proteinstructures.com/multiple-sequence-alignment/
#
################################################################################

import sys
import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

try:
    from Bio import AlignIO
    from Bio.Seq import Seq
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("Error: NumPy is required. Install with: pip install numpy")
    sys.exit(1)


################################################################################
# Constants
################################################################################

# Standard genetic code (DNA nucleotides)
IUPAC_DNA = set('ACGT')
IUPAC_DNA_AMBIG = set('ACGTRYSWKMBDHVN-')

# Standard amino acids
IUPAC_PROTEIN = set('ACDEFGHIKLMNPQRSTVWY')
IUPAC_PROTEIN_AMBIG = set('ACDEFGHIKLMNPQRSTVWYXBZJ*-')

# Color scheme for amino acids (Clustal)
# Used for similarity comparison
CLUSTAL_GROUPS = {
    'AVLIM': 'strong hydrophobic',
    'FWY': 'aromatic',
    'ST': 'polar uncharged',
    'RHK': 'positive charge',
    'DENQ': 'negative charge and asparagine',
    'C': 'cysteine'
}


################################################################################
# Class: AlignmentStats
################################################################################

class AlignmentStats:
    """Calculate statistics for sequence alignments."""

    def __init__(self, alignment_file: str, format: str = 'fasta'):
        """
        Initialize with alignment file.

        Args:
            alignment_file: Path to alignment file
            format: Format of alignment (default: fasta)
        """
        self.alignment_file = alignment_file
        self.format = format
        self.alignment = None
        self.num_sequences = 0
        self.alignment_length = 0
        self.seq_type = None  # 'DNA' or 'protein'

    def load_alignment(self) -> bool:
        """
        Load alignment from file.

        Returns:
            True if successful, False otherwise
        """
        try:
            self.alignment = AlignIO.read(self.alignment_file, self.format)
            self.num_sequences = len(self.alignment)
            self.alignment_length = self.alignment.get_alignment_length()

            # Detect sequence type
            if self.num_sequences > 0:
                first_seq = str(self.alignment[0].seq).upper()
                # Remove gaps for type detection
                seq_no_gaps = first_seq.replace('-', '')
                if all(c in IUPAC_DNA_AMBIG for c in seq_no_gaps):
                    self.seq_type = 'DNA'
                else:
                    self.seq_type = 'protein'

            return True
        except Exception as e:
            print(f"Error loading alignment: {e}")
            return False

    def get_gap_percentage(self) -> float:
        """Calculate overall gap percentage in alignment."""
        total_chars = self.num_sequences * self.alignment_length
        gap_count = sum(str(record.seq).count('-') for record in self.alignment)
        return (gap_count / total_chars * 100) if total_chars > 0 else 0.0

    def get_position_gap_percentage(self) -> List[float]:
        """Calculate gap percentage for each position."""
        gap_percentages = []
        for pos in range(self.alignment_length):
            column = self.alignment[:, pos]
            gap_count = column.count('-')
            gap_pct = (gap_count / len(column) * 100) if len(column) > 0 else 0.0
            gap_percentages.append(gap_pct)
        return gap_percentages

    def get_per_sequence_stats(self) -> Dict[str, Dict]:
        """
        Calculate statistics for each sequence.

        Returns:
            Dictionary with per-sequence statistics
        """
        stats = {}
        for record in self.alignment:
            seq_str = str(record.seq)
            gap_count = seq_str.count('-')
            total_len = len(seq_str)

            stats[record.id] = {
                'length': total_len,
                'gaps': gap_count,
                'gap_percentage': (gap_count / total_len * 100) if total_len > 0 else 0.0,
                'non_gap_length': total_len - gap_count
            }
        return stats

    def calculate_shannon_entropy(self) -> List[float]:
        """
        Calculate Shannon entropy for each position.

        Shannon entropy measures sequence variation at each position:
        H = -sum(p_i * log2(p_i)) for each amino acid/nucleotide i
        High entropy = high variation, Low entropy = high conservation

        Returns:
            List of entropy values for each position
        """
        entropies = []

        for pos in range(self.alignment_length):
            column = str(self.alignment[:, pos]).upper()

            # Count characters excluding gaps
            char_counts = {}
            total_non_gap = 0

            for char in column:
                if char != '-':
                    char_counts[char] = char_counts.get(char, 0) + 1
                    total_non_gap += 1

            # Calculate entropy
            if total_non_gap == 0:
                entropy = 0.0
            else:
                entropy = 0.0
                for count in char_counts.values():
                    if count > 0:
                        p = count / total_non_gap
                        entropy -= p * math.log2(p)

            entropies.append(entropy)

        return entropies

    def calculate_conservation_score(self) -> List[float]:
        """
        Calculate conservation score for each position.

        Conservation = 1 - (entropy / max_entropy)
        Returns value between 0 (highly variable) and 1 (fully conserved)

        Returns:
            List of conservation scores for each position
        """
        entropies = self.calculate_shannon_entropy()

        # Maximum entropy depends on sequence type and alignment
        # For DNA: log2(4) = 2.0, for protein: log2(20) = 4.32
        if self.seq_type == 'DNA':
            max_entropy = math.log2(4)  # 4 nucleotides
        else:
            max_entropy = math.log2(20)  # 20 amino acids

        conservation_scores = []
        for entropy in entropies:
            if max_entropy > 0:
                score = 1 - (entropy / max_entropy)
                # Clamp to [0, 1]
                score = max(0, min(1, score))
            else:
                score = 0.0

            conservation_scores.append(score)

        return conservation_scores

    def calculate_identity(self) -> float:
        """
        Calculate overall sequence identity in alignment.

        Identity = percentage of positions where all sequences are identical

        Returns:
            Percentage of identical positions
        """
        identical_positions = 0

        for pos in range(self.alignment_length):
            column = str(self.alignment[:, pos]).upper()
            # Remove gaps
            column_no_gap = column.replace('-', '')

            if len(column_no_gap) > 0:
                # Check if all non-gap characters are the same
                if len(set(column_no_gap)) == 1:
                    identical_positions += 1

        return (identical_positions / self.alignment_length * 100) if self.alignment_length > 0 else 0.0

    def get_summary_statistics(self) -> Dict:
        """Get comprehensive alignment statistics."""
        return {
            'file': self.alignment_file,
            'format': self.format,
            'seq_type': self.seq_type,
            'num_sequences': self.num_sequences,
            'alignment_length': self.alignment_length,
            'total_characters': self.num_sequences * self.alignment_length,
            'gap_percentage': round(self.get_gap_percentage(), 2),
            'overall_identity': round(self.calculate_identity(), 2),
            'mean_entropy': round(np.mean(self.calculate_shannon_entropy()), 3),
            'mean_conservation': round(np.mean(self.calculate_conservation_score()), 3)
        }

    def print_summary(self):
        """Print alignment summary statistics."""
        stats = self.get_summary_statistics()

        print("\n" + "=" * 70)
        print("ALIGNMENT STATISTICS SUMMARY")
        print("=" * 70)
        print(f"File: {stats['file']}")
        print(f"Format: {stats['format']}")
        print(f"Sequence Type: {stats['seq_type']}")
        print(f"\nAlignment Dimensions:")
        print(f"  Number of sequences: {stats['num_sequences']}")
        print(f"  Alignment length: {stats['alignment_length']} bp")
        print(f"  Total characters: {stats['total_characters']:,}")
        print(f"\nAlignment Quality Metrics:")
        print(f"  Overall gap percentage: {stats['gap_percentage']}%")
        print(f"  Overall sequence identity: {stats['overall_identity']}%")
        print(f"  Mean Shannon entropy: {stats['mean_entropy']} bits")
        print(f"  Mean conservation score: {stats['mean_conservation']} (0=variable, 1=conserved)")
        print("=" * 70 + "\n")

    def print_per_sequence_stats(self):
        """Print statistics for each sequence."""
        stats = self.get_per_sequence_stats()

        print("\n" + "=" * 80)
        print("PER-SEQUENCE STATISTICS")
        print("=" * 80)
        print(f"{'Sequence ID':<30} {'Length':>8} {'Gaps':>8} {'Gap %':>8} {'Non-gap':>8}")
        print("-" * 80)

        for seq_id, seq_stats in stats.items():
            print(f"{seq_id:<30} {seq_stats['length']:>8} {seq_stats['gaps']:>8} "
                  f"{seq_stats['gap_percentage']:>7.1f}% {seq_stats['non_gap_length']:>8}")

        print("=" * 80 + "\n")

    def export_position_statistics_csv(self, output_file: str):
        """
        Export position-by-position statistics to CSV file.

        Args:
            output_file: Path to output CSV file
        """
        gap_pcts = self.get_position_gap_percentage()
        entropies = self.calculate_shannon_entropy()
        conservation = self.calculate_conservation_score()

        try:
            with open(output_file, 'w') as f:
                f.write("Position,GapPercentage,ShannonEntropy,ConservationScore\n")

                for pos in range(self.alignment_length):
                    f.write(f"{pos + 1},{gap_pcts[pos]:.2f},{entropies[pos]:.3f},{conservation[pos]:.3f}\n")

            print(f"Position statistics exported to: {output_file}")
        except Exception as e:
            print(f"Error exporting statistics: {e}")
            sys.exit(1)


################################################################################
# Main function
################################################################################

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Calculate alignment quality statistics',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic statistics
  %(prog)s -i alignment.fasta

  # Export position statistics to CSV
  %(prog)s -i alignment.fasta -o stats.csv

  # Process Stockholm format
  %(prog)s -i alignment.sto -f stockholm
        '''
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input alignment file')
    parser.add_argument('-o', '--output', default=None,
                        help='Output CSV file for position statistics (optional)')
    parser.add_argument('-f', '--format', default='fasta',
                        choices=['fasta', 'stockholm', 'clustal', 'phylip'],
                        help='Alignment format (default: fasta)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print per-sequence statistics')

    args = parser.parse_args()

    # Check if input file exists
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    # Create stats calculator and load alignment
    stats = AlignmentStats(args.input, args.format)

    if not stats.load_alignment():
        sys.exit(1)

    # Print summary statistics
    stats.print_summary()

    # Print per-sequence statistics if verbose
    if args.verbose:
        stats.print_per_sequence_stats()

    # Export position statistics if requested
    if args.output:
        stats.export_position_statistics_csv(args.output)


if __name__ == '__main__':
    main()
