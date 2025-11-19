#!/usr/bin/env python3

################################################################################
# Alignment Trimming and Cleaning Script
#
# Purpose: Remove poorly aligned and ambiguous regions from sequence alignments.
#          Implements trimming strategies for improving phylogenetic analysis.
#
# Features:
#   - Gap-based trimming (remove positions with high gap content)
#   - Similarity-based trimming (remove low conservation regions)
#   - Sequence-based trimming (remove sequences with excessive gaps)
#   - Entropy-based trimming (remove high-entropy regions)
#   - Sliding window approach for contiguous block selection
#   - Multiple trimming strategies (strict/moderate/lenient)
#
# Requirements:
#   - Python 3.6+
#   - BioPython (pip install biopython)
#   - NumPy (pip install numpy)
#
# Usage:
#   ./trim_alignment.py -i aligned.fasta -o trimmed.fasta [OPTIONS]
#
# Examples:
#   # Remove positions with >50% gaps
#   ./trim_alignment.py -i alignment.fasta -o trimmed.fasta --max-gap 50
#
#   # Strict trimming (>10% gaps allowed)
#   ./trim_alignment.py -i alignment.fasta -o trimmed.fasta --strict
#
#   # Remove sequences with >30% gaps
#   ./trim_alignment.py -i alignment.fasta -o trimmed.fasta --seq-gap 30
#
#   # Entropy-based trimming (remove variable regions)
#   ./trim_alignment.py -i alignment.fasta -o trimmed.fasta --min-conservation 0.5
#
# References:
#   - trimAl: https://pmc.ncbi.nlm.nih.gov/articles/PMC2712344/
#   - Gblocks comparison: https://pmc.ncbi.nlm.nih.gov/articles/PMC4538881/
#   - pytrimal: https://github.com/althonos/pytrimal
#   - Gap penalty considerations: https://www.biostars.org/p/6522/
#
################################################################################

import sys
import argparse
import math
from pathlib import Path
from typing import List, Tuple, Set

try:
    from Bio import AlignIO, SeqIO
    from Bio.SeqRecord import SeqRecord
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
# Class: AlignmentTrimmer
################################################################################

class AlignmentTrimmer:
    """Trim and clean sequence alignments."""

    def __init__(self, alignment_file: str, format: str = 'fasta'):
        """
        Initialize trimmer with alignment file.

        Args:
            alignment_file: Path to alignment file
            format: Format of alignment (default: fasta)
        """
        self.alignment_file = alignment_file
        self.format = format
        self.alignment = None
        self.original_length = 0
        self.num_sequences = 0

    def load_alignment(self) -> bool:
        """Load alignment from file."""
        try:
            self.alignment = AlignIO.read(self.alignment_file, self.format)
            self.original_length = self.alignment.get_alignment_length()
            self.num_sequences = len(self.alignment)
            return True
        except Exception as e:
            print(f"Error loading alignment: {e}")
            return False

    def get_position_gap_fraction(self) -> List[float]:
        """
        Calculate gap fraction for each alignment position.

        Returns:
            List of gap fractions (0-1) for each position
        """
        gap_fractions = []
        for pos in range(self.alignment.get_alignment_length()):
            column = str(self.alignment[:, pos])
            gap_count = column.count('-')
            gap_fraction = gap_count / len(column)
            gap_fractions.append(gap_fraction)
        return gap_fractions

    def get_position_conservation(self) -> List[float]:
        """
        Calculate Shannon entropy-based conservation for each position.

        Returns:
            List of conservation scores (0-1) for each position
        """
        conservation_scores = []

        for pos in range(self.alignment.get_alignment_length()):
            column = str(self.alignment[:, pos]).upper()

            # Count characters excluding gaps
            char_counts = {}
            total_non_gap = 0

            for char in column:
                if char != '-':
                    char_counts[char] = char_counts.get(char, 0) + 1
                    total_non_gap += 1

            # Calculate Shannon entropy
            if total_non_gap == 0:
                entropy = 0.0
            else:
                entropy = 0.0
                for count in char_counts.values():
                    if count > 0:
                        p = count / total_non_gap
                        entropy -= p * math.log2(p)

            # Convert entropy to conservation (1 - normalized entropy)
            # Maximum entropy: log2(4) for DNA, log2(20) for protein
            max_entropy = math.log2(4)  # conservative estimate
            if max_entropy > 0:
                conservation = 1 - (entropy / max_entropy)
                conservation = max(0, min(1, conservation))  # Clamp to [0, 1]
            else:
                conservation = 0.0

            conservation_scores.append(conservation)

        return conservation_scores

    def get_sequence_gap_fraction(self) -> List[Tuple[str, float]]:
        """
        Calculate gap fraction for each sequence.

        Returns:
            List of (sequence_id, gap_fraction) tuples
        """
        results = []
        for record in self.alignment:
            seq = str(record.seq)
            gap_count = seq.count('-')
            gap_fraction = gap_count / len(seq) if len(seq) > 0 else 0
            results.append((record.id, gap_fraction))
        return results

    def trim_by_gap_threshold(self, max_gap_fraction: float) -> Set[int]:
        """
        Identify positions to keep based on gap threshold.

        Args:
            max_gap_fraction: Maximum allowed gap fraction (0-1)

        Returns:
            Set of position indices to keep
        """
        gap_fractions = self.get_position_gap_fraction()
        positions_to_keep = {i for i, gap_frac in enumerate(gap_fractions)
                            if gap_frac <= max_gap_fraction}
        return positions_to_keep

    def trim_by_conservation(self, min_conservation: float) -> Set[int]:
        """
        Identify positions to keep based on conservation threshold.

        Args:
            min_conservation: Minimum conservation score (0-1)

        Returns:
            Set of position indices to keep
        """
        conservation = self.get_position_conservation()
        positions_to_keep = {i for i, cons in enumerate(conservation)
                            if cons >= min_conservation}
        return positions_to_keep

    def get_contiguous_blocks(self, positions: Set[int], min_block_size: int = 1) -> List[Tuple[int, int]]:
        """
        Find contiguous blocks of positions to keep.

        Args:
            positions: Set of position indices to keep
            min_block_size: Minimum block size to report

        Returns:
            List of (start, end) tuples for contiguous blocks
        """
        if not positions:
            return []

        sorted_positions = sorted(positions)
        blocks = []
        block_start = sorted_positions[0]
        prev_pos = sorted_positions[0]

        for pos in sorted_positions[1:]:
            if pos != prev_pos + 1:  # Gap detected
                if prev_pos - block_start + 1 >= min_block_size:
                    blocks.append((block_start, prev_pos + 1))
                block_start = pos
            prev_pos = pos

        # Add final block
        if prev_pos - block_start + 1 >= min_block_size:
            blocks.append((block_start, prev_pos + 1))

        return blocks

    def trim_sequences_by_gap(self, max_gap_fraction: float) -> Set[int]:
        """
        Identify sequences to keep based on gap threshold.

        Args:
            max_gap_fraction: Maximum allowed gap fraction (0-1)

        Returns:
            Set of sequence indices to keep
        """
        seq_gaps = self.get_sequence_gap_fraction()
        sequences_to_keep = {i for i, (_, gap_frac) in enumerate(seq_gaps)
                            if gap_frac <= max_gap_fraction}
        return sequences_to_keep

    def create_trimmed_alignment(self, position_indices: Set[int],
                                 sequence_indices: Set[int] = None) -> 'Bio.Align.MultipleSeqAlignment':
        """
        Create trimmed alignment from position and sequence indices.

        Args:
            position_indices: Set of positions to keep
            sequence_indices: Set of sequences to keep (None = all)

        Returns:
            Trimmed alignment object
        """
        if sequence_indices is None:
            sequence_indices = set(range(self.num_sequences))

        # Extract trimmed sequences
        trimmed_records = []
        sorted_positions = sorted(position_indices)

        for seq_idx in sorted(sequence_indices):
            record = self.alignment[seq_idx]
            seq = str(record.seq)

            # Extract positions
            trimmed_seq = ''.join(seq[pos] for pos in sorted_positions)

            # Create new record
            new_record = SeqRecord(
                Seq(trimmed_seq),
                id=record.id,
                description=record.description
            )
            trimmed_records.append(new_record)

        # Create alignment
        from Bio.Align import MultipleSeqAlignment
        trimmed_alignment = MultipleSeqAlignment(trimmed_records)

        return trimmed_alignment

    def apply_preset_strategy(self, strategy: str) -> Tuple[Set[int], Set[int]]:
        """
        Apply predefined trimming strategy.

        Args:
            strategy: 'strict', 'moderate', or 'lenient'

        Returns:
            Tuple of (position_indices, sequence_indices) to keep
        """
        if strategy == 'strict':
            # Remove positions with >10% gaps, sequences with >30% gaps
            pos_indices = self.trim_by_gap_threshold(0.10)
            seq_indices = self.trim_sequences_by_gap(0.30)
            return pos_indices, seq_indices

        elif strategy == 'moderate':
            # Remove positions with >30% gaps, sequences with >50% gaps
            pos_indices = self.trim_by_gap_threshold(0.30)
            seq_indices = self.trim_sequences_by_gap(0.50)
            return pos_indices, seq_indices

        elif strategy == 'lenient':
            # Remove positions with >50% gaps, sequences with >70% gaps
            pos_indices = self.trim_by_gap_threshold(0.50)
            seq_indices = self.trim_sequences_by_gap(0.70)
            return pos_indices, seq_indices

        else:
            raise ValueError(f"Unknown strategy: {strategy}")

    def print_trimming_report(self, position_indices: Set[int],
                             sequence_indices: Set[int]):
        """Print trimming statistics."""
        positions_removed = self.original_length - len(position_indices)
        sequences_removed = self.num_sequences - len(sequence_indices)

        print("\n" + "=" * 70)
        print("ALIGNMENT TRIMMING REPORT")
        print("=" * 70)
        print(f"Original alignment:")
        print(f"  Sequences: {self.num_sequences}")
        print(f"  Positions: {self.original_length}")
        print(f"\nTrimmed alignment:")
        print(f"  Sequences: {len(sequence_indices)} (removed {sequences_removed})")
        print(f"  Positions: {len(position_indices)} (removed {positions_removed})")
        print(f"\nRetention rates:")
        print(f"  Positions: {len(position_indices)/self.original_length*100:.1f}%")
        print(f"  Sequences: {len(sequence_indices)/self.num_sequences*100:.1f}%")

        # Sequences removed
        if sequences_removed > 0:
            removed_seqs = [self.alignment[i].id for i in range(self.num_sequences)
                           if i not in sequence_indices]
            print(f"\nSequences removed:")
            for seq_id in removed_seqs[:5]:  # Show first 5
                print(f"  - {seq_id}")
            if len(removed_seqs) > 5:
                print(f"  ... and {len(removed_seqs) - 5} more")

        print("=" * 70 + "\n")

    def save_alignment(self, alignment, output_file: str, format: str = 'fasta') -> bool:
        """
        Save trimmed alignment to file.

        Args:
            alignment: Alignment object
            output_file: Output file path
            format: Output format

        Returns:
            True if successful
        """
        try:
            AlignIO.write(alignment, output_file, format)
            print(f"Trimmed alignment saved to: {output_file}")
            return True
        except Exception as e:
            print(f"Error saving alignment: {e}")
            return False


################################################################################
# Main function
################################################################################

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Trim and clean sequence alignments',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Trimming strategies:
  strict     Remove positions with >10% gaps, sequences with >30% gaps
  moderate   Remove positions with >30% gaps, sequences with >50% gaps
  lenient    Remove positions with >50% gaps, sequences with >70% gaps

Custom parameters override presets.

Examples:
  # Apply strict trimming strategy
  %(prog)s -i alignment.fasta -o trimmed.fasta --strict

  # Remove positions with >50% gaps and sequences with >30% gaps
  %(prog)s -i alignment.fasta -o trimmed.fasta --max-gap 50 --seq-gap 30

  # Keep only well-conserved positions (conservation > 0.7)
  %(prog)s -i alignment.fasta -o trimmed.fasta --min-conservation 0.7

  # Lenient trimming
  %(prog)s -i alignment.fasta -o trimmed.fasta --lenient
        '''
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input alignment file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output alignment file')
    parser.add_argument('-f', '--format', default='fasta',
                        choices=['fasta', 'stockholm', 'clustal', 'phylip'],
                        help='Alignment format (default: fasta)')

    # Strategy options
    strategy_group = parser.add_argument_group('Trimming strategies (presets)')
    strategy_group.add_argument('--strict', action='store_true',
                               help='Apply strict trimming')
    strategy_group.add_argument('--moderate', action='store_true',
                               help='Apply moderate trimming')
    strategy_group.add_argument('--lenient', action='store_true',
                               help='Apply lenient trimming')

    # Custom parameters
    param_group = parser.add_argument_group('Custom parameters')
    param_group.add_argument('--max-gap', type=float, metavar='PCT',
                            help='Max gap percentage in positions (0-100)')
    param_group.add_argument('--seq-gap', type=float, metavar='PCT',
                            help='Max gap percentage in sequences (0-100)')
    param_group.add_argument('--min-conservation', type=float, metavar='SCORE',
                            help='Minimum conservation score (0-1)')

    parser.add_argument('--verbose', action='store_true',
                        help='Print detailed trimming statistics')

    args = parser.parse_args()

    # Check if input file exists
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    # Determine trimming strategy
    strategy = None
    if args.strict:
        strategy = 'strict'
    elif args.moderate:
        strategy = 'moderate'
    elif args.lenient:
        strategy = 'lenient'

    # Check for conflicting custom parameters
    has_custom_params = any([args.max_gap is not None, args.seq_gap is not None,
                            args.min_conservation is not None])

    if strategy and has_custom_params:
        print("Warning: Custom parameters override preset strategy")

    # Load alignment
    trimmer = AlignmentTrimmer(args.input, args.format)
    if not trimmer.load_alignment():
        sys.exit(1)

    print(f"Loaded alignment: {trimmer.num_sequences} sequences, {trimmer.original_length} bp")

    # Apply trimming
    if has_custom_params:
        # Custom trimming
        position_indices = set(range(trimmer.original_length))
        sequence_indices = set(range(trimmer.num_sequences))

        if args.max_gap is not None:
            max_gap = args.max_gap / 100.0
            position_indices = position_indices.intersection(
                trimmer.trim_by_gap_threshold(max_gap)
            )
            print(f"Trimming positions with >{args.max_gap}% gaps...")

        if args.seq_gap is not None:
            max_seq_gap = args.seq_gap / 100.0
            sequence_indices = sequence_indices.intersection(
                trimmer.trim_sequences_by_gap(max_seq_gap)
            )
            print(f"Trimming sequences with >{args.seq_gap}% gaps...")

        if args.min_conservation is not None:
            position_indices = position_indices.intersection(
                trimmer.trim_by_conservation(args.min_conservation)
            )
            print(f"Trimming positions with conservation <{args.min_conservation}...")

    elif strategy:
        # Preset strategy
        print(f"Applying {strategy} trimming strategy...")
        position_indices, sequence_indices = trimmer.apply_preset_strategy(strategy)

    else:
        # Default: moderate strategy
        print("Applying moderate trimming strategy (default)...")
        position_indices, sequence_indices = trimmer.apply_preset_strategy('moderate')

    # Create trimmed alignment
    trimmed_alignment = trimmer.create_trimmed_alignment(position_indices, sequence_indices)

    # Print report
    trimmer.print_trimming_report(position_indices, sequence_indices)

    # Save trimmed alignment
    if trimmer.save_alignment(trimmed_alignment, args.output, args.format):
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()
