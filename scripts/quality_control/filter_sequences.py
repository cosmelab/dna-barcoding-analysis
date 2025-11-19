#!/usr/bin/env python3
"""
filter_sequences.py - Filter DNA sequences based on quality metrics

Educational example: This script demonstrates sequence filtering, which removes
entire sequences that don't meet quality standards.

This is different from trimming - instead of removing low-quality bases from
the ends of sequences, filtering removes entire sequences that are:
- Too short (< minimum length)
- Have too many low-quality bases
- Have too high ambiguity (N characters)
- Have low average quality

When to use filtering:
- Quality control of raw sequencing data
- Before assembly or alignment
- To improve downstream analysis reliability
- To set sequence quality standards

Common filtering criteria:
- Minimum length: 300-400 bp (depends on application)
- Minimum average quality: Q20-Q30
- Maximum N characters: 0-5%
- Maximum low-quality bases: 5-10%

Usage:
    # Filter FASTQ with default parameters
    python filter_sequences.py input.fastq --output filtered.fastq

    # Strict filtering (require Q25, min 400 bp)
    python filter_sequences.py input.fastq --min-quality 25 --min-length 400

    # Filter with multiple criteria
    python filter_sequences.py input.fastq \\
        --min-quality 20 \\
        --min-length 300 \\
        --max-n 5 \\
        --output filtered.fastq

Requirements:
    BioPython: pip install biopython
"""

import sys
import argparse
from Bio import SeqIO


class SequenceFilter:
    """
    A class to filter sequences based on multiple quality criteria.

    This encapsulates the filtering logic and makes it reusable.
    """

    def __init__(self, min_length=0, min_quality=0, max_n_percent=100):
        """
        Initialize the filter with criteria.

        Args:
            min_length (int): Minimum sequence length in bp
            min_quality (float): Minimum average PHRED quality score
            max_n_percent (float): Maximum percentage of N characters allowed (0-100)
        """
        self.min_length = min_length
        self.min_quality = min_quality
        self.max_n_percent = max_n_percent

        # Statistics tracking
        self.stats = {
            'total_sequences': 0,
            'passed_sequences': 0,
            'failed_length': 0,
            'failed_quality': 0,
            'failed_ambiguity': 0,
            'bases_in': 0,
            'bases_out': 0,
        }

    def check_length(self, record):
        """
        Check if sequence meets minimum length requirement.

        Args:
            record: SeqRecord object
            returns: True if passes, False otherwise
        """
        return len(record.seq) >= self.min_length

    def check_quality(self, record):
        """
        Check if sequence meets minimum average quality requirement.

        Args:
            record: SeqRecord object with quality scores
            returns: Tuple (passes: bool, avg_quality: float)
        """
        if 'phred_quality' not in record.letter_annotations:
            # No quality scores - cannot evaluate
            return True, None

        qualities = record.letter_annotations['phred_quality']

        if len(qualities) == 0:
            return False, 0.0

        avg_quality = sum(qualities) / len(qualities)

        return avg_quality >= self.min_quality, avg_quality

    def check_ambiguity(self, record):
        """
        Check if sequence meets maximum N character threshold.

        N represents an ambiguous base (couldn't be called as A/T/G/C).
        High N counts suggest poor sequencing quality.

        Args:
            record: SeqRecord object
            returns: Tuple (passes: bool, n_percent: float)
        """
        seq_str = str(record.seq).upper()
        n_count = seq_str.count('N')
        n_percent = (n_count / len(seq_str)) * 100 if len(seq_str) > 0 else 0

        return n_percent <= self.max_n_percent, n_percent

    def filter(self, record):
        """
        Apply all filter criteria to a sequence record.

        Args:
            record: SeqRecord object

        Returns:
            Tuple (passes: bool, reason: str or None)
                  reason is None if passes, otherwise explains why it failed
        """
        self.stats['total_sequences'] += 1
        self.stats['bases_in'] += len(record.seq)

        # Check length
        if not self.check_length(record):
            self.stats['failed_length'] += 1
            return False, 'length'

        # Check quality
        passes_quality, avg_quality = self.check_quality(record)
        if not passes_quality:
            self.stats['failed_quality'] += 1
            return False, 'quality'

        # Check ambiguity
        passes_ambiguity, n_percent = self.check_ambiguity(record)
        if not passes_ambiguity:
            self.stats['failed_ambiguity'] += 1
            return False, 'ambiguity'

        # All checks passed
        self.stats['passed_sequences'] += 1
        self.stats['bases_out'] += len(record.seq)
        return True, None

    def get_statistics(self):
        """
        Get a summary of filtering statistics.

        Returns:
            dict: Statistics dictionary
        """
        return self.stats


def filter_sequences(input_file, filter_obj, input_format='fastq'):
    """
    Apply filters to all sequences in an input file.

    Args:
        input_file (str): Path to input sequence file
        filter_obj (SequenceFilter): Configured filter object
        input_format (str): 'fastq' or 'fasta'

    Returns:
        list: List of (record, filter_result) tuples for passed sequences
    """
    results = []

    try:
        with open(input_file, 'r') as f:
            for record in SeqIO.parse(f, input_format):
                passes, reason = filter_obj.filter(record)

                if passes:
                    results.append(record)
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error processing sequences: {e}", file=sys.stderr)
        return None

    return results


def save_filtered_sequences(records, output_file, output_format='fastq'):
    """
    Save filtered sequences to output file.

    Args:
        records (list): List of SeqRecord objects
        output_file (str): Path to output file
        output_format (str): 'fastq' or 'fasta'
    """
    if not records:
        print("No sequences passed filtering", file=sys.stderr)
        return False

    try:
        SeqIO.write(records, output_file, output_format)
        print(f"\nFiltered sequences saved to {output_file}")
        return True
    except Exception as e:
        print(f"Error saving sequences: {e}", file=sys.stderr)
        return False


def print_statistics(stats):
    """
    Print detailed filtering statistics.

    Args:
        stats (dict): Statistics dictionary from SequenceFilter
    """
    total = stats['total_sequences']

    print(f"\n{'='*70}")
    print(f"FILTERING STATISTICS")
    print(f"{'='*70}")
    print(f"\nSequence counts:")
    print(f"  Total sequences: {total}")
    print(f"  Passed filtering: {stats['passed_sequences']} ({100*stats['passed_sequences']/max(1,total):.1f}%)")
    print(f"  Failed - too short: {stats['failed_length']}")
    print(f"  Failed - low quality: {stats['failed_quality']}")
    print(f"  Failed - too many Ns: {stats['failed_ambiguity']}")

    print(f"\nBase counts:")
    print(f"  Bases in (input): {stats['bases_in']:,}")
    print(f"  Bases out (output): {stats['bases_out']:,}")
    print(f"  Bases filtered: {stats['bases_in'] - stats['bases_out']:,}")

    if stats['passed_sequences'] > 0:
        avg_length = stats['bases_out'] / stats['passed_sequences']
        print(f"  Average length (output): {avg_length:.0f} bp")


def main():
    """
    Main function - Parse arguments and filter sequences.
    """
    parser = argparse.ArgumentParser(
        description='Filter DNA sequences based on quality metrics',
        epilog='Example: python filter_sequences.py input.fastq --min-quality 20 --min-length 300 --output filtered.fastq'
    )

    parser.add_argument(
        'input_file',
        help='Input sequence file (FASTQ or FASTA)'
    )

    parser.add_argument(
        '--min-length', '-l',
        type=int,
        default=0,
        help='Minimum sequence length in bp (default: 0)'
    )

    parser.add_argument(
        '--min-quality', '-q',
        type=float,
        default=0,
        help='Minimum average quality score (default: 0)'
    )

    parser.add_argument(
        '--max-n', '-n',
        type=float,
        default=100,
        help='Maximum percentage of N characters allowed (default: 100)'
    )

    parser.add_argument(
        '--format', '-f',
        choices=['fastq', 'fasta'],
        default='fastq',
        help='Input file format (default: fastq)'
    )

    parser.add_argument(
        '--output', '-o',
        help='Output file (defaults to input_filtered.fastq/fasta)'
    )

    args = parser.parse_args()

    # Validate parameters
    if args.min_length < 0:
        print("Error: Minimum length must be >= 0", file=sys.stderr)
        sys.exit(1)

    if not 0 <= args.min_quality <= 40:
        print("Error: Minimum quality must be between 0 and 40", file=sys.stderr)
        sys.exit(1)

    if not 0 <= args.max_n <= 100:
        print("Error: Maximum N percentage must be between 0 and 100", file=sys.stderr)
        sys.exit(1)

    # Determine output filename
    if args.output is None:
        base = args.input_file.rsplit('.', 1)[0]
        args.output = f"{base}_filtered.{args.format}"

    # Print filter settings
    print(f"{'='*70}")
    print(f"FILTERING SETTINGS")
    print(f"{'='*70}")
    print(f"Input file: {args.input_file}")
    print(f"Output file: {args.output}")
    print(f"\nFilter criteria:")
    print(f"  Minimum length: {args.min_length} bp")
    print(f"  Minimum average quality: {args.min_quality}")
    print(f"  Maximum N percentage: {args.max_n}%\n")

    # Create filter object
    sequence_filter = SequenceFilter(
        min_length=args.min_length,
        min_quality=args.min_quality,
        max_n_percent=args.max_n
    )

    # Filter sequences
    filtered_records = filter_sequences(args.input_file, sequence_filter, args.format)

    if filtered_records is None:
        sys.exit(1)

    # Save filtered sequences
    if not save_filtered_sequences(filtered_records, args.output, args.format):
        sys.exit(1)

    # Print statistics
    print_statistics(sequence_filter.get_statistics())


if __name__ == '__main__':
    main()
