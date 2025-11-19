#!/usr/bin/env python3
"""
trim_quality.py - Trim low-quality bases from sequence ends

Educational example: This script demonstrates quality-based trimming,
which removes low-confidence basecalls from the beginning and end of
a DNA sequence.

Why trim low quality bases?
- DNA sequencers often have worse accuracy at the start and end of reads
- Low quality bases increase error rates in downstream analyses
- Trimming improves alignment and assembly quality

How it works:
1. Reads in FASTA format sequences with associated quality scores (QUAL or FASTQ)
2. Finds the first base with quality >= threshold (from left)
3. Finds the last base with quality >= threshold (from right)
4. Keeps only the sequence between these points
5. Saves trimmed sequences to output file

Quality score scale (PHRED quality):
- Q20 = 1% error rate (99% accuracy) - commonly used minimum
- Q25 = 0.3% error rate
- Q30 = 0.1% error rate

Usage:
    # Trim FASTQ file with Q20 threshold
    python trim_quality.py input.fastq --quality 20 --output trimmed.fastq

    # Trim FASTA+QUAL file
    python trim_quality.py input.fasta --qual input.qual --quality 20

    # Strict trimming with Q30
    python trim_quality.py input.fastq --quality 30

Requirements:
    BioPython: pip install biopython
"""

import sys
import argparse
from collections import defaultdict
from Bio import SeqIO


def trim_sequence(record, quality_threshold):
    """
    Trim low-quality bases from both ends of a sequence.

    This function finds the first and last position with quality >= threshold,
    then returns a trimmed record containing only that region.

    Args:
        record: A SeqRecord object (must have quality scores)
        quality_threshold (int): Minimum quality score to keep

    Returns:
        tuple: (trimmed_record, trim_start, trim_end, bases_removed)
               or (None, 0, 0, 0) if sequence is too short after trimming
    """
    # Get quality scores from the record
    if 'phred_quality' not in record.letter_annotations:
        print(f"Warning: {record.id} has no quality scores, skipping", file=sys.stderr)
        return None, 0, 0, 0

    qualities = record.letter_annotations['phred_quality']

    # Check that we have sequence and quality data
    if len(qualities) == 0:
        return None, 0, 0, 0

    # Find the FIRST position with quality >= threshold (from left)
    # Start at position 0 by default (no trimming from left)
    trim_start = 0
    for i, qual in enumerate(qualities):
        if qual >= quality_threshold:
            trim_start = i
            break
    else:
        # All bases are below threshold
        return None, 0, 0, len(qualities)

    # Find the LAST position with quality >= threshold (from right)
    # We search backwards from the end
    trim_end = len(qualities)
    for i in range(len(qualities) - 1, -1, -1):
        if qualities[i] >= quality_threshold:
            trim_end = i + 1  # +1 because slicing is exclusive at the end
            break

    # Check if trimmed sequence would be too short (less than 50 bp)
    trimmed_length = trim_end - trim_start
    if trimmed_length < 50:
        # Sequence is too short to be useful
        return None, trim_start, trim_end, len(qualities)

    # Create a new trimmed record
    trimmed_record = record[trim_start:trim_end]

    # Calculate how many bases were removed
    bases_removed = trim_start + (len(qualities) - trim_end)

    return trimmed_record, trim_start, trim_end, bases_removed


def process_sequences(input_file, quality_threshold, input_format='fastq', qual_file=None):
    """
    Process a sequence file and trim all sequences.

    Args:
        input_file (str): Path to input sequence file
        quality_threshold (int): Minimum quality score
        input_format (str): 'fastq' or 'fasta'
        qual_file (str): Path to .qual file (only for FASTA+QUAL format)

    Returns:
        list: List of (trimmed_record, stats_dict) tuples
    """
    results = []
    stats = {
        'total_sequences': 0,
        'passed_sequences': 0,
        'failed_sequences': 0,
        'bases_before': 0,
        'bases_after': 0,
    }

    try:
        # Parse sequences
        sequences = SeqIO.to_dict(SeqIO.parse(input_file, input_format))

        # If using separate QUAL file, merge quality information
        if qual_file and input_format == 'fasta':
            qual_dict = SeqIO.to_dict(SeqIO.parse(qual_file, 'qual'))

            for seq_id, record in sequences.items():
                if seq_id in qual_dict:
                    qual_record = qual_dict[seq_id]
                    # Add quality scores to the sequence record
                    record.letter_annotations['phred_quality'] = qual_record.letter_annotations['phred_quality']

        # Process each sequence
        for record in sequences.values():
            stats['total_sequences'] += 1
            stats['bases_before'] += len(record.seq)

            trimmed, start, end, removed = trim_sequence(record, quality_threshold)

            if trimmed:
                stats['passed_sequences'] += 1
                stats['bases_after'] += len(trimmed.seq)
                results.append((trimmed, {
                    'original_length': len(record.seq),
                    'trimmed_length': len(trimmed.seq),
                    'trim_start': start,
                    'trim_end': end,
                    'bases_removed': removed,
                }))
            else:
                stats['failed_sequences'] += 1

    except FileNotFoundError as e:
        print(f"Error: Could not open file {e.filename}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error processing sequences: {e}", file=sys.stderr)
        return None

    stats['total_results'] = results
    return stats


def save_sequences(results, output_file, output_format='fastq'):
    """
    Save trimmed sequences to output file.

    Args:
        results (dict): Results dictionary from process_sequences()
        output_file (str): Path to output file
        output_format (str): 'fastq' or 'fasta'
    """
    if not results or 'total_results' not in results:
        print("No sequences to save", file=sys.stderr)
        return

    trimmed_records = [r[0] for r in results['total_results']]

    try:
        SeqIO.write(trimmed_records, output_file, output_format)
        print(f"\nTrimmed sequences saved to {output_file}")
    except Exception as e:
        print(f"Error saving sequences: {e}", file=sys.stderr)


def print_statistics(stats):
    """
    Print trimming statistics.

    Args:
        stats (dict): Statistics dictionary from process_sequences()
    """
    if not stats:
        return

    print(f"\n{'='*70}")
    print(f"TRIMMING STATISTICS")
    print(f"{'='*70}")
    print(f"Total sequences: {stats['total_sequences']}")
    print(f"Passed trimming: {stats['passed_sequences']}")
    print(f"Failed trimming (too short): {stats['failed_sequences']}")
    print(f"Pass rate: {100*stats['passed_sequences']/max(1, stats['total_sequences']):.1f}%")

    print(f"\nBase statistics:")
    print(f"Bases before trimming: {stats['bases_before']:,}")
    print(f"Bases after trimming: {stats['bases_after']:,}")
    print(f"Bases removed: {stats['bases_before'] - stats['bases_after']:,}")

    if stats['total_sequences'] > 0:
        avg_before = stats['bases_before'] / stats['total_sequences']
        avg_after = stats['bases_after'] / max(1, stats['passed_sequences'])
        print(f"\nAverage sequence length before: {avg_before:.0f} bp")
        print(f"Average sequence length after: {avg_after:.0f} bp")


def main():
    """
    Main function - Parse arguments and trim sequences.
    """
    parser = argparse.ArgumentParser(
        description='Trim low-quality bases from sequence ends',
        epilog='Example: python trim_quality.py input.fastq --quality 20 --output trimmed.fastq'
    )

    parser.add_argument(
        'input_file',
        help='Input sequence file (FASTQ or FASTA)'
    )

    parser.add_argument(
        '--quality', '-q',
        type=int,
        default=20,
        help='Minimum quality threshold (default: 20, range: 0-40)'
    )

    parser.add_argument(
        '--format', '-f',
        choices=['fastq', 'fasta'],
        default='fastq',
        help='Input file format (default: fastq)'
    )

    parser.add_argument(
        '--qual',
        help='Quality file for FASTA format (separate .qual file)'
    )

    parser.add_argument(
        '--output', '-o',
        help='Output file (defaults to input_trimmed.fastq/fasta)'
    )

    args = parser.parse_args()

    # Validate quality threshold
    if not 0 <= args.quality <= 40:
        print("Error: Quality threshold must be between 0 and 40", file=sys.stderr)
        sys.exit(1)

    # Determine output filename
    if args.output is None:
        base = args.input_file.rsplit('.', 1)[0]
        args.output = f"{base}_trimmed.{args.format}"

    # Process sequences
    print(f"Trimming sequences with quality >= {args.quality}...")
    print(f"Input: {args.input_file}")
    print(f"Output: {args.output}\n")

    stats = process_sequences(
        args.input_file,
        args.quality,
        input_format=args.format,
        qual_file=args.qual
    )

    if stats is None:
        sys.exit(1)

    # Save trimmed sequences
    save_sequences(stats, args.output, output_format=args.format)

    # Print statistics
    print_statistics(stats)


if __name__ == '__main__':
    main()
