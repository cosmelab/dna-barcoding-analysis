#!/usr/bin/env python3
"""
FASTA File Manipulation Tools

This module provides functions to parse, split, merge, and filter FASTA files.
FASTA format is a text-based format for nucleotide and protein sequences.

Example FASTA format:
    >sequence_id description
    ATCGATCGATCGATCG
    ATCGATCGATCGATCG
    >another_sequence
    GCTAGCTAGCTAGCTA

Author: DNA Barcoding Analysis Course
License: MIT
"""

import sys
from pathlib import Path
from typing import List, Tuple, Dict, Iterator


def read_fasta(filename: str) -> Iterator[Tuple[str, str]]:
    """
    Parse a FASTA file and yield (header, sequence) tuples.

    This function reads a FASTA file line by line and yields sequence records.
    Memory-efficient for large files as it uses generators.

    Args:
        filename (str): Path to the FASTA file

    Yields:
        Tuple[str, str]: (header without '>', sequence)

    Example:
        >>> for header, seq in read_fasta('sequences.fasta'):
        ...     print(f"{header}: {len(seq)} bp")
    """
    header = None
    sequence = []

    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.rstrip('\n')

                if line.startswith('>'):
                    # If we have a previous sequence, yield it
                    if header is not None:
                        yield header, ''.join(sequence)

                    # Start new sequence
                    header = line[1:]  # Remove '>' character
                    sequence = []
                else:
                    # Add to current sequence
                    if header is not None:
                        sequence.append(line)

            # Don't forget the last sequence
            if header is not None:
                yield header, ''.join(sequence)

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)


def write_fasta(filename: str, sequences: List[Tuple[str, str]],
                line_width: int = 80) -> int:
    """
    Write sequences to a FASTA file.

    Sequences are wrapped at a specified line width for readability.

    Args:
        filename (str): Output file path
        sequences (List[Tuple[str, str]]): List of (header, sequence) tuples
        line_width (int): Characters per line (default: 80)

    Returns:
        int: Number of sequences written

    Example:
        >>> sequences = [('seq1', 'ATCGATCG'), ('seq2', 'GCTAGCTA')]
        >>> write_fasta('output.fasta', sequences)
        2
    """
    try:
        with open(filename, 'w') as f:
            count = 0
            for header, sequence in sequences:
                # Write header
                f.write(f">{header}\n")

                # Write sequence in chunks
                for i in range(0, len(sequence), line_width):
                    f.write(sequence[i:i+line_width] + '\n')

                count += 1

        print(f"Successfully wrote {count} sequences to '{filename}'", file=sys.stderr)
        return count

    except IOError as e:
        print(f"Error writing to file: {e}", file=sys.stderr)
        sys.exit(1)


def split_fasta(input_file: str, output_prefix: str,
                sequences_per_file: int = 1000) -> Dict[str, int]:
    """
    Split a large FASTA file into smaller chunks.

    Useful for processing large sequence databases in parallel.

    Args:
        input_file (str): Path to input FASTA file
        output_prefix (str): Prefix for output files (e.g., 'split' -> split_001.fasta)
        sequences_per_file (int): Number of sequences per output file

    Returns:
        Dict[str, int]: Mapping of output filenames to sequence counts

    Example:
        >>> results = split_fasta('large.fasta', 'chunks', 500)
        >>> print(results)
        {'chunks_001.fasta': 500, 'chunks_002.fasta': 500, 'chunks_003.fasta': 200}
    """
    output_files = {}
    current_file_num = 1
    current_sequences = []
    total_seqs = 0

    try:
        for header, sequence in read_fasta(input_file):
            current_sequences.append((header, sequence))

            # Write file when we reach the limit
            if len(current_sequences) >= sequences_per_file:
                output_file = f"{output_prefix}_{current_file_num:03d}.fasta"
                write_fasta(output_file, current_sequences)
                output_files[output_file] = len(current_sequences)

                current_sequences = []
                current_file_num += 1
                total_seqs += sequences_per_file

        # Write remaining sequences
        if current_sequences:
            output_file = f"{output_prefix}_{current_file_num:03d}.fasta"
            write_fasta(output_file, current_sequences)
            output_files[output_file] = len(current_sequences)

        return output_files

    except Exception as e:
        print(f"Error splitting FASTA file: {e}", file=sys.stderr)
        sys.exit(1)


def merge_fasta(input_files: List[str], output_file: str) -> int:
    """
    Merge multiple FASTA files into a single file.

    Args:
        input_files (List[str]): List of input FASTA file paths
        output_file (str): Path for merged output file

    Returns:
        int: Total number of sequences merged

    Example:
        >>> merge_fasta(['file1.fasta', 'file2.fasta'], 'merged.fasta')
        1000
    """
    all_sequences = []

    try:
        for input_file in input_files:
            print(f"Reading {input_file}...", file=sys.stderr)
            for header, sequence in read_fasta(input_file):
                all_sequences.append((header, sequence))

        count = write_fasta(output_file, all_sequences)
        return count

    except Exception as e:
        print(f"Error merging FASTA files: {e}", file=sys.stderr)
        sys.exit(1)


def filter_fasta(input_file: str, output_file: str,
                min_length: int = 0, max_length: int = float('inf'),
                header_pattern: str = None) -> int:
    """
    Filter FASTA sequences by length or header pattern.

    Args:
        input_file (str): Input FASTA file
        output_file (str): Output FASTA file
        min_length (int): Minimum sequence length (default: 0)
        max_length (int): Maximum sequence length (default: no limit)
        header_pattern (str): Regex pattern to match in header (optional)

    Returns:
        int: Number of sequences that passed filters

    Example:
        >>> filter_fasta('seqs.fasta', 'filtered.fasta', min_length=100)
        250
    """
    import re

    filtered_sequences = []
    pattern = re.compile(header_pattern) if header_pattern else None

    try:
        for header, sequence in read_fasta(input_file):
            seq_len = len(sequence)

            # Check length criteria
            if seq_len < min_length or seq_len > max_length:
                continue

            # Check header pattern if provided
            if pattern and not pattern.search(header):
                continue

            filtered_sequences.append((header, sequence))

        count = write_fasta(output_file, filtered_sequences)
        return count

    except Exception as e:
        print(f"Error filtering FASTA file: {e}", file=sys.stderr)
        sys.exit(1)


def count_sequences(filename: str) -> int:
    """
    Count the number of sequences in a FASTA file.

    Args:
        filename (str): Input FASTA file

    Returns:
        int: Number of sequences

    Example:
        >>> count_sequences('data.fasta')
        1042
    """
    count = 0
    try:
        for _ in read_fasta(filename):
            count += 1
        return count
    except Exception as e:
        print(f"Error counting sequences: {e}", file=sys.stderr)
        return 0


def get_sequence_by_id(filename: str, seq_id: str) -> str:
    """
    Extract a single sequence by its ID from a FASTA file.

    Args:
        filename (str): Input FASTA file
        seq_id (str): Sequence ID to search for

    Returns:
        str: The sequence, or empty string if not found

    Example:
        >>> seq = get_sequence_by_id('data.fasta', 'seq001')
        >>> print(seq)
        ATCGATCGATCGATCG
    """
    try:
        for header, sequence in read_fasta(filename):
            # Check if seq_id matches the header (exact match or first word)
            if header == seq_id or header.split()[0] == seq_id:
                return sequence
        return ""
    except Exception as e:
        print(f"Error retrieving sequence: {e}", file=sys.stderr)
        return ""


if __name__ == "__main__":
    # Example usage and basic tests
    print("FASTA Tools Module")
    print("This module provides functions for manipulating FASTA files.")
    print("\nAvailable functions:")
    print("  - read_fasta(filename): Iterator over (header, sequence) tuples")
    print("  - write_fasta(filename, sequences): Write sequences to file")
    print("  - split_fasta(input, prefix, count): Split file into chunks")
    print("  - merge_fasta(inputs, output): Merge multiple files")
    print("  - filter_fasta(input, output, min_len, max_len, pattern): Filter sequences")
    print("  - count_sequences(filename): Count total sequences")
    print("  - get_sequence_by_id(filename, seq_id): Get specific sequence")
