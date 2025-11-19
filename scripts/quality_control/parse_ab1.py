#!/usr/bin/env python3
"""
parse_ab1.py - Parse ABI .ab1 DNA sequencing files

Educational example: This script demonstrates how to read and extract
information from ABI Sanger sequencing chromatogram files (.ab1 format).

The .ab1 format is a binary file format used by Applied Biosystems
DNA sequencers to store:
  - The DNA sequence (basecalls)
  - Quality scores for each base
  - Raw trace data (fluorescence intensity for each dye channel)
  - Sequencing metadata (machine model, date, sample name, etc.)

This script uses BioPython's SeqIO module, which handles the binary
parsing automatically. This is much easier than parsing .ab1 files manually!

Usage:
    python parse_ab1.py sample.ab1
    python parse_ab1.py sample.ab1 --save-fasta sample.fasta
    python parse_ab1.py sample.ab1 --save-qual sample.qual

Requirements:
    BioPython: pip install biopython
"""

import sys
import argparse
from pathlib import Path
from Bio import SeqIO


def parse_ab1_file(ab1_file):
    """
    Parse a single .ab1 file using BioPython.

    Args:
        ab1_file (str): Path to the .ab1 file

    Returns:
        list: List of SeqRecord objects (usually just one for single-read files)
    """
    # SeqIO.parse() is a generator that yields SeqRecord objects
    # The file must be opened in binary mode ('rb') because .ab1 is binary
    records = []

    try:
        with open(ab1_file, 'rb') as f:
            # The 'abi' format tells BioPython it's an ABI file
            for record in SeqIO.parse(f, 'abi'):
                records.append(record)
    except FileNotFoundError:
        print(f"Error: File '{ab1_file}' not found", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error reading {ab1_file}: {e}", file=sys.stderr)
        return None

    return records


def print_sequence_info(record):
    """
    Print extracted information about a sequence record.

    Args:
        record: A SeqRecord object from BioPython
    """
    # Basic sequence information
    print(f"\n{'='*70}")
    print(f"SEQUENCE INFORMATION")
    print(f"{'='*70}")
    print(f"Sequence ID: {record.id}")
    print(f"Description: {record.description}")
    print(f"Sequence length: {len(record.seq)} bp")
    print(f"Sequence: {str(record.seq)[:100]}...")

    # Quality scores are stored in annotations as a list of integers
    if 'phred_quality' in record.letter_annotations:
        quals = record.letter_annotations['phred_quality']
        print(f"\n{'='*70}")
        print(f"QUALITY INFORMATION")
        print(f"{'='*70}")
        print(f"Number of quality scores: {len(quals)}")
        print(f"Min quality: {min(quals)}")
        print(f"Max quality: {max(quals)}")
        print(f"Average quality: {sum(quals)/len(quals):.1f}")
        print(f"Quality scores (first 50): {quals[:50]}")


def print_metadata(record):
    """
    Print metadata extracted from the .ab1 file.

    The annotations dictionary contains various metadata from the file,
    including machine information, sample details, and processing settings.

    Args:
        record: A SeqRecord object from BioPython
    """
    print(f"\n{'='*70}")
    print(f"METADATA INFORMATION")
    print(f"{'='*70}")

    # Common metadata keys in AB1 files:
    metadata_keys = [
        'machine_model',
        'sample_well',
        'run_start',
        'run_finish',
        'polymer',
        'dye',
    ]

    for key in metadata_keys:
        if key in record.annotations:
            print(f"{key}: {record.annotations[key]}")

    # Show all available annotations if there are others
    other_keys = set(record.annotations.keys()) - set(metadata_keys)
    if other_keys:
        print(f"\nOther available metadata: {', '.join(sorted(other_keys))}")


def save_fasta(record, output_file):
    """
    Save the sequence to a FASTA file.

    FASTA format is a simple text format for DNA sequences:
    >sequence_id description
    ATCGATCGATCGATCG...

    Args:
        record: A SeqRecord object
        output_file (str): Path to save the FASTA file
    """
    SeqIO.write(record, output_file, 'fasta')
    print(f"\nSequence saved to FASTA: {output_file}")


def save_qual(record, output_file):
    """
    Save quality scores to a QUAL file.

    QUAL format is a simple text format for quality scores:
    >sequence_id description
    30 28 32 31 29 30...

    This format is commonly used with FASTA for sequence quality information.

    Args:
        record: A SeqRecord object
        output_file (str): Path to save the QUAL file
    """
    # Extract quality scores from the record
    if 'phred_quality' in record.letter_annotations:
        quals = record.letter_annotations['phred_quality']

        # Format as QUAL file
        with open(output_file, 'w') as f:
            # Write header line (same as FASTA header)
            f.write(f">{record.id} {record.description}\n")

            # Write quality scores space-separated, 50 per line
            for i in range(0, len(quals), 50):
                qual_line = ' '.join(str(q) for q in quals[i:i+50])
                f.write(qual_line + '\n')

        print(f"Quality scores saved to QUAL: {output_file}")
    else:
        print("Warning: No quality scores found in sequence", file=sys.stderr)


def main():
    """
    Main function - Parse arguments and process .ab1 file.
    """
    parser = argparse.ArgumentParser(
        description='Parse ABI .ab1 DNA sequencing files',
        epilog='Example: python parse_ab1.py sample.ab1 --save-fasta output.fasta'
    )

    parser.add_argument(
        'ab1_file',
        help='Path to the .ab1 file to parse'
    )

    parser.add_argument(
        '--save-fasta',
        metavar='FILE',
        help='Save sequence to FASTA file'
    )

    parser.add_argument(
        '--save-qual',
        metavar='FILE',
        help='Save quality scores to QUAL file'
    )

    args = parser.parse_args()

    # Parse the .ab1 file
    print(f"Parsing {args.ab1_file}...")
    records = parse_ab1_file(args.ab1_file)

    if not records:
        sys.exit(1)

    # Process the first record (most .ab1 files have only one)
    record = records[0]

    # Print information about the sequence
    print_sequence_info(record)

    # Print metadata
    print_metadata(record)

    # Save outputs if requested
    if args.save_fasta:
        save_fasta(record, args.save_fasta)

    if args.save_qual:
        save_qual(record, args.save_qual)


if __name__ == '__main__':
    main()
