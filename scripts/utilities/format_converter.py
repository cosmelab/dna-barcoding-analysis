#!/usr/bin/env python3
"""
Sequence Format Converter

This module converts between common bioinformatics sequence formats:
- FASTA: Simple text format with header and sequence
- FASTQ: FASTA + quality scores
- GenBank: Comprehensive format with annotations and metadata
- PHYLIP: Phylogenetics format with fixed-width sequences

Format specifications:

FASTA:
    >sequence_id description
    ATCGATCGATCGATCG
    ATCGATCGATCGATCG

FASTQ:
    @sequence_id description
    ATCGATCGATCGATCG
    +
    IIIIIIIIIIIIIIII

GenBank: Complex format with features, source, organism info

PHYLIP:
    4 100
    seq1       ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    seq2       ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    seq3       ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    seq4       ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG

Author: DNA Barcoding Analysis Course
License: MIT
"""

import sys
import re
from pathlib import Path
from typing import List, Tuple, Dict, Iterator
from dataclasses import dataclass


@dataclass
class Sequence:
    """
    Represents a biological sequence with metadata.

    Attributes:
        header (str): Sequence identifier and description
        sequence (str): The actual sequence
        quality (str): Quality scores (optional, for FASTQ)
        features (dict): Additional metadata/annotations (optional)
    """
    header: str
    sequence: str
    quality: str = None
    features: dict = None


class FASTAHandler:
    """Read and write FASTA format files."""

    @staticmethod
    def read(filename: str) -> Iterator[Sequence]:
        """
        Parse a FASTA file.

        Args:
            filename (str): Path to FASTA file

        Yields:
            Sequence: Sequence objects
        """
        try:
            with open(filename, 'r') as f:
                header = None
                sequence = []

                for line in f:
                    line = line.rstrip('\n')

                    if line.startswith('>'):
                        if header is not None:
                            yield Sequence(header, ''.join(sequence))

                        header = line[1:]
                        sequence = []
                    else:
                        if header is not None:
                            sequence.append(line)

                if header is not None:
                    yield Sequence(header, ''.join(sequence))

        except FileNotFoundError:
            print(f"Error: File '{filename}' not found.", file=sys.stderr)
            sys.exit(1)

    @staticmethod
    def write(filename: str, sequences: List[Sequence], line_width: int = 80):
        """
        Write sequences to FASTA file.

        Args:
            filename (str): Output file path
            sequences (List[Sequence]): List of Sequence objects
            line_width (int): Characters per line (default: 80)
        """
        try:
            with open(filename, 'w') as f:
                for seq in sequences:
                    f.write(f">{seq.header}\n")

                    # Write sequence in chunks
                    for i in range(0, len(seq.sequence), line_width):
                        f.write(seq.sequence[i:i+line_width] + '\n')

            print(f"Wrote {len(sequences)} sequences to '{filename}'", file=sys.stderr)

        except IOError as e:
            print(f"Error writing file: {e}", file=sys.stderr)
            sys.exit(1)


class FASTQHandler:
    """
    Read and write FASTQ format files.

    FASTQ includes quality scores for each base.
    Quality encoding can be Phred33 (Illumina, most common) or Phred64.
    """

    @staticmethod
    def read(filename: str) -> Iterator[Sequence]:
        """
        Parse a FASTQ file.

        FASTQ format:
        Line 1: @header
        Line 2: sequence
        Line 3: + (separator)
        Line 4: quality scores

        Args:
            filename (str): Path to FASTQ file

        Yields:
            Sequence: Sequence objects with quality scores
        """
        try:
            with open(filename, 'r') as f:
                lines = []
                for line in f:
                    lines.append(line.rstrip('\n'))

                    # Process every 4 lines (one FASTQ record)
                    if len(lines) == 4:
                        header = lines[0][1:]  # Remove @ prefix
                        sequence = lines[1]
                        quality = lines[3]

                        yield Sequence(header, sequence, quality)
                        lines = []

                # Check for incomplete record
                if lines:
                    print("Warning: Incomplete FASTQ record at end of file", file=sys.stderr)

        except FileNotFoundError:
            print(f"Error: File '{filename}' not found.", file=sys.stderr)
            sys.exit(1)

    @staticmethod
    def write(filename: str, sequences: List[Sequence]):
        """
        Write sequences to FASTQ file.

        Args:
            filename (str): Output file path
            sequences (List[Sequence]): List of Sequence objects with quality scores
        """
        try:
            with open(filename, 'w') as f:
                for seq in sequences:
                    if seq.quality is None:
                        raise ValueError(f"Sequence '{seq.header}' lacks quality scores for FASTQ")

                    f.write(f"@{seq.header}\n")
                    f.write(f"{seq.sequence}\n")
                    f.write(f"+\n")
                    f.write(f"{seq.quality}\n")

            print(f"Wrote {len(sequences)} sequences to '{filename}'", file=sys.stderr)

        except IOError as e:
            print(f"Error writing file: {e}", file=sys.stderr)
            sys.exit(1)

    @staticmethod
    def phred_to_quality(phred_score: int, encoding: str = 'Phred33') -> str:
        """
        Convert numeric Phred score to quality character.

        Phred score represents probability of error:
        Q = -10 * log10(error_probability)

        Phred33: ASCII 33-126 (! to ~)
        Phred64: ASCII 64-126 (@ to ~)

        Args:
            phred_score (int): Quality score (0-40)
            encoding (str): 'Phred33' or 'Phred64'

        Returns:
            str: Single quality character
        """
        if encoding == 'Phred33':
            offset = 33
        elif encoding == 'Phred64':
            offset = 64
        else:
            raise ValueError("Encoding must be 'Phred33' or 'Phred64'")

        return chr(phred_score + offset)

    @staticmethod
    def quality_to_phred(quality_char: str, encoding: str = 'Phred33') -> int:
        """
        Convert quality character to numeric Phred score.

        Args:
            quality_char (str): Quality character
            encoding (str): 'Phred33' or 'Phred64'

        Returns:
            int: Phred score (0-40)
        """
        if encoding == 'Phred33':
            offset = 33
        elif encoding == 'Phred64':
            offset = 64
        else:
            raise ValueError("Encoding must be 'Phred33' or 'Phred64'")

        return ord(quality_char) - offset


class PHYLIPHandler:
    """
    Read and write PHYLIP format files.

    PHYLIP format is used for phylogenetic analysis.
    Sequences must be aligned (same length).
    """

    @staticmethod
    def read(filename: str) -> Tuple[int, int, List[Sequence]]:
        """
        Parse a PHYLIP file.

        First line: num_species num_sites
        Following lines: species_name (10 chars) sequence

        Args:
            filename (str): Path to PHYLIP file

        Returns:
            Tuple: (num_species, num_sites, sequences)
        """
        try:
            sequences = []

            with open(filename, 'r') as f:
                # Read header line
                header = f.readline().rstrip('\n').split()
                num_species = int(header[0])
                num_sites = int(header[1])

                # Read sequences
                for line in f:
                    if not line.strip():
                        continue

                    # PHYLIP format: first 10 chars are name, rest is sequence
                    parts = line.split()
                    if len(parts) >= 2:
                        name = parts[0]
                        # Join remaining parts as sequence (may be split across columns)
                        sequence = ''.join(parts[1:])

                        sequences.append(Sequence(name, sequence))

            return num_species, num_sites, sequences

        except FileNotFoundError:
            print(f"Error: File '{filename}' not found.", file=sys.stderr)
            sys.exit(1)

    @staticmethod
    def write(filename: str, sequences: List[Sequence]):
        """
        Write sequences to PHYLIP format.

        Args:
            filename (str): Output file path
            sequences (List[Sequence]): List of Sequence objects
        """
        try:
            if not sequences:
                raise ValueError("No sequences to write")

            # Verify all sequences have same length
            seq_length = len(sequences[0].sequence)
            for seq in sequences:
                if len(seq.sequence) != seq_length:
                    raise ValueError("All sequences must have the same length for PHYLIP format")

            with open(filename, 'w') as f:
                # Write header
                f.write(f"{len(sequences)} {seq_length}\n")

                # Write sequences
                for seq in sequences:
                    # Limit name to 10 characters
                    name = seq.header[:10]
                    # Pad name to 10 characters
                    name = name.ljust(10)
                    f.write(f"{name} {seq.sequence}\n")

            print(f"Wrote {len(sequences)} aligned sequences to '{filename}'", file=sys.stderr)

        except IOError as e:
            print(f"Error writing file: {e}", file=sys.stderr)
            sys.exit(1)


class FormatConverter:
    """Convert between different sequence formats."""

    @staticmethod
    def fastq_to_fasta(input_file: str, output_file: str):
        """
        Convert FASTQ to FASTA (discarding quality scores).

        Args:
            input_file (str): Input FASTQ file
            output_file (str): Output FASTA file
        """
        sequences = list(FASTQHandler.read(input_file))
        FASTAHandler.write(output_file, sequences)
        print(f"Converted {len(sequences)} sequences from FASTQ to FASTA", file=sys.stderr)

    @staticmethod
    def fasta_to_fastq(input_file: str, output_file: str, quality_score: int = 30):
        """
        Convert FASTA to FASTQ (generating quality scores).

        Args:
            input_file (str): Input FASTA file
            output_file (str): Output FASTQ file
            quality_score (int): Default quality score for all bases
        """
        sequences = list(FASTAHandler.read(input_file))

        # Generate quality scores
        quality_char = FASTQHandler.phred_to_quality(quality_score)

        for seq in sequences:
            seq.quality = quality_char * len(seq.sequence)

        FASTQHandler.write(output_file, sequences)
        print(f"Converted {len(sequences)} sequences from FASTA to FASTQ", file=sys.stderr)

    @staticmethod
    def fasta_to_phylip(input_file: str, output_file: str):
        """
        Convert FASTA to PHYLIP (requires aligned sequences).

        Args:
            input_file (str): Input FASTA file
            output_file (str): Output PHYLIP file
        """
        sequences = list(FASTAHandler.read(input_file))

        # Check if sequences are aligned
        if sequences:
            seq_length = len(sequences[0].sequence)
            for seq in sequences:
                if len(seq.sequence) != seq_length:
                    print("Warning: Sequences are not aligned! PHYLIP format requires aligned sequences.",
                          file=sys.stderr)
                    return

        PHYLIPHandler.write(output_file, sequences)
        print(f"Converted {len(sequences)} sequences to PHYLIP format", file=sys.stderr)

    @staticmethod
    def phylip_to_fasta(input_file: str, output_file: str):
        """
        Convert PHYLIP to FASTA.

        Args:
            input_file (str): Input PHYLIP file
            output_file (str): Output FASTA file
        """
        num_species, num_sites, sequences = PHYLIPHandler.read(input_file)
        FASTAHandler.write(output_file, sequences)
        print(f"Converted {num_species} sequences from PHYLIP to FASTA", file=sys.stderr)

    @staticmethod
    def detect_format(filename: str) -> str:
        """
        Attempt to detect sequence format from file content.

        Args:
            filename (str): Path to sequence file

        Returns:
            str: Detected format ('FASTA', 'FASTQ', or 'UNKNOWN')
        """
        try:
            with open(filename, 'r') as f:
                first_char = f.read(1)

                if first_char == '>':
                    return 'FASTA'
                elif first_char == '@':
                    return 'FASTQ'
                else:
                    return 'UNKNOWN'

        except FileNotFoundError:
            print(f"Error: File '{filename}' not found.", file=sys.stderr)
            return 'UNKNOWN'


def auto_convert(input_file: str, output_file: str, output_format: str = None):
    """
    Automatically convert between formats.

    If output format is not specified, it's inferred from output file extension.

    Args:
        input_file (str): Input sequence file
        output_file (str): Output file path
        output_format (str): Output format ('FASTA', 'FASTQ', 'PHYLIP')
                            If None, inferred from file extension

    Example:
        >>> auto_convert('seqs.fastq', 'seqs.fasta')  # Auto-detects FASTA format
        >>> auto_convert('seqs.fasta', 'seqs.phy', 'PHYLIP')
    """
    # Determine output format from extension if not specified
    if output_format is None:
        ext = Path(output_file).suffix.lower()
        format_map = {
            '.fasta': 'FASTA',
            '.fa': 'FASTA',
            '.fastq': 'FASTQ',
            '.fq': 'FASTQ',
            '.phy': 'PHYLIP',
            '.phylip': 'PHYLIP',
        }
        output_format = format_map.get(ext, 'UNKNOWN')

    if output_format == 'UNKNOWN':
        print(f"Error: Could not determine output format", file=sys.stderr)
        sys.exit(1)

    # Detect input format
    input_format = FormatConverter.detect_format(input_file)

    if input_format == 'UNKNOWN':
        print(f"Warning: Could not detect input format", file=sys.stderr)

    # Perform conversion
    converter = FormatConverter()

    if input_format == 'FASTQ' and output_format == 'FASTA':
        converter.fastq_to_fasta(input_file, output_file)
    elif input_format == 'FASTA' and output_format == 'FASTQ':
        converter.fasta_to_fastq(input_file, output_file)
    elif input_format == 'FASTA' and output_format == 'PHYLIP':
        converter.fasta_to_phylip(input_file, output_file)
    elif input_format == 'PHYLIP' and output_format == 'FASTA':
        converter.phylip_to_fasta(input_file, output_file)
    elif input_format == output_format:
        # Same format, just copy
        import shutil
        shutil.copy(input_file, output_file)
        print(f"Copied {input_file} to {output_file}", file=sys.stderr)
    else:
        print(f"Error: Conversion from {input_format} to {output_format} not implemented",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    print("Format Converter Module")
    print("This module converts between sequence formats.")
    print("\nSupported formats:")
    print("  - FASTA: Simple text format")
    print("  - FASTQ: FASTA with quality scores")
    print("  - PHYLIP: Phylogenetics format (aligned)")
    print("\nExample usage:")
    print("  from format_converter import auto_convert")
    print("  auto_convert('seqs.fastq', 'seqs.fasta')")
