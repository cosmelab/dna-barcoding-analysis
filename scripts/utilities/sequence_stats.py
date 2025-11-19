#!/usr/bin/env python3
"""
Sequence Statistics and Composition Analysis

This module calculates various statistics about DNA/RNA sequences including:
- GC content (percentage of G and C nucleotides)
- Length statistics
- Nucleotide composition
- Codon usage (for coding sequences)
- Ambiguous nucleotide detection

Author: DNA Barcoding Analysis Course
License: MIT
"""

import sys
from pathlib import Path
from typing import Dict, Tuple, List
from collections import Counter
import statistics


class SequenceStats:
    """
    Calculate statistics for DNA/RNA sequences.

    DNA/RNA bases:
    - A (Adenine)
    - T (Thymine, in DNA) / U (Uracil, in RNA)
    - G (Guanine)
    - C (Cytosine)
    - N (Ambiguous/Unknown)
    """

    # Standard IUPAC nucleotide codes
    DNA_BASES = set('ATCGN')
    RNA_BASES = set('AUCGN')
    IUPAC_CODES = {
        'A': 'Adenine',
        'T': 'Thymine (DNA)',
        'U': 'Uracil (RNA)',
        'G': 'Guanine',
        'C': 'Cytosine',
        'R': 'Purine (A or G)',
        'Y': 'Pyrimidine (C or T)',
        'W': 'Weak (A or T)',
        'S': 'Strong (G or C)',
        'K': 'Keto (G or T)',
        'M': 'Amino (A or C)',
        'B': 'Not A (C, G, or T)',
        'D': 'Not C (A, G, or T)',
        'H': 'Not G (A, C, or T)',
        'V': 'Not T (A, C, or G)',
        'N': 'Any nucleotide',
    }

    def __init__(self, sequence: str):
        """
        Initialize with a sequence.

        Args:
            sequence (str): DNA/RNA sequence (can be upper or lower case)
        """
        self.original_sequence = sequence
        self.sequence = sequence.upper()
        self.length = len(self.sequence)

    def gc_content(self) -> float:
        """
        Calculate GC content as percentage.

        GC content is important because:
        - It affects DNA stability (G-C bonds are stronger)
        - It influences sequence complexity
        - It's used in phylogenetic analysis
        - It can indicate mutation patterns

        Returns:
            float: GC content as percentage (0-100)

        Example:
            >>> stats = SequenceStats("ATCGATCGATCGATCG")
            >>> print(f"GC content: {stats.gc_content():.1f}%")
            GC content: 50.0%
        """
        if self.length == 0:
            return 0.0

        gc_count = self.sequence.count('G') + self.sequence.count('C')
        return (gc_count / self.length) * 100

    def nucleotide_composition(self) -> Dict[str, int]:
        """
        Count occurrences of each nucleotide.

        Returns:
            Dict[str, int]: Dictionary with nucleotide counts

        Example:
            >>> stats = SequenceStats("ATCGATCG")
            >>> composition = stats.nucleotide_composition()
            >>> print(composition)
            {'A': 2, 'T': 2, 'C': 2, 'G': 2}
        """
        composition = {}
        counter = Counter(self.sequence)

        # Include all standard bases even if count is 0
        for base in 'ATCGN':
            composition[base] = counter.get(base, 0)

        # Include ambiguous codes if present
        for base in counter:
            if base not in composition:
                composition[base] = counter[base]

        return dict(sorted(composition.items()))

    def nucleotide_percentages(self) -> Dict[str, float]:
        """
        Calculate percentage composition of each nucleotide.

        Returns:
            Dict[str, float]: Dictionary with nucleotide percentages

        Example:
            >>> stats = SequenceStats("ATCGATCGATCG")
            >>> percentages = stats.nucleotide_percentages()
            >>> print(f"A: {percentages['A']:.1f}%")
            A: 33.3%
        """
        if self.length == 0:
            return {}

        composition = self.nucleotide_composition()
        return {
            base: (count / self.length) * 100
            for base, count in composition.items()
        }

    def at_content(self) -> float:
        """
        Calculate AT content as percentage.

        AT content is complementary to GC content.

        Returns:
            float: AT content as percentage (0-100)

        Example:
            >>> stats = SequenceStats("ATCGATCG")
            >>> print(f"AT content: {stats.at_content():.1f}%")
            AT content: 50.0%
        """
        if self.length == 0:
            return 0.0

        at_count = self.sequence.count('A') + self.sequence.count('T')
        return (at_count / self.length) * 100

    def ambiguous_nucleotides(self) -> List[str]:
        """
        Identify any non-standard nucleotides in the sequence.

        Returns:
            List[str]: List of ambiguous nucleotide positions and codes

        Example:
            >>> stats = SequenceStats("ATCGNATCGN")
            >>> ambiguous = stats.ambiguous_nucleotides()
            >>> print(ambiguous)
            ['N at position 4', 'N at position 9']
        """
        ambiguous = []
        standard_bases = set('ATCG')

        for i, base in enumerate(self.sequence):
            if base not in standard_bases:
                ambiguous.append(f"{base} at position {i + 1}")

        return ambiguous

    def codon_usage(self) -> Dict[str, int]:
        """
        Calculate codon usage for a coding sequence.

        Codons are 3-nucleotide sequences that code for amino acids.
        This function extracts all codons from the sequence.

        Returns:
            Dict[str, int]: Dictionary with codon frequencies

        Example:
            >>> stats = SequenceStats("ATGATGATG")
            >>> codons = stats.codon_usage()
            >>> print(codons)
            {'ATG': 3}
        """
        codon_dict = {}

        # Extract codons (non-overlapping, reading frame 0)
        for i in range(0, len(self.sequence) - 2, 3):
            codon = self.sequence[i:i+3]

            # Only count complete codons
            if len(codon) == 3:
                codon_dict[codon] = codon_dict.get(codon, 0) + 1

        return dict(sorted(codon_dict.items()))

    def is_valid_dna(self) -> bool:
        """
        Check if sequence contains only valid DNA characters.

        Valid characters: A, T, C, G, N, and IUPAC ambiguity codes

        Returns:
            bool: True if valid DNA sequence

        Example:
            >>> stats = SequenceStats("ATCG")
            >>> print(stats.is_valid_dna())
            True
        """
        valid_chars = set('ATCGNRYWSKMBDHV')
        return all(base in valid_chars for base in self.sequence)

    def is_valid_rna(self) -> bool:
        """
        Check if sequence contains only valid RNA characters.

        Valid characters: A, U, C, G, N, and IUPAC ambiguity codes

        Returns:
            bool: True if valid RNA sequence
        """
        valid_chars = set('AUCGNRYWSKMBDHV')
        return all(base in valid_chars for base in self.sequence)

    def summary(self) -> Dict[str, any]:
        """
        Generate a complete summary of sequence statistics.

        Returns:
            Dict: Dictionary containing all major statistics

        Example:
            >>> stats = SequenceStats("ATCGATCGATCG")
            >>> summary = stats.summary()
            >>> for key, value in summary.items():
            ...     print(f"{key}: {value}")
        """
        composition = self.nucleotide_composition()
        percentages = self.nucleotide_percentages()

        return {
            'length': self.length,
            'gc_content': round(self.gc_content(), 2),
            'at_content': round(self.at_content(), 2),
            'composition': composition,
            'percentages': {k: round(v, 2) for k, v in percentages.items()},
            'is_valid_dna': self.is_valid_dna(),
            'is_valid_rna': self.is_valid_rna(),
            'ambiguous_nucleotides': self.ambiguous_nucleotides(),
            'num_codons': self.length // 3,
        }

    def print_summary(self):
        """
        Print a nicely formatted summary of sequence statistics.
        """
        summary = self.summary()

        print(f"\n{'='*50}")
        print(f"Sequence Statistics Summary")
        print(f"{'='*50}")
        print(f"Sequence Length: {summary['length']} bp")
        print(f"GC Content: {summary['gc_content']:.1f}%")
        print(f"AT Content: {summary['at_content']:.1f}%")
        print(f"\nNucleotide Composition:")
        for base, count in summary['composition'].items():
            if count > 0:
                pct = summary['percentages'][base]
                print(f"  {base}: {count:>5} ({pct:>6.2f}%)")

        if summary['ambiguous_nucleotides']:
            print(f"\nAmbiguous Nucleotides: {len(summary['ambiguous_nucleotides'])}")
            for amb in summary['ambiguous_nucleotides'][:5]:  # Show first 5
                print(f"  - {amb}")
            if len(summary['ambiguous_nucleotides']) > 5:
                print(f"  ... and {len(summary['ambiguous_nucleotides']) - 5} more")

        print(f"\nSequence Type:")
        print(f"  Valid DNA: {summary['is_valid_dna']}")
        print(f"  Valid RNA: {summary['is_valid_rna']}")
        print(f"{'='*50}\n")


def analyze_sequence_file(filename: str, seq_id_filter: str = None):
    """
    Analyze all sequences in a FASTA file.

    Args:
        filename (str): Path to FASTA file
        seq_id_filter (str): Optional - only analyze sequences matching this ID

    Example:
        >>> analyze_sequence_file('sequences.fasta')
    """
    from fasta_tools import read_fasta

    try:
        count = 0
        total_length = 0
        gc_contents = []

        for header, sequence in read_fasta(filename):
            # Skip if ID filter is set and doesn't match
            if seq_id_filter and seq_id_filter not in header:
                continue

            count += 1
            stats = SequenceStats(sequence)
            total_length += len(sequence)
            gc_contents.append(stats.gc_content())

            print(f"\n{header}")
            print(f"  Length: {len(sequence)} bp")
            print(f"  GC Content: {stats.gc_content():.1f}%")
            composition = stats.nucleotide_composition()
            print(f"  Composition: A={composition['A']} T={composition['T']} "
                  f"G={composition['G']} C={composition['C']}")

        # Print summary statistics
        if count > 0:
            print(f"\n{'='*50}")
            print(f"Summary for {filename}")
            print(f"{'='*50}")
            print(f"Total sequences: {count}")
            print(f"Total length: {total_length} bp")
            print(f"Average length: {total_length/count:.0f} bp")
            print(f"Average GC content: {statistics.mean(gc_contents):.1f}%")
            if count > 1:
                print(f"GC content range: {min(gc_contents):.1f}% - {max(gc_contents):.1f}%")

    except Exception as e:
        print(f"Error analyzing file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    # Example usage
    print("Sequence Statistics Module")
    print("This module calculates statistics for DNA/RNA sequences.")
    print("\nExample usage:")
    print("  from sequence_stats import SequenceStats")
    print("  stats = SequenceStats('ATCGATCGATCG')")
    print("  print(f'GC content: {stats.gc_content():.1f}%')")
    print("  stats.print_summary()")
