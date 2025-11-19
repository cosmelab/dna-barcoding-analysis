#!/usr/bin/env python3
"""
Sequence Subsampling and Filtering Tools

This module provides functions for subsampling (randomly selecting) sequences
from large datasets. Useful for:
- Creating smaller test datasets
- Reducing computational burden for analysis
- Creating representative samples
- Cross-validation in machine learning

Methods supported:
- Random sampling: Select random sequences
- Stratified sampling: Sample proportionally from groups
- Length-based sampling: Sample by sequence length ranges
- Quality-based sampling: Sample by sequence quality

Author: DNA Barcoding Analysis Course
License: MIT
"""

import sys
import random
from pathlib import Path
from typing import List, Tuple, Dict, Set
import statistics


class SequenceSubsampler:
    """
    Tools for subsampling sequences from large datasets.
    """

    def __init__(self, seed: int = None):
        """
        Initialize the subsampler with optional random seed.

        Args:
            seed (int): Random seed for reproducibility (optional)

        Example:
            >>> subsampler = SequenceSubsampler(seed=42)
            # Results will be reproducible
        """
        if seed is not None:
            random.seed(seed)
        self.seed = seed

    @staticmethod
    def read_fasta_with_metadata(filename: str) -> List[Tuple[str, str, int]]:
        """
        Read FASTA file and return sequences with metadata.

        Args:
            filename (str): Path to FASTA file

        Returns:
            List[Tuple[str, str, int]]: List of (header, sequence, length) tuples

        Example:
            >>> seqs = SequenceSubsampler.read_fasta_with_metadata('data.fasta')
            >>> print(f"Read {len(seqs)} sequences")
        """
        from fasta_tools import read_fasta

        sequences = []
        try:
            for header, sequence in read_fasta(filename):
                sequences.append((header, sequence, len(sequence)))
            return sequences
        except Exception as e:
            print(f"Error reading FASTA file: {e}", file=sys.stderr)
            sys.exit(1)

    @staticmethod
    def write_fasta(filename: str, sequences: List[Tuple[str, str, int]]) -> int:
        """
        Write subsampled sequences to FASTA file.

        Args:
            filename (str): Output file path
            sequences (List): List of sequence tuples

        Returns:
            int: Number of sequences written
        """
        from fasta_tools import write_fasta

        # Convert back to (header, sequence) format
        seq_list = [(header, seq) for header, seq, _ in sequences]
        return write_fasta(filename, seq_list)

    @staticmethod
    def random_sample(sequences: List[Tuple[str, str, int]],
                     n: int = None, fraction: float = None) -> List[Tuple[str, str, int]]:
        """
        Randomly sample sequences.

        Either n (number) or fraction (percentage) must be specified.

        Args:
            sequences (List): List of sequence tuples
            n (int): Number of sequences to sample
            fraction (float): Fraction of sequences to sample (0.0-1.0)

        Returns:
            List: Subsampled sequences

        Example:
            >>> all_seqs = SequenceSubsampler.read_fasta_with_metadata('data.fasta')
            >>> sample = SequenceSubsampler.random_sample(all_seqs, n=100)
            >>> print(f"Sampled {len(sample)} sequences")

            >>> # Or sample 50% of sequences
            >>> sample = SequenceSubsampler.random_sample(all_seqs, fraction=0.5)
        """
        if n is None and fraction is None:
            raise ValueError("Either 'n' or 'fraction' must be specified")

        if n is not None and fraction is not None:
            raise ValueError("Specify either 'n' or 'fraction', not both")

        # Calculate sample size
        if fraction is not None:
            if fraction < 0 or fraction > 1:
                raise ValueError("Fraction must be between 0 and 1")
            n = int(len(sequences) * fraction)

        if n > len(sequences):
            print(f"Warning: Requested {n} sequences but only {len(sequences)} available",
                  file=sys.stderr)
            n = len(sequences)

        return random.sample(sequences, n)

    @staticmethod
    def length_based_sample(sequences: List[Tuple[str, str, int]],
                           min_length: int = None,
                           max_length: int = None,
                           n: int = None,
                           fraction: float = None) -> List[Tuple[str, str, int]]:
        """
        Sample sequences within a length range.

        Args:
            sequences (List): List of sequence tuples
            min_length (int): Minimum sequence length
            max_length (int): Maximum sequence length
            n (int): Number of sequences to sample
            fraction (float): Fraction of filtered sequences to sample

        Returns:
            List: Filtered and subsampled sequences

        Example:
            >>> seqs = SequenceSubsampler.read_fasta_with_metadata('data.fasta')
            >>> # Get 100 sequences between 400-600 bp
            >>> sample = SequenceSubsampler.length_based_sample(
            ...     seqs, min_length=400, max_length=600, n=100)
        """
        # Filter by length
        if min_length is not None:
            sequences = [s for s in sequences if s[2] >= min_length]

        if max_length is not None:
            sequences = [s for s in sequences if s[2] <= max_length]

        print(f"Found {len(sequences)} sequences in length range", file=sys.stderr)

        # Now subsample
        if n is not None or fraction is not None:
            return SequenceSubsampler.random_sample(sequences, n=n, fraction=fraction)
        else:
            return sequences

    @staticmethod
    def stratified_sample(sequences: List[Tuple[str, str, int]],
                         n: int = None,
                         fraction: float = None,
                         num_strata: int = 5) -> List[Tuple[str, str, int]]:
        """
        Stratified random sampling by sequence length.

        Stratification ensures representation across different sequence lengths.
        Useful for ensuring diversity in subsamples.

        Args:
            sequences (List): List of sequence tuples
            n (int): Total number of sequences to sample
            fraction (float): Fraction of total sequences to sample
            num_strata (int): Number of length strata (default: 5)

        Returns:
            List: Stratified subsample

        Example:
            >>> seqs = SequenceSubsampler.read_fasta_with_metadata('data.fasta')
            >>> # Get 200 sequences stratified by length
            >>> sample = SequenceSubsampler.stratified_sample(seqs, n=200, num_strata=5)
        """
        if n is None and fraction is None:
            raise ValueError("Either 'n' or 'fraction' must be specified")

        # Calculate total sample size
        if fraction is not None:
            n = int(len(sequences) * fraction)

        # Sort by length
        sorted_seqs = sorted(sequences, key=lambda x: x[2])

        # Divide into strata
        strata_size = len(sequences) // num_strata
        strata = []

        for i in range(num_strata):
            start = i * strata_size
            if i == num_strata - 1:
                # Last stratum gets remaining sequences
                end = len(sequences)
            else:
                end = (i + 1) * strata_size

            strata.append(sorted_seqs[start:end])

        # Sample proportionally from each stratum
        sample = []
        for stratum in strata:
            stratum_sample_size = max(1, int((len(stratum) / len(sequences)) * n))
            if stratum_sample_size <= len(stratum):
                sample.extend(random.sample(stratum, stratum_sample_size))
            else:
                sample.extend(stratum)

        # Adjust if we have too many or too few
        if len(sample) > n:
            sample = random.sample(sample, n)
        elif len(sample) < n:
            # Add remaining from any stratum
            all_sampled = set(id(s) for s in sample)
            remaining = [s for s in sequences if id(s) not in all_sampled]
            needed = n - len(sample)
            if remaining:
                sample.extend(random.sample(remaining, min(needed, len(remaining))))

        return sample[:n]

    @staticmethod
    def grouped_sample(sequences: List[Tuple[str, str, int]],
                      group_func,
                      n: int = None,
                      fraction: float = None) -> List[Tuple[str, str, int]]:
        """
        Sample sequences grouped by a custom function.

        Allows flexible grouping based on any attribute of the header.

        Args:
            sequences (List): List of sequence tuples
            group_func: Function that takes header and returns group name
            n (int): Number of sequences to sample
            fraction (float): Fraction to sample

        Returns:
            List: Grouped and subsampled sequences

        Example:
            >>> seqs = SequenceSubsampler.read_fasta_with_metadata('data.fasta')
            >>> # Group by first part of header (before _)
            >>> group_func = lambda h: h.split('_')[0]
            >>> sample = SequenceSubsampler.grouped_sample(
            ...     seqs, group_func, fraction=0.1)
        """
        # Group sequences
        groups = {}
        for seq_tuple in sequences:
            header = seq_tuple[0]
            group = group_func(header)

            if group not in groups:
                groups[group] = []
            groups[group].append(seq_tuple)

        print(f"Grouped into {len(groups)} groups", file=sys.stderr)
        for group, seqs in groups.items():
            print(f"  {group}: {len(seqs)} sequences", file=sys.stderr)

        # Sample from each group proportionally
        if fraction is not None:
            n_total = int(len(sequences) * fraction)
        else:
            n_total = n

        sample = []
        for group, seqs in groups.items():
            group_fraction = len(seqs) / len(sequences)
            group_n = max(1, int(n_total * group_fraction))
            group_sample = random.sample(seqs, min(group_n, len(seqs)))
            sample.extend(group_sample)

        return sample[:n_total]

    @staticmethod
    def remove_duplicates(sequences: List[Tuple[str, str, int]],
                         by_sequence: bool = True) -> Tuple[List[Tuple[str, str, int]], int]:
        """
        Remove duplicate sequences from list.

        Args:
            sequences (List): List of sequence tuples
            by_sequence (bool): If True, compare by sequence content.
                               If False, compare by header.

        Returns:
            Tuple[List, int]: (deduplicated sequences, number removed)

        Example:
            >>> seqs = SequenceSubsampler.read_fasta_with_metadata('data.fasta')
            >>> unique_seqs, num_dupes = SequenceSubsampler.remove_duplicates(seqs)
            >>> print(f"Removed {num_dupes} duplicates")
        """
        seen = set()
        unique = []
        duplicates = 0

        for seq_tuple in sequences:
            if by_sequence:
                # Compare by sequence content
                key = seq_tuple[1]  # The sequence
            else:
                # Compare by header
                key = seq_tuple[0]  # The header

            if key not in seen:
                seen.add(key)
                unique.append(seq_tuple)
            else:
                duplicates += 1

        return unique, duplicates

    @staticmethod
    def subsample_by_identity(sequences: List[Tuple[str, str, int]],
                             max_identity: float = 0.95) -> List[Tuple[str, str, int]]:
        """
        Subsample to reduce sequence identity (simple clustering).

        Keeps sequences that are less similar to already selected sequences.
        Uses simple Hamming distance as similarity measure.

        Args:
            sequences (List): List of sequence tuples
            max_identity (float): Maximum identity threshold (0.0-1.0)

        Returns:
            List: Reduced sequence set

        Example:
            >>> seqs = SequenceSubsampler.read_fasta_with_metadata('data.fasta')
            >>> # Get non-redundant sequences with <95% identity
            >>> nr_seqs = SequenceSubsampler.subsample_by_identity(
            ...     seqs, max_identity=0.95)
        """
        def hamming_distance(seq1: str, seq2: str) -> float:
            """Calculate Hamming distance as fraction of differing positions."""
            if len(seq1) != len(seq2):
                # For different lengths, pad shorter sequence
                min_len = min(len(seq1), len(seq2))
                seq1 = seq1[:min_len]
                seq2 = seq2[:min_len]

            if min_len == 0:
                return 1.0

            differences = sum(1 for a, b in zip(seq1, seq2) if a != b)
            return differences / min_len

        selected = []

        for new_seq in sequences:
            # Check similarity to already selected sequences
            is_similar = False

            for selected_seq in selected:
                distance = hamming_distance(new_seq[1], selected_seq[1])
                identity = 1 - distance

                if identity > max_identity:
                    is_similar = True
                    break

            if not is_similar:
                selected.append(new_seq)

        return selected


def subsample_fasta(input_file: str, output_file: str,
                   n: int = None, fraction: float = None,
                   seed: int = None):
    """
    Subsample sequences from a FASTA file.

    Args:
        input_file (str): Input FASTA file
        output_file (str): Output FASTA file
        n (int): Number of sequences to sample
        fraction (float): Fraction of sequences to sample
        seed (int): Random seed for reproducibility

    Example:
        >>> subsample_fasta('large.fasta', 'sample.fasta', n=100, seed=42)
    """
    subsampler = SequenceSubsampler(seed=seed)

    # Read sequences
    sequences = subsampler.read_fasta_with_metadata(input_file)
    print(f"Read {len(sequences)} sequences from '{input_file}'", file=sys.stderr)

    # Subsample
    sample = subsampler.random_sample(sequences, n=n, fraction=fraction)
    print(f"Selected {len(sample)} sequences", file=sys.stderr)

    # Write output
    subsampler.write_fasta(output_file, sample)


def subsample_fasta_stratified(input_file: str, output_file: str,
                              n: int = None, fraction: float = None,
                              seed: int = None):
    """
    Subsample sequences using stratification by length.

    Ensures representation across different sequence length ranges.

    Args:
        input_file (str): Input FASTA file
        output_file (str): Output FASTA file
        n (int): Number of sequences to sample
        fraction (float): Fraction of sequences to sample
        seed (int): Random seed for reproducibility

    Example:
        >>> subsample_fasta_stratified('large.fasta', 'sample.fasta', fraction=0.1)
    """
    subsampler = SequenceSubsampler(seed=seed)

    # Read sequences
    sequences = subsampler.read_fasta_with_metadata(input_file)
    print(f"Read {len(sequences)} sequences", file=sys.stderr)

    # Calculate length statistics
    lengths = [s[2] for s in sequences]
    print(f"Sequence length: {min(lengths)}-{max(lengths)} bp "
          f"(mean: {statistics.mean(lengths):.0f})", file=sys.stderr)

    # Stratified sampling
    sample = subsampler.stratified_sample(sequences, n=n, fraction=fraction)
    print(f"Selected {len(sample)} sequences (stratified)", file=sys.stderr)

    # Write output
    subsampler.write_fasta(output_file, sample)


def get_sequence_statistics(sequences: List[Tuple[str, str, int]]) -> Dict:
    """
    Calculate statistics about sequences.

    Args:
        sequences (List): List of sequence tuples

    Returns:
        Dict: Statistics dictionary

    Example:
        >>> seqs = SequenceSubsampler.read_fasta_with_metadata('data.fasta')
        >>> stats = get_sequence_statistics(seqs)
        >>> print(f"Total sequences: {stats['count']}")
    """
    lengths = [s[2] for s in sequences]

    return {
        'count': len(sequences),
        'total_length': sum(lengths),
        'min_length': min(lengths) if lengths else 0,
        'max_length': max(lengths) if lengths else 0,
        'mean_length': statistics.mean(lengths) if lengths else 0,
        'median_length': statistics.median(lengths) if lengths else 0,
        'stdev_length': statistics.stdev(lengths) if len(lengths) > 1 else 0,
    }


if __name__ == "__main__":
    print("Sequence Subsampling Module")
    print("This module provides tools for subsampling sequences.")
    print("\nAvailable methods:")
    print("  - random_sample: Select random sequences")
    print("  - stratified_sample: Sample stratified by length")
    print("  - grouped_sample: Sample grouped by custom function")
    print("  - length_based_sample: Sample within length range")
    print("  - remove_duplicates: Remove duplicate sequences")
    print("  - subsample_by_identity: Reduce sequence redundancy")
    print("\nExample usage:")
    print("  from subsample_sequences import SequenceSubsampler")
    print("  subsampler = SequenceSubsampler(seed=42)")
    print("  seqs = subsampler.read_fasta_with_metadata('data.fasta')")
    print("  sample = subsampler.random_sample(seqs, n=100)")
