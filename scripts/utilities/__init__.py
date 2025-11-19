"""
DNA Barcoding Analysis - Bioinformatics Utility Scripts

This package provides general-purpose bioinformatics tools for sequence analysis.

Available modules:
    - fasta_tools: Parse, split, merge FASTA files
    - sequence_stats: Calculate GC content, length, composition
    - reverse_complement: Reverse complement sequences and related operations
    - format_converter: Convert between FASTA/FASTQ/GenBank/PHYLIP formats
    - subsample_sequences: Random subsampling and sequence filtering

Example:
    >>> from utilities.fasta_tools import read_fasta
    >>> for header, sequence in read_fasta('data.fasta'):
    ...     print(f"{header}: {len(sequence)} bp")

Author: DNA Barcoding Analysis Course
License: MIT
"""

__version__ = "1.0.0"

# Import main classes and functions for convenience
try:
    from .fasta_tools import (
        read_fasta,
        write_fasta,
        split_fasta,
        merge_fasta,
        filter_fasta,
        count_sequences,
        get_sequence_by_id,
    )
except ImportError:
    pass

try:
    from .sequence_stats import SequenceStats, analyze_sequence_file
except ImportError:
    pass

try:
    from .reverse_complement import (
        reverse_complement,
        complement,
        reverse,
        is_palindrome,
        find_palindromes,
    )
except ImportError:
    pass

try:
    from .format_converter import (
        FormatConverter,
        FASTAHandler,
        FASTQHandler,
        PHYLIPHandler,
        auto_convert,
    )
except ImportError:
    pass

try:
    from .subsample_sequences import (
        SequenceSubsampler,
        subsample_fasta,
        subsample_fasta_stratified,
    )
except ImportError:
    pass

__all__ = [
    # fasta_tools
    'read_fasta',
    'write_fasta',
    'split_fasta',
    'merge_fasta',
    'filter_fasta',
    'count_sequences',
    'get_sequence_by_id',
    # sequence_stats
    'SequenceStats',
    'analyze_sequence_file',
    # reverse_complement
    'reverse_complement',
    'complement',
    'reverse',
    'is_palindrome',
    'find_palindromes',
    # format_converter
    'FormatConverter',
    'FASTAHandler',
    'FASTQHandler',
    'PHYLIPHandler',
    'auto_convert',
    # subsample_sequences
    'SequenceSubsampler',
    'subsample_fasta',
    'subsample_fasta_stratified',
]
