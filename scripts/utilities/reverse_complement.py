#!/usr/bin/env python3
"""
Reverse Complement and Sequence Manipulation Tools

This module provides functions for manipulating DNA sequences, including:
- Reverse complement (used in DNA analysis to find the same sequence on opposite strand)
- Simple complement
- Reverse
- Sequence validation
- Palindrome detection

The reverse complement is essential for DNA barcoding because:
1. DNA is double-stranded (forward and reverse complement strands)
2. Sequences may be given from either strand
3. BLAST searches will find matches on both strands

Example:
    Original:     5'-ATCGATCG-3'
    Complement:   3'-TAGCTAGC-5'
    Rev Complement: 5'-CGATCGAT-3'

Author: DNA Barcoding Analysis Course
License: MIT
"""

import sys
from typing import Tuple, List


# Define complement mappings
COMPLEMENT_DNA = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
    'N': 'N',  # Ambiguous base
    'R': 'Y',  # Purine -> Pyrimidine
    'Y': 'R',  # Pyrimidine -> Purine
    'W': 'W',  # Weak (A or T)
    'S': 'S',  # Strong (G or C)
    'K': 'M',  # Keto -> Amino
    'M': 'K',  # Amino -> Keto
    'B': 'V',  # Not A -> Not T
    'D': 'H',  # Not C -> Not G
    'H': 'D',  # Not G -> Not C
    'V': 'B',  # Not T -> Not A
}

COMPLEMENT_RNA = {
    'A': 'U',
    'U': 'A',
    'G': 'C',
    'C': 'G',
    'N': 'N',
    'R': 'Y',
    'Y': 'R',
    'W': 'W',
    'S': 'S',
    'K': 'M',
    'M': 'K',
    'B': 'V',
    'D': 'H',
    'H': 'D',
    'V': 'B',
}


def complement(sequence: str, molecule_type: str = 'DNA') -> str:
    """
    Return the complement of a sequence.

    In DNA:
    - A pairs with T
    - T pairs with A
    - G pairs with C
    - C pairs with G

    In RNA:
    - A pairs with U
    - U pairs with A
    - G pairs with C
    - C pairs with G

    Args:
        sequence (str): Input sequence
        molecule_type (str): 'DNA' or 'RNA' (default: 'DNA')

    Returns:
        str: Complement sequence

    Example:
        >>> comp = complement("ATCGATCG")
        >>> print(comp)
        TAGCTAGC
    """
    sequence = sequence.upper()
    complement_map = COMPLEMENT_DNA if molecule_type.upper() == 'DNA' else COMPLEMENT_RNA

    complemented = []
    for base in sequence:
        if base in complement_map:
            complemented.append(complement_map[base])
        else:
            # Unknown base, keep as N
            complemented.append('N')

    return ''.join(complemented)


def reverse(sequence: str) -> str:
    """
    Return the reverse of a sequence.

    Note: This is NOT the same as complement or reverse complement.

    Args:
        sequence (str): Input sequence

    Returns:
        str: Reversed sequence

    Example:
        >>> rev = reverse("ATCGATCG")
        >>> print(rev)
        GCTAGCAT
    """
    return sequence[::-1]


def reverse_complement(sequence: str, molecule_type: str = 'DNA') -> str:
    """
    Return the reverse complement of a sequence.

    This is the same sequence as read from the opposite strand in the 5' to 3' direction.
    Essential for DNA analysis since sequences can be found on either strand.

    The process:
    1. Get complement of the sequence
    2. Reverse it

    Or equivalently:
    1. Reverse the sequence
    2. Get complement of the result

    Args:
        sequence (str): Input sequence
        molecule_type (str): 'DNA' or 'RNA' (default: 'DNA')

    Returns:
        str: Reverse complement sequence

    Example:
        >>> rc = reverse_complement("ATCGATCG")
        >>> print(rc)
        CGATCGAT

        >>> # Round trip - reverse complement twice equals original
        >>> original = "ATCGATCG"
        >>> assert reverse_complement(reverse_complement(original)) == original
    """
    # Method 1: complement then reverse
    complemented = complement(sequence, molecule_type)
    return reverse(complemented)

    # Method 2 (equivalent): reverse then complement
    # reversed_seq = reverse(sequence)
    # return complement(reversed_seq, molecule_type)


def reverse_complement_pair(forward: str, reverse_seq: str,
                          molecule_type: str = 'DNA') -> Tuple[str, str]:
    """
    Verify that two sequences are reverse complements of each other.

    Args:
        forward (str): First sequence
        reverse_seq (str): Second sequence
        molecule_type (str): 'DNA' or 'RNA'

    Returns:
        Tuple[str, str]: (forward, reverse_complement of forward)
                        Useful for validating sequence pairs

    Example:
        >>> forward = "ATCGATCG"
        >>> reverse = "CGATCGAT"
        >>> f, r = reverse_complement_pair(forward, reverse)
        >>> print(f"Forward: {f}")
        >>> print(f"Expected RC: {r}")
    """
    rc = reverse_complement(forward, molecule_type)
    return forward, rc


def is_palindrome(sequence: str) -> bool:
    """
    Check if a sequence is a palindrome (same as its reverse complement).

    Palindromic sequences have biological significance (recognition sites for enzymes).

    Args:
        sequence (str): Input sequence

    Returns:
        bool: True if sequence equals its reverse complement

    Example:
        >>> is_palindrome("GAATTC")  # EcoRI restriction site
        True
        >>> is_palindrome("ATCGATCG")
        False
    """
    sequence = sequence.upper()
    rc = reverse_complement(sequence)
    return sequence == rc


def find_palindromes(sequence: str, min_length: int = 4) -> List[Tuple[int, str]]:
    """
    Find all palindromic subsequences of a given minimum length.

    Palindromes are important in molecular biology as they're often
    restriction enzyme recognition sites.

    Args:
        sequence (str): Input sequence
        min_length (int): Minimum length of palindrome (default: 4)

    Returns:
        List[Tuple[int, str]]: List of (position, palindrome) tuples

    Example:
        >>> pals = find_palindromes("ATGAATTCATCG")
        >>> for pos, pal in pals:
        ...     print(f"Position {pos}: {pal}")
        Position 3: GAATTC
    """
    sequence = sequence.upper()
    palindromes = []

    # Search for palindromes of increasing length
    for length in range(min_length, len(sequence) + 1):
        for i in range(len(sequence) - length + 1):
            subseq = sequence[i:i+length]
            if is_palindrome(subseq):
                palindromes.append((i, subseq))

    return palindromes


def validate_sequence(sequence: str, molecule_type: str = 'DNA') -> Tuple[bool, str]:
    """
    Validate that a sequence contains only valid bases.

    Args:
        sequence (str): Input sequence
        molecule_type (str): 'DNA' or 'RNA'

    Returns:
        Tuple[bool, str]: (is_valid, error_message)

    Example:
        >>> is_valid, msg = validate_sequence("ATCGATCG")
        >>> print(f"Valid: {is_valid}")
        Valid: True

        >>> is_valid, msg = validate_sequence("ATCXATCG")
        >>> print(f"Valid: {is_valid}, Message: {msg}")
        Valid: False, Message: Invalid character 'X' at position 3
    """
    sequence = sequence.upper()

    if molecule_type.upper() == 'DNA':
        valid_chars = set('ATCGNRYWSKMBDHV')
        mol_name = 'DNA'
    else:
        valid_chars = set('AUCGNRYWSKMBDHV')
        mol_name = 'RNA'

    for i, base in enumerate(sequence):
        if base not in valid_chars:
            return False, f"Invalid {mol_name} character '{base}' at position {i+1}"

    return True, "Valid sequence"


def rotate_sequence(sequence: str, positions: int = 1) -> str:
    """
    Rotate a sequence by moving characters from the start to the end.

    Useful for finding different reading frames or rotated sequences.

    Args:
        sequence (str): Input sequence
        positions (int): Number of positions to rotate (default: 1)

    Returns:
        str: Rotated sequence

    Example:
        >>> seq = "ATCGATCG"
        >>> print(rotate_sequence(seq, 1))
        TCGATCGA
        >>> print(rotate_sequence(seq, 3))
        GATCGATC
    """
    if len(sequence) == 0:
        return sequence

    # Normalize positions to sequence length
    positions = positions % len(sequence)
    return sequence[positions:] + sequence[:positions]


def generate_all_orientations(sequence: str) -> dict:
    """
    Generate all 4 possible orientations of a DNA sequence.

    Useful for comprehensive sequence analysis:
    1. Forward (original)
    2. Reverse (backwards)
    3. Complement (opposite strand, same direction)
    4. Reverse Complement (opposite strand, backwards)

    Args:
        sequence (str): Input sequence

    Returns:
        dict: Dictionary with keys: 'forward', 'reverse', 'complement', 'reverse_complement'

    Example:
        >>> seqs = generate_all_orientations("ATCG")
        >>> for orientation, seq in seqs.items():
        ...     print(f"{orientation}: {seq}")
        forward: ATCG
        reverse: GCTA
        complement: TAGC
        reverse_complement: CGAT
    """
    return {
        'forward': sequence.upper(),
        'reverse': reverse(sequence.upper()),
        'complement': complement(sequence),
        'reverse_complement': reverse_complement(sequence),
    }


def process_fasta_reverse_complement(input_file: str, output_file: str):
    """
    Read a FASTA file and write the reverse complement of each sequence.

    Args:
        input_file (str): Input FASTA file
        output_file (str): Output FASTA file

    Example:
        >>> process_fasta_reverse_complement('seqs.fasta', 'seqs_rc.fasta')
    """
    from fasta_tools import read_fasta, write_fasta

    try:
        sequences = []
        for header, sequence in read_fasta(input_file):
            rc_seq = reverse_complement(sequence)
            # Add "_RC" to header to indicate reverse complement
            new_header = f"{header}_RC"
            sequences.append((new_header, rc_seq))

        count = write_fasta(output_file, sequences)
        print(f"Generated reverse complements for {count} sequences", file=sys.stderr)

    except Exception as e:
        print(f"Error processing FASTA file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    # Example usage
    print("Reverse Complement Module")
    print("This module provides DNA sequence manipulation tools.")
    print("\nExample usage:")
    print("  from reverse_complement import reverse_complement")
    print("  rc = reverse_complement('ATCGATCG')")
    print("  print(f'Reverse complement: {rc}')")
    print("\nKey concepts:")
    print("  - Complement: A<->T, G<->C")
    print("  - Reverse: Read backwards")
    print("  - Reverse complement: Both complement AND reverse")
    print("\nWhy it matters:")
    print("  - DNA is double-stranded")
    print("  - Same sequence appears on both strands")
    print("  - Reverse complement is the sequence from the other strand")
