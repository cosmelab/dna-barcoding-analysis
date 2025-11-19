#!/usr/bin/env python3
"""
Calculate Pairwise Genetic Distances from Phylogenetic Trees and Sequences

PURPOSE:
    Compute evolutionary distances between sequences from:
    1. Phylogenetic tree branch lengths
    2. Direct sequence alignment analysis

PHYLOGENETIC CONCEPTS:
    - GENETIC DISTANCE: How different are two sequences?
    - BRANCH LENGTH: The distance in evolutionary time
    - SUBSTITUTION RATE: Number of changes per site
    - PAIRWISE DISTANCE: Distance between each pair of sequences

WHY THIS MATTERS FOR DNA BARCODING:
    - Species differ by ~3% (COI gene)
    - Individuals of same species ~0.3% difference
    - Genetic distance helps threshold species separation
    - Critical for species identification

DISTANCE TYPES:
    1. TREE DISTANCES: From branch lengths in phylogenetic tree
       - Uses evolutionary model implicit in tree
       - Faster to calculate
       - Requires tree file

    2. SEQUENCE DISTANCES: Calculated from aligned sequences directly
       - Hamming distance: Simple count of differences
       - P-distance: Proportion of different sites
       - Jukes-Cantor: Corrects for multiple substitutions
       - Kimura 2-parameter: Accounts for transition/transversion bias

REQUIREMENTS:
    - Python 3.6+
    - Biopython
    - pandas (for nice output tables)
    - numpy

INSTALLATION:
    conda install biopython pandas numpy

USAGE:
    # Calculate from tree branch lengths
    python calculate_distances.py -t tree.treefile

    # Calculate from sequence alignment
    python calculate_distances.py -a alignment.fasta

    # Save to CSV file
    python calculate_distances.py -a alignment.fasta -o distances.csv

    # Show only distances > 0.01 (1% difference)
    python calculate_distances.py -a alignment.fasta --min-distance 0.01

EXAMPLES:
    # Tree-based distances (if tree already exists)
    python calculate_distances.py -t COI_iqtree.treefile

    # Sequence-based distances with visualization
    python calculate_distances.py -a aligned_coi.fasta -o coi_distances.csv

    # Find species threshold (usually ~3% for DNA barcoding)
    python calculate_distances.py -a aligned_coi.fasta --threshold-analysis

AUTHOR: Educational Materials for DNA Barcoding
VERSION: 1.0
"""

import sys
import os
import argparse
from pathlib import Path
import numpy as np

from Bio import Phylo, SeqIO, AlignIO
from collections import defaultdict

###############################################################################
# PHYLOGENETIC CONCEPTS
###############################################################################

DISTANCE_CONCEPTS = """
UNDERSTANDING GENETIC DISTANCES:

1. WHAT IS GENETIC DISTANCE?
   - Measure of how different two sequences are
   - Units: Usually proportions (0.0 to ~0.5)
   - 0.01 = 1% different (1 change per 100 sites)
   - 0.1 = 10% different (10 changes per 100 sites)

2. DISTANCE TYPES:

   a) Hamming Distance:
      - Simplest: count positions where sequences differ
      - Don't account for multiple substitutions
      - Good for closely related sequences

   b) P-distance (Proportion):
      - Hamming distance / total length
      - Same as percent difference
      - Example: 10 differences in 100 positions = 0.10

   c) Jukes-Cantor:
      - Corrects for multiple substitutions at same site
      - Assumes all changes equally likely
      - Works for moderately divergent sequences

   d) Kimura 2-parameter:
      - Distinguishes transitions (A↔G, C↔T)
      - From transversions (all other changes)
      - More accurate for real DNA sequences

3. DNA BARCODING THRESHOLDS:
   - Within species (same species): ~0.003 (0.3%)
   - Between species (different species): ~0.03 (3%)
   - This creates a "barcoding gap"
   - Species ID: if distance < 0.03 to reference = same species

4. WHY DOES DISTANCE MATTER?

   Example: Beetle identification
   - You have your unknown beetle sequence
   - Calculate distance to all reference sequences
   - If closest match is 0.01 = probably same species
   - If closest match is 0.05 = different species

5. TREE DISTANCES VS SEQUENCE DISTANCES:

   Tree distances:
   - Use branch lengths from tree
   - Incorporate evolutionary model
   - More biologically meaningful
   - Slower if many sequences

   Sequence distances:
   - Calculate directly from alignment
   - No model needed
   - Faster
   - Good for quick assessment

6. CORRECTING FOR MULTIPLE SUBSTITUTIONS:

   Problem: Sites can mutate multiple times
   - Simple count misses hidden changes
   - Two sequences at a site might both have mutated
   - Appear same but actually different

   Solution: Correction models
   - Account for unobserved changes
   - Different models for different data
   - Mathematical correction to raw distance
"""

###############################################################################
# FUNCTION: Validate input files
###############################################################################

def validate_file(filepath):
    """Check if file exists and is readable."""
    if not os.path.exists(filepath):
        print(f"ERROR: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    if os.path.getsize(filepath) == 0:
        print(f"ERROR: File is empty: {filepath}", file=sys.stderr)
        sys.exit(1)

    return True


###############################################################################
# FUNCTION: Calculate sequence-based distances
###############################################################################

def calculate_p_distance(seq1, seq2):
    """
    Calculate proportion distance (p-distance) between two sequences.

    P-distance = number of differences / sequence length

    Args:
        seq1, seq2: Bio.Seq objects

    Returns:
        Float: proportion difference (0.0 to 1.0)

    Explanation:
        - Simplest distance metric
        - What proportion of sites differ?
        - Example: 10 differences in 100 sites = 0.10
    """

    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be same length")

    differences = sum(s1 != s2 for s1, s2 in zip(seq1, seq2))
    return differences / len(seq1)


def calculate_jukes_cantor(seq1, seq2):
    """
    Calculate Jukes-Cantor corrected distance.

    Corrects p-distance for multiple substitutions.
    Assumes: All substitution rates equal, 75% unchanged = same as start

    Formula: d = -3/4 * ln(1 - 4/3 * p)
    where p = p-distance

    Args:
        seq1, seq2: Bio.Seq objects

    Returns:
        Float: Corrected distance

    Explanation:
        - Raw difference count can underestimate divergence
        - Multiple hits at same site look like no change
        - Jukes-Cantor corrects for this
        - More accurate for moderately divergent sequences
    """

    p = calculate_p_distance(seq1, seq2)

    # Avoid domain error for ln
    if p >= 0.75:
        return float('inf')  # Divergence too high for this model

    # Jukes-Cantor formula
    try:
        jc = -0.75 * np.log(1 - (4/3) * p)
        return max(0, jc)  # Distance can't be negative
    except:
        return float('inf')


def calculate_kimura_2p(seq1, seq2):
    """
    Calculate Kimura 2-parameter distance.

    Distinguishes transitions (A↔G, C↔T) from transversions.
    More accurate than Jukes-Cantor for real DNA data.

    Args:
        seq1, seq2: Bio.Seq objects

    Returns:
        Float: Corrected distance

    Explanation:
        - Transitions (Ts): A↔G, C↔T (similar purine/pyrimidine)
        - Transversions (Tv): all other changes
        - Ts happen more often than Tv in real DNA
        - Kimura accounts for this difference
    """

    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be same length")

    transitions = 0
    transversions = 0

    for s1, s2 in zip(seq1, seq2):
        if s1 == s2 or s1 not in 'ACGT' or s2 not in 'ACGT':
            continue

        # Check if transition or transversion
        bases = {s1, s2}

        # Transitions: purine↔purine or pyrimidine↔pyrimidine
        if bases in [{'A', 'G'}, {'C', 'T'}]:
            transitions += 1
        else:
            transversions += 1

    total_changes = transitions + transversions
    if total_changes == 0:
        return 0.0

    n = len(seq1)
    p = transitions / n
    q = transversions / n

    # Kimura formula
    try:
        term1 = 1 - 2*p - q
        term2 = 1 - 2*q

        if term1 <= 0 or term2 <= 0:
            return float('inf')

        d = -0.5 * np.log(term1) - 0.25 * np.log(term2)
        return max(0, d)
    except:
        return float('inf')


###############################################################################
# FUNCTION: Calculate all pairwise distances from alignment
###############################################################################

def calculate_distance_matrix(alignment, method='p-distance'):
    """
    Calculate all pairwise distances from aligned sequences.

    Args:
        alignment: Bio.AlignIO alignment object
        method: 'p-distance', 'jukes-cantor', or 'kimura-2p'

    Returns:
        Dictionary: {(seq1_name, seq2_name): distance}
        List: All non-zero distances (for statistics)
    """

    # Choose distance function
    if method == 'p-distance':
        distance_func = calculate_p_distance
    elif method == 'jukes-cantor':
        distance_func = calculate_jukes_cantor
    elif method == 'kimura-2p':
        distance_func = calculate_kimura_2p
    else:
        raise ValueError(f"Unknown method: {method}")

    # Calculate all pairwise distances
    distances = {}
    distance_values = []

    sequences = list(alignment)
    for i, record1 in enumerate(sequences):
        for record2 in sequences[i+1:]:
            try:
                d = distance_func(record1.seq, record2.seq)
                distances[(record1.id, record2.id)] = d

                if d != float('inf') and not np.isnan(d):
                    distance_values.append(d)
            except Exception as e:
                print(f"Warning: Could not calculate distance between {record1.id} and {record2.id}: {e}",
                      file=sys.stderr)
                distances[(record1.id, record2.id)] = float('nan')

    return distances, distance_values


###############################################################################
# FUNCTION: Calculate distances from tree
###############################################################################

def calculate_tree_distances(treefile):
    """
    Calculate pairwise distances from tree branch lengths.

    Args:
        treefile: Path to Newick tree file

    Returns:
        Dictionary: {(seq1_name, seq2_name): distance}
        List: All distances (for statistics)
    """

    # Load tree
    tree = Phylo.read(treefile, 'newick')

    # Get all terminal nodes
    terminals = tree.get_terminals()

    # Calculate distances
    distances = {}
    distance_values = []

    for i, term1 in enumerate(terminals):
        for term2 in terminals[i+1:]:
            d = tree.distance(term1, term2)
            distances[(term1.name, term2.name)] = d
            distance_values.append(d)

    return distances, distance_values


###############################################################################
# FUNCTION: Calculate statistics
###############################################################################

def calculate_statistics(distances):
    """
    Calculate summary statistics for distance values.

    Args:
        distances: List of distance values

    Returns:
        Dictionary with statistics
    """

    distances = [d for d in distances if d != float('inf') and not np.isnan(d)]

    if not distances:
        return None

    stats = {
        'count': len(distances),
        'mean': np.mean(distances),
        'median': np.median(distances),
        'std': np.std(distances),
        'min': np.min(distances),
        'max': np.max(distances),
        'q25': np.percentile(distances, 25),
        'q75': np.percentile(distances, 75),
    }

    return stats


###############################################################################
# FUNCTION: Print distance matrix
###############################################################################

def print_distance_matrix(distances):
    """
    Print distance matrix in table format.

    Args:
        distances: Dictionary {(seq1, seq2): distance}
    """

    print("\n" + "="*70)
    print("PAIRWISE DISTANCE MATRIX")
    print("="*70)

    # Format: Seq1 vs Seq2 = Distance
    for (seq1, seq2), distance in sorted(distances.items()):
        if distance == float('inf'):
            print(f"{seq1:30s} vs {seq2:30s}: [undefined]")
        elif np.isnan(distance):
            print(f"{seq1:30s} vs {seq2:30s}: [NaN]")
        else:
            print(f"{seq1:30s} vs {seq2:30s}: {distance:.6f}")

    print("="*70)


###############################################################################
# FUNCTION: Save to CSV
###############################################################################

def save_distances_csv(distances, output_file):
    """Save distance matrix to CSV file."""
    try:
        with open(output_file, 'w') as f:
            f.write("Sequence_1,Sequence_2,Distance\n")
            for (seq1, seq2), distance in sorted(distances.items()):
                f.write(f"{seq1},{seq2},{distance}\n")
        print(f"✓ Distances saved to: {output_file}")
    except Exception as e:
        print(f"ERROR saving to CSV: {e}", file=sys.stderr)
        sys.exit(1)


###############################################################################
# FUNCTION: Threshold analysis
###############################################################################

def threshold_analysis(distances):
    """
    Analyze distance distribution for species threshold.

    DNA barcoding typically has:
    - Intraspecific: <0.003 (0.3%)
    - Interspecific: >0.03 (3%)
    - "Barcoding gap": space between them
    """

    print("\n" + "="*70)
    print("SPECIES THRESHOLD ANALYSIS")
    print("="*70)

    distances_list = [d for d in distances.values()
                      if d != float('inf') and not np.isnan(d)]

    if not distances_list:
        print("No valid distances to analyze")
        return

    stats = calculate_statistics(distances_list)

    print(f"\nDistance Statistics:")
    print(f"  Number of pairwise comparisons: {stats['count']}")
    print(f"  Mean distance:                  {stats['mean']:.6f} ({stats['mean']*100:.3f}%)")
    print(f"  Median distance:                {stats['median']:.6f} ({stats['median']*100:.3f}%)")
    print(f"  Std deviation:                  {stats['std']:.6f}")
    print(f"  Range:                          {stats['min']:.6f} to {stats['max']:.6f}")

    print(f"\nDNA Barcoding Interpretation:")
    print(f"  Species threshold (typical):    0.03 (3%)")
    print(f"  Intraspecific variation:        <0.003 (0.3%)")
    print(f"  Interspecific variation:        >0.03 (3%)")

    # Count distances in ranges
    threshold = 0.03
    intra_threshold = 0.003

    intra = sum(1 for d in distances_list if d < intra_threshold)
    inter = sum(1 for d in distances_list if d > threshold)
    gap = sum(1 for d in distances_list if intra_threshold <= d <= threshold)

    print(f"\nDistance Distribution:")
    print(f"  Distances < {intra_threshold} (intraspecific):   {intra} ({intra*100/len(distances_list):.1f}%)")
    print(f"  Distances {intra_threshold}-{threshold} (gap):       {gap} ({gap*100/len(distances_list):.1f}%)")
    print(f"  Distances > {threshold} (interspecific):   {inter} ({inter*100/len(distances_list):.1f}%)")

    if gap > 0:
        print(f"\n✓ Good barcoding gap detected!")
    else:
        print(f"\n⚠ Warning: May not be clear barcode gap (overlapping distances)")

    print("="*70)


###############################################################################
# FUNCTION: Command-line interface
###############################################################################

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '-a', '--alignment',
        metavar='FILE',
        help='Aligned FASTA file'
    )
    input_group.add_argument(
        '-t', '--tree',
        metavar='FILE',
        help='Phylogenetic tree file (Newick format)'
    )

    parser.add_argument(
        '-o', '--output',
        metavar='FILE',
        help='Save distances to CSV file'
    )

    parser.add_argument(
        '-m', '--method',
        choices=['p-distance', 'jukes-cantor', 'kimura-2p'],
        default='kimura-2p',
        help='Distance calculation method (default: kimura-2p)'
    )

    parser.add_argument(
        '--threshold-analysis',
        action='store_true',
        help='Analyze distance distribution for species threshold'
    )

    parser.add_argument(
        '--concepts',
        action='store_true',
        help='Print explanation of genetic distances'
    )

    return parser.parse_args()


###############################################################################
# MAIN
###############################################################################

def main():
    # Parse arguments
    args = parse_arguments()

    # Print concepts if requested
    if args.concepts:
        print(DISTANCE_CONCEPTS)
        return

    # Load data and calculate distances
    if args.alignment:
        validate_file(args.alignment)
        print(f"Loading alignment from: {args.alignment}")

        try:
            alignment = AlignIO.read(args.alignment, 'fasta')
            print(f"✓ Loaded {len(alignment)} sequences")
            print(f"  Sequence length: {alignment.get_alignment_length()}")

            print(f"\nCalculating {args.method} distances...")
            distances, distance_values = calculate_distance_matrix(alignment, args.method)
            print(f"✓ Calculated {len(distances)} pairwise distances")

        except Exception as e:
            print(f"ERROR loading alignment: {e}", file=sys.stderr)
            sys.exit(1)

    else:  # Tree input
        validate_file(args.tree)
        print(f"Loading tree from: {args.tree}")

        try:
            distances, distance_values = calculate_tree_distances(args.tree)
            print(f"✓ Calculated {len(distances)} pairwise distances from tree")
        except Exception as e:
            print(f"ERROR loading tree: {e}", file=sys.stderr)
            sys.exit(1)

    # Print distance matrix
    print_distance_matrix(distances)

    # Calculate and print statistics
    if distance_values:
        stats = calculate_statistics(distance_values)
        if stats:
            print(f"\nSummary Statistics:")
            print(f"  Mean:     {stats['mean']:.6f}")
            print(f"  Median:   {stats['median']:.6f}")
            print(f"  Std Dev:  {stats['std']:.6f}")
            print(f"  Range:    {stats['min']:.6f} to {stats['max']:.6f}")

    # Threshold analysis if requested
    if args.threshold_analysis:
        threshold_analysis(distances)

    # Save to file if requested
    if args.output:
        save_distances_csv(distances, args.output)

    # Print next steps
    print(f"\nFor species identification:")
    print(f"  Compare your unknown sequence distance to reference database")
    print(f"  If distance < 0.03 (3%) = likely same species")
    print(f"  If distance > 0.03 (3%) = likely different species")


if __name__ == '__main__':
    main()
