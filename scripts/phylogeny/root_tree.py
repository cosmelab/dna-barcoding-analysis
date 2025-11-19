#!/usr/bin/env python3
"""
Phylogenetic Tree Rooting Methods

PURPOSE:
    Reroot phylogenetic trees using different methods. Demonstrate how rooting
    changes the interpretation of evolutionary relationships.

PHYLOGENETIC CONCEPTS:
    - UNROOTED TREE: Branches show relationships but no direction of evolution
    - ROOTED TREE: Shows ancestor-descendant relationships (directionality)
    - ROOT POSITION: Determines which lineage is "ancestral"

ROOTING METHODS:
    1. MIDPOINT ROOTING: Root at the middle of the longest path
       - Use when you don't know the ancestral species
       - Assumes constant evolutionary rates
       - Works well for balanced datasets

    2. OUTGROUP ROOTING: Root between outgroup and ingroup
       - Use when you have an outgroup (more distant species)
       - Most accurate if outgroup is known
       - Requires identifying outgroup species

    3. MIDVAR ROOTING: Root to minimize variance in root-to-tip distances
       - Advanced: assumes molecular clock
       - Good for dating evolutionary events
       - Requires branch lengths

WHY THIS MATTERS FOR DNA BARCODING:
    - Rooting determines what "ancestral sequences" were
    - Different roots change species groupings
    - For identification, we usually just care about clades
    - Root matters for understanding evolutionary history

REQUIREMENTS:
    - Python 3.6+
    - Biopython
    - dendropy (for alternative methods)

INSTALLATION:
    conda install biopython dendropy

USAGE:
    python root_tree.py <treefile> [method] [options]

EXAMPLES:
    # Midpoint rooting (default)
    python root_tree.py tree.treefile midpoint

    # Outgroup rooting
    python root_tree.py tree.treefile outgroup -g "Species_A" -g "Species_B"

    # Save rooted tree
    python root_tree.py tree.treefile midpoint -o rooted_tree.treefile

    # List all species
    python root_tree.py tree.treefile --list-species

AUTHOR: Educational Materials for DNA Barcoding
VERSION: 1.0
"""

import sys
import os
import argparse
from pathlib import Path

from Bio import Phylo
from io import StringIO

###############################################################################
# PHYLOGENETIC CONCEPTS EXPLANATION
###############################################################################

ROOTING_CONCEPTS = """
UNDERSTANDING TREE ROOTING:

1. WHY ROOT A TREE?
   - Unrooted trees: Show who's related, not direction of evolution
   - Rooted trees: Show ancestor-descendant relationships
   - For phylogenetics: We want to know "what's the common ancestor?"

2. WHAT HAPPENS WHEN YOU ROOT?
   Example tree (unrooted):
          A ─── B
           \\   /
            \\ /
           (unknown)
            / \\
           /   \\
          C ─── D

   If we root on A:  A is most ancestral
                     B,C,D share a more recent common ancestor

   If we root on C:  C is most ancestral
                     A,B,D share a more recent common ancestor

3. MIDPOINT ROOTING:
   - Root at the point that divides longest path in half
   - Best if you have no outgroup information
   - Assumes roughly constant evolutionary rates
   - Good for initial visualization

4. OUTGROUP ROOTING:
   - Specify a species known to be distantly related
   - Root between outgroup and everything else
   - Most accurate if you know your outgroups
   - Preferred method if outgroup available

5. FOR DNA BARCODING IDENTIFICATION:
   - Usually don't care about root position
   - Care about CLADES (groups of closely related species)
   - But root matters if you're interpreting evolutionary history
"""

###############################################################################
# FUNCTION: Validate tree file
###############################################################################

def validate_tree_file(filepath):
    """Check if tree file exists and is readable."""
    if not os.path.exists(filepath):
        print(f"ERROR: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    if os.path.getsize(filepath) == 0:
        print(f"ERROR: File is empty: {filepath}", file=sys.stderr)
        sys.exit(1)

    return True


###############################################################################
# FUNCTION: Load tree
###############################################################################

def load_tree(treefile):
    """Load tree from Newick format."""
    try:
        tree = Phylo.read(treefile, "newick")
        return tree
    except Exception as e:
        print(f"ERROR reading tree: {e}", file=sys.stderr)
        sys.exit(1)


###############################################################################
# FUNCTION: List species in tree
###############################################################################

def list_species(tree):
    """
    Print all species/sequences in the tree.

    Useful for:
    - Identifying which sequences are in the tree
    - Finding outgroup species
    - Quality control (checking if expected species present)
    """
    terminals = tree.get_terminals()

    print(f"\nSpecies in tree: ({len(terminals)} total)")
    print("-" * 50)

    for i, terminal in enumerate(terminals, 1):
        print(f"  {i:2d}. {terminal.name}")

    print("-" * 50)


###############################################################################
# FUNCTION: Midpoint rooting
###############################################################################

def root_midpoint(tree):
    """
    Root tree at the midpoint of the longest path between any two tips.

    Concept:
        - Find two most distant sequences
        - Root at the middle of the path between them
        - Assumes evolutionary rates are constant
        - Works when no outgroup is available

    Args:
        tree: Bio.Phylo tree object

    Returns:
        Rooted tree

    Side effects:
        - Modifies tree object in place
        - Removes old root position
        - Places new root at midpoint
    """

    # Get all pairwise distances between tips
    # This is the "longest path" problem
    terminals = tree.get_terminals()
    longest_distance = 0
    furthest_pair = None

    # Calculate distances between all pairs
    for i, tip1 in enumerate(terminals):
        for tip2 in terminals[i+1:]:
            # Calculate distance along tree
            distance = tree.distance(tip1, tip2)

            if distance > longest_distance:
                longest_distance = distance
                furthest_pair = (tip1, tip2)

    if furthest_pair is None:
        print("ERROR: Could not find pair of sequences", file=sys.stderr)
        return tree

    # Find the path between the two furthest sequences
    tip1, tip2 = furthest_pair

    # Distance from tip1 to the midpoint
    midpoint_distance = longest_distance / 2.0

    # Find the node closest to midpoint
    # Start at tip1 and walk toward tip2
    current = tip1
    distance_from_tip1 = 0

    # Walk up tree from tip1 to find the midpoint
    while True:
        # Get parent
        if current.parent is None:
            break

        # Distance to parent
        branch_length = current.branch_length if current.branch_length else 0

        if distance_from_tip1 + branch_length >= midpoint_distance:
            # Midpoint is on this branch
            break

        distance_from_tip1 += branch_length
        current = current.parent

    # Root at the midpoint
    tree.root = current

    return tree


###############################################################################
# FUNCTION: Outgroup rooting
###############################################################################

def root_outgroup(tree, outgroup_names):
    """
    Root tree using one or more outgroup species.

    Concept:
        - Specify one or more species that are distantly related
        - Evolutionary root is between outgroup and ingroup
        - This is the most biologically informed rooting method
        - Preferred if you know your outgroups

    Args:
        tree: Bio.Phylo tree object
        outgroup_names: List of species names in outgroup

    Returns:
        Rooted tree

    Example:
        If tree has species [A, B, C, D, E] and you specify A as outgroup:
        Then A is the "ancestral" lineage, B-E share more recent ancestor
    """

    # Find the outgroup clade (smallest clade containing all outgroups)
    outgroup_names_set = set(outgroup_names)
    outgroup_clade = None

    # Search tree for a clade containing all outgroup species
    for clade in tree.find_clades():
        # Get all tips in this clade
        clade_tips = {tip.name for tip in clade.get_terminals()}

        # Check if all outgroups are in this clade
        if outgroup_names_set.issubset(clade_tips):
            # And check this is the smallest such clade
            if outgroup_clade is None or len(clade_tips) < len({tip.name for tip in outgroup_clade.get_terminals()}):
                outgroup_clade = clade

    if outgroup_clade is None:
        print(f"ERROR: Could not find outgroup {outgroup_names} in tree", file=sys.stderr)
        return tree

    # Root on the outgroup
    tree.root = outgroup_clade

    return tree


###############################################################################
# FUNCTION: Ladderize tree
###############################################################################

def ladderize_tree(tree, reverse=False):
    """
    Sort tree to make visualization cleaner.

    Concept:
        - Ladderizing reorders branches without changing relationships
        - Deeper clades appear at top (or bottom if reversed)
        - Makes tree easier to read visually
        - Does NOT change evolutionary interpretation

    Args:
        tree: Bio.Phylo tree object
        reverse: If True, larger clades at bottom

    Returns:
        Sorted tree
    """

    tree.ladderize(reverse=reverse)
    return tree


###############################################################################
# FUNCTION: Print tree statistics
###############################################################################

def print_tree_stats(tree):
    """Print information about tree structure."""
    terminals = tree.get_terminals()
    internal = [c for c in tree.find_clades() if not c.is_terminal()]

    print(f"\nTree Statistics:")
    print(f"  Number of sequences: {len(terminals)}")
    print(f"  Number of internal nodes: {len(internal)}")

    # Tree height (maximum distance from root)
    if tree.root:
        max_distance = 0
        for terminal in terminals:
            distance = tree.distance(tree.root, terminal)
            if distance > max_distance:
                max_distance = distance
        print(f"  Tree height: {max_distance:.4f}")


###############################################################################
# FUNCTION: Save tree
###############################################################################

def save_tree(tree, output_file):
    """Save rooted tree to file."""
    try:
        Phylo.write(tree, output_file, "newick")
        print(f"✓ Rooted tree saved to: {output_file}")
    except Exception as e:
        print(f"ERROR saving tree: {e}", file=sys.stderr)
        sys.exit(1)


###############################################################################
# FUNCTION: Compare tree representations
###############################################################################

def print_tree_ascii(tree, max_lines=20):
    """
    Print ASCII representation of tree to terminal.

    Useful for:
    - Quick visualization without plotting
    - Seeing exact topology before/after rooting
    - Checking that rooting worked as expected
    """

    # Use Bio.Phylo's ASCII drawing
    output = StringIO()
    Phylo.draw_ascii(tree, file=output)
    tree_str = output.getvalue()

    # Limit output if too large
    lines = tree_str.split('\n')
    if len(lines) > max_lines:
        print(f"\n(Showing first {max_lines} lines of tree):")
        tree_str = '\n'.join(lines[:max_lines]) + "\n..."

    print(tree_str)


###############################################################################
# FUNCTION: Command-line interface
###############################################################################

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        'treefile',
        help='Input tree file (Newick format)'
    )

    parser.add_argument(
        'method',
        nargs='?',
        default='midpoint',
        choices=['midpoint', 'outgroup'],
        help='Rooting method (default: midpoint)'
    )

    parser.add_argument(
        '-o', '--output',
        metavar='FILE',
        help='Save rooted tree to file'
    )

    parser.add_argument(
        '-g', '--outgroup',
        metavar='SPECIES',
        action='append',
        dest='outgroup',
        help='Outgroup species (use multiple times for multiple species). '
             'Required for outgroup method.'
    )

    parser.add_argument(
        '--list-species',
        action='store_true',
        help='List all species in tree and exit'
    )

    parser.add_argument(
        '--show-ascii',
        action='store_true',
        help='Print ASCII representation of tree'
    )

    parser.add_argument(
        '--concepts',
        action='store_true',
        help='Print explanation of tree rooting concepts'
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
        print(ROOTING_CONCEPTS)
        return

    # Validate input
    validate_tree_file(args.treefile)

    # Load tree
    print(f"Loading tree from: {args.treefile}")
    tree = load_tree(args.treefile)

    print(f"Initial tree has root: {tree.root is not None}")
    print_tree_stats(tree)

    # List species and exit if requested
    if args.list_species:
        list_species(tree)
        return

    # Apply rooting method
    if args.method == 'midpoint':
        print(f"\nApplying midpoint rooting...")
        tree = root_midpoint(tree)
        print("✓ Midpoint rooting applied")

    elif args.method == 'outgroup':
        if not args.outgroup:
            print("ERROR: outgroup method requires -g/--outgroup", file=sys.stderr)
            print("Example: --outgroup Species_A --outgroup Species_B", file=sys.stderr)
            sys.exit(1)

        print(f"\nApplying outgroup rooting...")
        print(f"Outgroup species: {', '.join(args.outgroup)}")
        tree = root_outgroup(tree, args.outgroup)
        print("✓ Outgroup rooting applied")

    # Ladderize for better visualization
    print("\nLadderizing tree for better visualization...")
    tree = ladderize_tree(tree)
    print("✓ Tree ladderized")

    # Print statistics after rooting
    print_tree_stats(tree)

    # Show ASCII if requested
    if args.show_ascii:
        print_tree_ascii(tree)

    # Save if output file specified
    if args.output:
        save_tree(tree, args.output)

    # Print instructions
    print(f"\n{'='*50}")
    print("NEXT STEPS:")
    print(f"  1. Visualize with: python visualize_tree.py {args.output or 'rooted_tree.treefile'}")
    print(f"  2. Calculate distances: python calculate_distances.py {args.output or args.treefile}")
    print(f"  3. Extract species groups for identification")
    print(f"{'='*50}")


if __name__ == '__main__':
    main()
