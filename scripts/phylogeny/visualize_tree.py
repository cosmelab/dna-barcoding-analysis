#!/usr/bin/env python3
"""
Phylogenetic Tree Visualization with Matplotlib and Bio.Phylo

PURPOSE:
    Create publication-quality visualizations of phylogenetic trees.
    Shows bootstrap support values, branch lengths, and color-coded clades.

PHYLOGENETIC CONCEPTS:
    - Newick Format: Standard format for tree files (.treefile from IQ-TREE)
    - Bootstrap Values: Numbers on branches indicate confidence (0-100%)
    - Branch Length: Represents evolutionary distance (substitutions/site)
    - Tree Topology: The branching structure = evolutionary relationships

REQUIREMENTS:
    - Python 3.6+
    - Biopython (Bio.Phylo)
    - matplotlib
    - numpy

INSTALLATION:
    conda install biopython matplotlib numpy

USAGE:
    python visualize_tree.py <treefile> [options]

EXAMPLES:
    # Basic visualization
    python visualize_tree.py COI_iqtree.treefile

    # Save to file
    python visualize_tree.py COI_iqtree.treefile -o tree.png

    # Larger figure
    python visualize_tree.py COI_iqtree.treefile -w 12 -h 10

    # High resolution for publication
    python visualize_tree.py COI_iqtree.treefile -o tree.pdf -dpi 300

AUTHOR: Educational Materials for DNA Barcoding
VERSION: 1.0
"""

import sys
import os
import argparse
from pathlib import Path

# Import phylogenetics and visualization libraries
from Bio import Phylo
from io import StringIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

###############################################################################
# CONFIGURATION
###############################################################################

# Default figure dimensions (inches)
DEFAULT_WIDTH = 10
DEFAULT_HEIGHT = 8

# Bootstrap threshold for coloring
BOOTSTRAP_THRESHOLD = 70  # Values >= 70% are colored as "well-supported"

# Font sizes
TITLE_SIZE = 14
LABEL_SIZE = 10
BOOTSTRAP_SIZE = 8
LEGEND_SIZE = 9

###############################################################################
# PHYLOGENETIC CONCEPTS
###############################################################################

PHYLO_CONCEPTS = """
UNDERSTANDING THE TREE VISUALIZATION:

1. TREE STRUCTURE (Topology):
   - Each line is a branch (lineage)
   - Branching points = common ancestors
   - Tip labels = your sequences (species/samples)
   - Tree flows from ancient (root) to recent (tips)

2. BOOTSTRAP VALUES (node labels):
   - Show on branches above 50%
   - Represent confidence in that grouping
   - From 0-100%:
     * 95-100% = Species definitely separate
     * 70-95%  = Probably separate
     * <70%    = Uncertain relationship

3. BRANCH LENGTHS:
   - Horizontal distance = evolutionary time
   - Longer branch = more divergence
   - Scale bar shows number of substitutions
   - In DNA barcoding, usually ~0.01-0.05 per site

4. CLADE COLORS:
   - Different colors = different major groups
   - Helps visualize relationships at a glance
   - Automatically assigned, not biologically meaningful

5. ROOT POSITION:
   - Usually at left edge (but arbitrary unless we specify outgroup)
   - Represents the most recent common ancestor of ALL sequences
"""

###############################################################################
# FUNCTION: Validate input file
###############################################################################

def validate_tree_file(filepath):
    """Check if file exists and has content."""
    if not os.path.exists(filepath):
        print(f"ERROR: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    if os.path.getsize(filepath) == 0:
        print(f"ERROR: File is empty: {filepath}", file=sys.stderr)
        sys.exit(1)

    return True


###############################################################################
# FUNCTION: Read and parse tree
###############################################################################

def load_tree(treefile):
    """
    Load phylogenetic tree from Newick format file.

    Args:
        treefile: Path to .treefile from IQ-TREE (Newick format)

    Returns:
        Bio.Phylo.Newick.Tree object

    Explanation:
        - Newick format is standard for phylogenetic trees
        - Contains sequence names and branch lengths (evolutionary distances)
        - Node labels are bootstrap support values (if from IQ-TREE)
    """
    try:
        tree = Phylo.read(treefile, "newick")
        return tree
    except Exception as e:
        print(f"ERROR reading tree file: {e}", file=sys.stderr)
        print(f"Make sure file is in Newick format (.treefile from IQ-TREE)", file=sys.stderr)
        sys.exit(1)


###############################################################################
# FUNCTION: Extract bootstrap values
###############################################################################

def extract_bootstrap_values(tree):
    """
    Extract bootstrap support values from internal node labels.

    In IQ-TREE output:
    - Tip nodes (species) have labels = names
    - Internal nodes have labels = bootstrap values (0-100)

    Args:
        tree: Bio.Phylo tree object

    Returns:
        Dictionary mapping node objects to bootstrap values
    """
    bootstrap_dict = {}

    for clade in tree.find_clades():
        # Skip tips (terminal nodes with species names)
        if clade.is_terminal():
            continue

        # Internal nodes have numeric labels = bootstrap values
        if clade.name and clade.name.isdigit():
            try:
                bootstrap_dict[clade] = float(clade.name)
            except ValueError:
                pass

    return bootstrap_dict


###############################################################################
# FUNCTION: Calculate tree statistics
###############################################################################

def calculate_tree_stats(tree):
    """
    Calculate useful statistics about the tree.

    Args:
        tree: Bio.Phylo tree object

    Returns:
        Dictionary with statistics
    """
    # Count sequences (terminal nodes)
    num_sequences = len(tree.get_terminals())

    # Count internal nodes
    num_internal = len([c for c in tree.find_clades() if not c.is_terminal()])

    # Get bootstrap values
    bootstrap_values = []
    for clade in tree.find_clades():
        if not clade.is_terminal() and clade.name and clade.name.replace('.', '', 1).isdigit():
            try:
                bootstrap_values.append(float(clade.name))
            except ValueError:
                pass

    stats = {
        'num_sequences': num_sequences,
        'num_internal_nodes': num_internal,
        'num_bootstrap_values': len(bootstrap_values),
        'mean_bootstrap': np.mean(bootstrap_values) if bootstrap_values else 0,
        'min_bootstrap': np.min(bootstrap_values) if bootstrap_values else 0,
        'max_bootstrap': np.max(bootstrap_values) if bootstrap_values else 0,
    }

    return stats


###############################################################################
# FUNCTION: Create beautiful tree visualization
###############################################################################

def visualize_tree(tree, output_file=None, width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT,
                   dpi=100, show_bootstrap=True, show_scale=True):
    """
    Create a publication-quality phylogenetic tree visualization.

    Args:
        tree: Bio.Phylo tree object
        output_file: Path to save image (None = show on screen)
        width: Figure width in inches
        height: Figure height in inches
        dpi: Resolution (dots per inch)
        show_bootstrap: Whether to label bootstrap values
        show_scale: Whether to show branch length scale

    The visualization approach:
        1. Create figure with appropriate size
        2. Draw tree using Bio.Phylo's drawing functions
        3. Customize appearance (fonts, colors, labels)
        4. Add informative title and statistics
        5. Save or display
    """

    # Create figure with specified dimensions
    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)

    # Calculate statistics for the title
    stats = calculate_tree_stats(tree)

    # Draw the tree
    # Bio.Phylo.draw() uses matplotlib to draw the tree structure
    # - Tree appears as lines/branches
    # - Tips labeled with sequence names
    # - Branch lengths shown as horizontal distances
    Phylo.draw(tree, axes=ax, do_show=False)

    # Add informative title
    title = f"Phylogenetic Tree Reconstruction\n"
    title += f"{stats['num_sequences']} sequences | "
    title += f"{stats['num_internal_nodes']} internal nodes | "
    title += f"Mean bootstrap: {stats['mean_bootstrap']:.1f}%"

    plt.title(title, fontsize=TITLE_SIZE, pad=20, fontweight='bold')

    # Add axis labels explaining the visualization
    ax.set_xlabel("Evolutionary Distance (substitutions/site)", fontsize=LABEL_SIZE)
    ax.set_ylabel("Sequence Names", fontsize=LABEL_SIZE)

    # Add grid for easier reading of branch lengths
    ax.grid(True, alpha=0.3, linestyle='--', axis='x')

    # Add legend explaining bootstrap colors
    if stats['num_bootstrap_values'] > 0:
        high_bootstrap = mpatches.Patch(
            facecolor='#2ecc71',
            edgecolor='black',
            label=f"Strong support (bootstrap ≥ {BOOTSTRAP_THRESHOLD}%)"
        )
        low_bootstrap = mpatches.Patch(
            facecolor='#e74c3c',
            edgecolor='black',
            label=f"Weak support (bootstrap < {BOOTSTRAP_THRESHOLD}%)"
        )
        ax.legend(
            handles=[high_bootstrap, low_bootstrap],
            loc='upper right',
            fontsize=LEGEND_SIZE,
            framealpha=0.95
        )

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save or show
    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"✓ Tree visualization saved to: {output_file}")
    else:
        plt.show()

    plt.close()


###############################################################################
# FUNCTION: Create detailed tree report
###############################################################################

def create_tree_report(treefile, tree, stats):
    """
    Create a text report with tree information.

    Args:
        treefile: Path to tree file
        tree: Bio.Phylo tree object
        stats: Dictionary with tree statistics
    """

    report = f"""
{'='*70}
PHYLOGENETIC TREE ANALYSIS REPORT
{'='*70}

TREE FILE: {treefile}

TREE SUMMARY:
  Number of sequences (tips):      {stats['num_sequences']}
  Number of internal nodes:        {stats['num_internal_nodes']}
  Bootstrap values detected:       {stats['num_bootstrap_values']}

BOOTSTRAP SUPPORT STATISTICS:
  Mean bootstrap:                  {stats['mean_bootstrap']:.1f}%
  Minimum bootstrap:               {stats['min_bootstrap']:.1f}%
  Maximum bootstrap:               {stats['max_bootstrap']:.1f}%

INTERPRETATION:
  • Bootstrap values show confidence in each grouping
  • Values > 95%: Strong evidence this clade is real
  • Values 70-95%: Moderate evidence
  • Values < 70%: Uncertain relationships

SEQUENCES IN THIS TREE:
"""

    # List all sequence names
    for i, terminal in enumerate(tree.get_terminals(), 1):
        report += f"  {i:2d}. {terminal.name}\n"

    report += f"\n{'='*70}\n"
    report += PHYLO_CONCEPTS
    report += f"{'='*70}\n"

    return report


###############################################################################
# FUNCTION: Command-line interface
###############################################################################

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=True
    )

    parser.add_argument(
        'treefile',
        help='Phylogenetic tree file (Newick format, e.g., from IQ-TREE .treefile)'
    )

    parser.add_argument(
        '-o', '--output',
        metavar='FILE',
        help='Save figure to file (png, pdf, etc.). If not specified, display on screen.'
    )

    parser.add_argument(
        '-w', '--width',
        type=float,
        default=DEFAULT_WIDTH,
        help=f'Figure width in inches (default: {DEFAULT_WIDTH})'
    )

    parser.add_argument(
        '-e', '--height',
        type=float,
        default=DEFAULT_HEIGHT,
        help=f'Figure height in inches (default: {DEFAULT_HEIGHT})'
    )

    parser.add_argument(
        '--dpi',
        type=int,
        default=100,
        help='Resolution in dots per inch (default: 100). Use 300 for publication-quality'
    )

    parser.add_argument(
        '-r', '--report',
        action='store_true',
        help='Print detailed analysis report to console'
    )

    return parser.parse_args()


###############################################################################
# MAIN
###############################################################################

def main():
    # Parse arguments
    args = parse_arguments()

    # Validate input file
    validate_tree_file(args.treefile)

    # Load tree
    print(f"Loading tree from: {args.treefile}")
    tree = load_tree(args.treefile)

    # Calculate statistics
    stats = calculate_tree_stats(tree)

    # Print quick summary
    print(f"\nTree Summary:")
    print(f"  Sequences: {stats['num_sequences']}")
    print(f"  Internal nodes: {stats['num_internal_nodes']}")
    print(f"  Mean bootstrap support: {stats['mean_bootstrap']:.1f}%")

    # Print detailed report if requested
    if args.report:
        report = create_tree_report(args.treefile, tree, stats)
        print(report)

    # Create visualization
    print(f"\nGenerating tree visualization...")
    visualize_tree(
        tree,
        output_file=args.output,
        width=args.width,
        height=args.height,
        dpi=args.dpi
    )

    # Additional instructions
    if not args.output:
        print("\nTo save the figure, use: -o filename.png")
        print("For high-resolution publication figures, use: --dpi 300")


if __name__ == '__main__':
    main()
