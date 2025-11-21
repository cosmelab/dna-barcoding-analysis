#!/usr/bin/env python3
"""
Phylogenetic Tree Construction
Uses IQ-TREE for tree building and Bio.Phylo for visualization
"""

import os
import sys
import subprocess
from pathlib import Path
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

def run_iqtree(alignment_file, output_prefix):
    """Run IQ-TREE for phylogenetic inference"""
    print(f"Running IQ-TREE (this may take a few minutes)...")

    try:
        result = subprocess.run(
            [
                'iqtree',
                '-s', str(alignment_file),
                '-pre', str(output_prefix),
                '-m', 'MFP',  # Model Finder Plus - auto-select best model
                '-bb', '1000',  # 1000 ultrafast bootstrap replicates
                '-nt', 'AUTO',  # Auto-detect number of threads
                '-redo'  # Redo analysis (overwrite existing files)
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        print(f"✓ Tree construction complete")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ IQ-TREE failed: {e.stderr}")
        return False
    except FileNotFoundError:
        print("✗ IQ-TREE not found. Please install IQ-TREE.")
        return False

def visualize_tree(tree_file, output_image):
    """Generate tree visualization with Bio.Phylo"""
    try:
        tree = Phylo.read(tree_file, "newick")

        # Create figure
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))

        # Draw tree
        Phylo.draw(tree, axes=ax, do_show=False)
        ax.set_title("Phylogenetic Tree (Maximum Likelihood)", fontsize=14, fontweight='bold')

        # Save figure
        plt.tight_layout()
        plt.savefig(output_image, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"✓ Tree visualization saved: {output_image}")
        return True
    except Exception as e:
        print(f"✗ Visualization failed: {e}")
        return False

def generate_html_report(tree_file, image_file, log_file, output_file):
    """Generate simple HTML report"""
    from datetime import datetime

    # Try to extract model info from log
    model_info = "See .iqtree file for details"
    if log_file.exists():
        try:
            with open(log_file) as f:
                for line in f:
                    if "Best-fit model:" in line:
                        model_info = line.split(":")[-1].strip()
                        break
        except:
            pass

    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Phylogenetic Tree Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #333; }}
        .summary {{ background-color: #f0f0f0; padding: 15px; margin: 20px 0; border-radius: 5px; }}
        .tree-img {{ max-width: 100%; border: 1px solid #ddd; margin: 20px 0; }}
        .info {{ background-color: #e3f2fd; padding: 10px; margin: 10px 0; border-left: 4px solid #2196F3; }}
    </style>
</head>
<body>
    <h1>Phylogenetic Tree Report</h1>
    <p>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>

    <div class="summary">
        <h2>Analysis Summary</h2>
        <p><strong>Method:</strong> Maximum Likelihood (IQ-TREE)</p>
        <p><strong>Bootstrap:</strong> 1000 ultrafast bootstrap replicates</p>
        <p><strong>Model:</strong> {model_info}</p>
        <p><strong>Tree file:</strong> {Path(tree_file).name}</p>
    </div>

    <h2>Phylogenetic Tree</h2>
    <img src="{Path(image_file).name}" alt="Phylogenetic Tree" class="tree-img">

    <div class="info">
        <h3>Opening in FigTree (Optional)</h3>
        <p>For a more detailed view with bootstrap values:</p>
        <ol>
            <li>Download and install <a href="http://tree.bio.ed.ac.uk/software/figtree/">FigTree</a></li>
            <li>Open the tree file: <code>{Path(tree_file).name}</code></li>
            <li>Enable "Node Labels" to see bootstrap support values</li>
        </ol>
    </div>

    <h2>Output Files</h2>
    <ul>
        <li><strong>{Path(tree_file).name}</strong> - Newick tree format (open in FigTree)</li>
        <li><strong>tree.png</strong> - Tree visualization</li>
        <li><strong>tree.iqtree</strong> - Detailed analysis report</li>
        <li><strong>tree.log</strong> - IQ-TREE run log</li>
    </ul>
</body>
</html>"""

    with open(output_file, 'w') as f:
        f.write(html)

def main():
    """Main phylogeny function"""
    if len(sys.argv) < 2:
        print("Usage: python build_tree.py <aligned_fasta> [output_directory]")
        print("Example: python build_tree.py results/aligned_sequences.fasta results/")
        sys.exit(1)

    input_fasta = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results")
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_fasta.exists():
        print(f"Input file not found: {input_fasta}")
        sys.exit(1)

    print(f"Input alignment: {input_fasta}")

    # Run IQ-TREE
    output_prefix = output_dir / "tree"
    if not run_iqtree(input_fasta, output_prefix):
        sys.exit(1)

    # IQ-TREE creates these files
    tree_file = output_prefix.with_suffix('.treefile')
    log_file = output_prefix.with_suffix('.log')
    iqtree_file = output_prefix.with_suffix('.iqtree')

    # Visualize tree
    print("\nGenerating tree visualization...")
    image_file = output_dir / "tree.png"
    visualize_tree(tree_file, image_file)

    # Generate HTML report
    print("\nGenerating HTML report...")
    html_file = output_dir / "phylogeny_report.html"
    generate_html_report(tree_file, image_file, log_file, html_file)
    print(f"HTML report: {html_file}")

    print(f"\nPhylogeny Analysis Complete!")
    print(f"\nOutput files in {output_dir}:")
    print(f"  - tree.treefile - Newick format (open in FigTree)")
    print(f"  - tree.png - Simple visualization")
    print(f"  - phylogeny_report.html - HTML report")
    print(f"  - tree.iqtree - Detailed analysis log")

if __name__ == "__main__":
    main()
