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
    """Generate tree visualization with Bio.Phylo with color-coded samples"""
    try:
        tree = Phylo.read(tree_file, "newick")

        # Identify which sequences are student samples vs references
        # References follow pattern: Genus_species_AccessionNumber
        # Anything else is a student sample
        def is_reference_sequence(name):
            """Check if this matches reference sequence pattern"""
            if not name:
                return False
            name = name.strip()

            # Reference pattern: Genus_species_AccessionNumber (e.g., Culex_pipiens_KP293422.1)
            # Must have at least 2 underscores and contain accession-like pattern
            parts = name.split('_')
            if len(parts) < 3:  # Need at least Genus_species_Accession
                return False

            # Check if last part contains accession number pattern (letters + numbers)
            last_part = parts[-1]
            # Accessions typically have format like: KP293422, KP293422.1, NC_036006, etc.
            has_letters = any(c.isalpha() for c in last_part)
            has_numbers = any(c.isdigit() for c in last_part)

            # Must have both letters and numbers to be an accession
            if has_letters and has_numbers:
                return True

            return False

        def is_student_sample(name):
            """Check if this is a student/tutorial sample (not a reference)"""
            if not name:
                return False
            # If it's not a reference, it's a student sample
            return not is_reference_sequence(name)

        # Color map for terminal nodes (leaves)
        def get_color_for_label(label):
            """Return color based on whether it's a sample or reference"""
            if is_student_sample(str(label)):
                return '#FF6B6B'  # Red/coral for student samples
            else:
                return '#4ECDC4'  # Teal for reference sequences

        # Create figure with more space
        fig, ax = plt.subplots(1, 1, figsize=(14, 10))

        # Draw tree first
        Phylo.draw(tree, axes=ax, do_show=False,
                   label_func=lambda x: x.name if x.name else '')

        # Now customize all text labels after drawing
        for text in ax.texts:
            label = text.get_text()
            if is_student_sample(label):
                # Student samples: Bold red labels, larger font
                text.set_fontweight('bold')
                text.set_fontsize(10)
                text.set_color('#FF6B6B')
            else:
                # Reference sequences: Normal teal labels, smaller font
                text.set_fontsize(8)
                text.set_color('#4ECDC4')

        # Add title
        ax.set_title("Phylogenetic Tree (Maximum Likelihood)\n" +
                     "Red = Your Samples  •  Teal = Reference Sequences",
                     fontsize=14, fontweight='bold', pad=20)

        # Add legend - position outside plot area to avoid blocking samples
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#FF6B6B', label='Your Samples (Student/Tutorial)'),
            Patch(facecolor='#4ECDC4', label='Reference Sequences (Known Species)')
        ]
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0, 1),
                  fontsize=10, framealpha=0.95, edgecolor='gray', fancybox=True)

        # Save figure
        plt.tight_layout()
        plt.savefig(output_image, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()

        print(f"✓ Tree visualization saved: {output_image}")
        return True
    except Exception as e:
        print(f"✗ Visualization failed: {e}")
        return False

def generate_html_report(tree_file, image_file, log_file, output_file):
    """Generate HTML report for phylogenetic tree using new design system"""
    from datetime import datetime
    from pathlib import Path

    # Try to extract model info from log
    model_info = "See .iqtree file for details"
    bootstrap_support = "1000 ultrafast"
    if log_file.exists():
        try:
            with open(log_file) as f:
                for line in f:
                    if "Best-fit model:" in line:
                        model_info = line.split(":")[-1].strip()
                        break
        except:
            pass

    # Read CSS files and embed them
    project_root = Path(__file__).parent.parent.parent
    base_css = (project_root / "tracking/styles/base.css").read_text()
    components_css = (project_root / "tracking/styles/components.css").read_text()
    reports_css = (project_root / "tracking/styles/reports.css").read_text()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Phylogenetic Tree Report - DNA Barcoding</title>

    <!-- Embedded CSS for reliable loading -->
    <style>
{base_css}

{components_css}

{reports_css}
    </style>
</head>
<body>
    <!-- Report Header -->
    <header class="report-header">
        <h1>🌳 Phylogenetic Tree Report</h1>
        <div class="progress-badge">Step 4 of 5</div>
        <div class="report-date">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</div>
    </header>

    <!-- Summary Dashboard -->
    <section class="summary-dashboard">
        <div class="metric-card metric-success">
            <div class="metric-value">IQ-TREE</div>
            <div class="metric-label">Analysis Method</div>
            <div class="metric-icon">⚙️</div>
        </div>

        <div class="metric-card metric-primary">
            <div class="metric-value">ML</div>
            <div class="metric-label">Maximum Likelihood</div>
            <div class="metric-icon">🧬</div>
        </div>

        <div class="metric-card metric-primary">
            <div class="metric-value">1000</div>
            <div class="metric-label">Bootstrap Replicates</div>
            <div class="metric-icon">🔁</div>
        </div>

        <div class="metric-card metric-success">
            <div class="metric-value">MFP</div>
            <div class="metric-label">Model Selection</div>
            <div class="metric-icon">📊</div>
        </div>
    </section>

    <!-- Main Content -->
    <main class="report-content">
        <div class="content-section">
            <h2>Phylogenetic Tree</h2>
            <p>Maximum likelihood tree showing evolutionary relationships</p>

            <div class="info-box info-tip">
                <strong>💡 Reading the Tree:</strong>
                <ul>
                    <li><strong>Branch lengths:</strong> Evolutionary distance (substitutions per site)</li>
                    <li><strong>Bootstrap values:</strong> Statistical support (≥70% considered reliable)</li>
                    <li><strong>Topology:</strong> Shows relationships, not time</li>
                    <li><strong>Your samples:</strong> Compare with reference species to identify unknown specimens</li>
                </ul>
            </div>

            <!-- Tree Visualization -->
            <div class="figure-container">
                <img src="{Path(image_file).name}" alt="Phylogenetic Tree" style="max-width: 100%; height: auto; border: 1px solid var(--border-color); border-radius: var(--radius-md); background: white;">
                <div class="figure-caption">
                    <strong>Figure 1:</strong> Maximum likelihood phylogenetic tree constructed using IQ-TREE with {bootstrap_support} bootstrap replicates.
                    Best-fit model: {model_info}
                </div>
            </div>

            <!-- Analysis Details -->
            <h3 style="margin-top: 2rem;">Analysis Details</h3>
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Parameter</th>
                        <th>Value</th>
                        <th>Description</th>
                    </tr>
                </thead>
                <tbody>
                    <tr class="row-pass">
                        <td><strong>Method</strong></td>
                        <td><span class="badge badge-primary">Maximum Likelihood</span></td>
                        <td>Statistical method for inferring phylogeny</td>
                    </tr>
                    <tr class="row-pass">
                        <td><strong>Software</strong></td>
                        <td><span class="badge badge-primary">IQ-TREE</span></td>
                        <td>State-of-the-art phylogenetic inference program</td>
                    </tr>
                    <tr class="row-pass">
                        <td><strong>Bootstrap</strong></td>
                        <td><span class="badge badge-success">1000 ultrafast replicates</span></td>
                        <td>Measures confidence in tree topology</td>
                    </tr>
                    <tr class="row-pass">
                        <td><strong>Model Selection</strong></td>
                        <td><span class="badge badge-primary">ModelFinder Plus (MFP)</span></td>
                        <td>Automatic selection of best-fit substitution model</td>
                    </tr>
                    <tr class="row-pass">
                        <td><strong>Best-fit Model</strong></td>
                        <td><code>{model_info}</code></td>
                        <td>DNA substitution model used for tree inference</td>
                    </tr>
                </tbody>
            </table>

            <!-- Output Files -->
            <h3 style="margin-top: 2rem;">Output Files</h3>
            <div class="info-box info-secondary">
                <strong>📁 Generated Files:</strong>
                <ul>
                    <li><code>{Path(tree_file).name}</code> - Newick tree format (text file, open in FigTree)</li>
                    <li><code>tree.png</code> - Tree visualization (this image)</li>
                    <li><code>tree.iqtree</code> - Detailed analysis report with statistics</li>
                    <li><code>tree.log</code> - IQ-TREE execution log</li>
                </ul>
            </div>

            <!-- Advanced Viewing -->
            <div class="info-box info-tip" style="margin-top: 1rem;">
                <strong>🔬 Advanced Viewing with FigTree:</strong>
                <p>For interactive tree exploration with bootstrap values and branch lengths:</p>
                <ol>
                    <li>Download <a href="http://tree.bio.ed.ac.uk/software/figtree/" target="_blank">FigTree</a> (free phylogenetic tree viewer)</li>
                    <li>Open <code>{Path(tree_file).name}</code> in FigTree</li>
                    <li>Enable "Node Labels" to display bootstrap support values</li>
                    <li>Adjust "Branch Labels" to show evolutionary distances</li>
                    <li>Export high-resolution images for publications</li>
                </ol>
            </div>
        </div>
    </main>

    <!-- Footer with Help -->
    <footer class="report-footer">
        <div class="help-section">
            <h3>Understanding Phylogenetic Trees</h3>

            <h4>What is a phylogenetic tree?</h4>
            <p>A phylogenetic tree is a diagram showing the evolutionary relationships between organisms based on DNA sequence similarity. Closely related species appear near each other on the tree.</p>

            <h4>How to interpret the tree</h4>
            <ul>
                <li><strong>Tips (leaves):</strong> Current species or samples being analyzed</li>
                <li><strong>Branches:</strong> Evolutionary lineages; longer branches = more changes</li>
                <li><strong>Nodes:</strong> Common ancestors where lineages split</li>
                <li><strong>Root:</strong> Most recent common ancestor of all sequences</li>
            </ul>

            <h4>What are bootstrap values?</h4>
            <p>Bootstrap values (0-100%) measure confidence in tree branching patterns:</p>
            <ul>
                <li><strong>≥95%:</strong> Very strong support</li>
                <li><strong>70-95%:</strong> Good support</li>
                <li><strong>&lt;70%:</strong> Weak support (uncertain relationship)</li>
            </ul>

            <h4>What is Maximum Likelihood?</h4>
            <p>Maximum Likelihood (ML) is a statistical method that finds the tree topology and branch lengths most likely to have produced the observed sequences, given a model of DNA evolution.</p>

            <h4>Why use IQ-TREE?</h4>
            <p>IQ-TREE is the fastest and most accurate ML phylogenetic software available. Key features:</p>
            <ul>
                <li>Automatic model selection (ModelFinder)</li>
                <li>Ultrafast bootstrap (1000 replicates in seconds)</li>
                <li>Handles large datasets efficiently</li>
                <li>Widely used in published research</li>
            </ul>

            <h4>How to use this tree</h4>
            <ol>
                <li>Find your samples on the tree (student sequences)</li>
                <li>Look at which reference species cluster with your samples</li>
                <li>Check bootstrap values for branch reliability</li>
                <li>Compare results with BLAST identification (next step)</li>
            </ol>
        </div>
    </footer>
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
