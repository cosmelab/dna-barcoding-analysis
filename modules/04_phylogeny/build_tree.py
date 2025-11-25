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

            # Accessions can have underscores (NC_054318) or dots (KP293422.1)
            # Check if the last 1-2 parts look like an accession:
            # 1. Last part has letters+numbers (KP293422, KP293422.1)
            # 2. Last part is only numbers AND second-to-last has letters (NC_054318 → NC + 054318)
            # 3. Last part is only numbers with dot (036006.1) AND prior part has letters

            last_part = parts[-1]
            has_letters_last = any(c.isalpha() for c in last_part)
            has_numbers_last = any(c.isdigit() for c in last_part)

            # Case 1: Last part has both letters and numbers (e.g., KP293422.1, MH376751.1)
            if has_letters_last and has_numbers_last:
                return True

            # Case 2: Last part is numeric-only, check if second-to-last has letters
            # This catches: NC_054318, NC_036006.1, etc.
            if has_numbers_last and len(parts) >= 3:
                second_to_last = parts[-2]
                if any(c.isalpha() for c in second_to_last):
                    return True

            return False

        def is_student_sample(name):
            """Check if this is a student/tutorial sample (not a reference)"""
            if not name:
                return False
            # If it's not a reference, it's a student sample
            return not is_reference_sequence(name)

        # Extract all unique genera from reference sequences
        genera = set()
        for terminal in tree.get_terminals():
            if terminal.name and is_reference_sequence(terminal.name):
                parts = terminal.name.strip().split('_')
                if len(parts) >= 1:
                    genus = parts[0]
                    genera.add(genus)

        # Create color palette for genera (avoid red tones - reserved for student samples)
        genus_colors = {
            'Aedes': '#AA96DA',      # Purple
            'Anopheles': '#4ECDC4',  # Teal
            'Culex': '#95E1D3',      # Mint
            'Deinocerites': '#A8E6CF',   # Seafoam
            'Psorophora': '#8B7FD8',  # Lavender
            'Uranotaenia': '#9BDEAC', # Light green
            'Wyeomyia': '#5DADE2',   # Sky blue
            'Ochlerotatus': '#85C1E2', # Light blue
            'Toxorhynchites': '#AED6F1', # Pale blue
            'Mansonia': '#B4E7CE',   # Pale green
            'Culiseta': '#6DC5D1',   # Turquoise
        }

        # Color map for terminal nodes (leaves)
        def get_color_for_label(label):
            """Return color based on genus or if it's a student sample"""
            if is_student_sample(str(label)):
                return '#FF6B6B'  # Bright red for student samples
            else:
                # Extract genus from reference sequence
                parts = str(label).strip().split('_')
                if len(parts) >= 1:
                    genus = parts[0]
                    return genus_colors.get(genus, '#808080')  # Gray as fallback
                return '#808080'  # Gray fallback

        # Create figure with more space
        fig, ax = plt.subplots(1, 1, figsize=(14, 10))

        # Draw tree first
        Phylo.draw(tree, axes=ax, do_show=False,
                   label_func=lambda x: x.name if x.name else '')

        # Get all terminal node names for filtering
        terminal_names = {str(term.name).strip() for term in tree.get_terminals() if term.name}

        # Now customize ONLY terminal labels (not bootstrap values or other text)
        for text in ax.texts:
            label = text.get_text().strip()

            # Skip if this text is not a terminal node label
            if label not in terminal_names:
                continue

            if is_student_sample(label):
                # Student samples: Bold red labels, larger font
                text.set_fontweight('bold')
                text.set_fontsize(10)
                text.set_color('#FF6B6B')
            else:
                # Reference sequences: Color by genus
                parts = label.split('_')
                if len(parts) >= 1:
                    genus = parts[0]
                    color = genus_colors.get(genus, '#808080')
                    text.set_fontsize(8)
                    text.set_color(color)
                else:
                    text.set_fontsize(8)
                    text.set_color('#808080')

        # Add title
        ax.set_title("Phylogenetic Tree (Maximum Likelihood)\n" +
                     "Red = Your Samples  •  Colors = Reference Genera",
                     fontsize=14, fontweight='bold', pad=20)

        # Build legend dynamically based on genera present in tree
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#FF6B6B', label='Your Samples (Student/Tutorial)', edgecolor='darkred', linewidth=1.5)
        ]

        # Add genera in alphabetical order
        for genus in sorted(genera):
            color = genus_colors.get(genus, '#808080')
            legend_elements.append(
                Patch(facecolor=color, label=f'{genus} (genus)', edgecolor='gray', linewidth=0.5)
            )

        # Position legend outside plot area to avoid blocking samples
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1),
                  fontsize=9, framealpha=0.95, edgecolor='gray', fancybox=True,
                  title='Species Groups', title_fontsize=10)

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
                    <li><strong>Your samples (red, bold):</strong> Compare with reference species to identify unknown specimens</li>
                    <li><strong>Reference sequences (colored by genus):</strong> Each mosquito genus has its own color for easy visual grouping</li>
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

            <!-- Optional: Tree Abbreviation -->
            <div class="info-box info-secondary" style="margin-top: 1rem;">
                <strong>📝 Optional: Shorten Species Names</strong>
                <p>If your tree labels are crowded or you want a cleaner appearance for publications:</p>
                <pre style="background: var(--bg-code); padding: 0.5rem; border-radius: var(--radius-sm); margin: 0.5rem 0;">python scripts/abbreviate_tree_names.py {Path(tree_file)}</pre>
                <p style="margin-bottom: 0; font-size: 0.9rem;">This creates an abbreviated version with shorter names (e.g., <code>Culex_pipiens_KP293422.1</code> → <code>cu_pip_KP29</code>) plus a reference table. Your sample names stay unchanged. See <code>docs/optional_tree_abbreviation.md</code> for details.</p>
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
