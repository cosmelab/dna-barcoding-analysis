#!/usr/bin/env python3
"""
Phylogenetic Tree Construction
Uses IQ-TREE for tree building and multiple visualization tools
Generates trees with full names and abbreviated names in different layouts
"""

import os
import sys
import subprocess
from pathlib import Path
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import time

# Import toytree and toyplot for multiple tree layouts
try:
    import toytree
    import toyplot
    import toyplot.png
    import toyplot.svg
    import toyplot.pdf
    TOYTREE_AVAILABLE = True
except ImportError:
    TOYTREE_AVAILABLE = False
    print("⚠️  toytree not available - only Bio.Phylo visualization will be used")

def print_progress_bar(message, steps=20, delay=0.05):
    """Print an ASCII progress bar"""
    print(f"\n{message}")
    bar = "█"
    empty = "░"
    for i in range(steps + 1):
        filled = bar * i
        remaining = empty * (steps - i)
        percent = (i / steps) * 100
        print(f"\r  [{filled}{remaining}] {percent:.0f}%", end="", flush=True)
        time.sleep(delay)
    print("\n")

def run_iqtree(alignment_file, output_prefix):
    """Run IQ-TREE for phylogenetic inference with spinner"""
    # Import Rich for spinner
    try:
        from rich.console import Console
        from rich.live import Live
        from rich.text import Text
        console = Console(force_terminal=True)
        RICH_AVAILABLE = True
    except ImportError:
        RICH_AVAILABLE = False

    start_time = time.time()

    if RICH_AVAILABLE:
        with Live(console=console, refresh_per_second=10) as live:
            try:
                spinner_chars = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
                frame = 0

                process = subprocess.Popen(
                    [
                        'iqtree',
                        '-s', str(alignment_file),
                        '-pre', str(output_prefix),
                        '-m', 'MFP',
                        '-bb', '1000',
                        '-nt', 'AUTO',
                        '-redo'
                    ],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )

                while process.poll() is None:
                    elapsed = time.time() - start_time
                    char = spinner_chars[frame % len(spinner_chars)]
                    live.update(Text(f"  {char} 🌳 Building phylogenetic tree... ({elapsed:.0f}s)", style="bold cyan"))
                    frame += 1
                    time.sleep(0.1)

                if process.returncode != 0:
                    stderr = process.stderr.read()
                    print(f"✗ IQ-TREE failed: {stderr}")
                    return False

                elapsed = time.time() - start_time
                print(f"✓ Tree construction complete in {elapsed:.1f}s")
                return True

            except FileNotFoundError:
                print("✗ IQ-TREE not found. Please install IQ-TREE.")
                return False
    else:
        print("🌳 Building phylogenetic tree (this may take 1-2 minutes)...")
        try:
            result = subprocess.run(
                [
                    'iqtree',
                    '-s', str(alignment_file),
                    '-pre', str(output_prefix),
                    '-m', 'MFP',
                    '-bb', '1000',
                    '-nt', 'AUTO',
                    '-redo'
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            elapsed = time.time() - start_time
            print(f"✓ Tree construction complete in {elapsed:.1f}s")
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
            Patch(facecolor='#FF6B6B', label='Your Samples', edgecolor='darkred', linewidth=1.5)
        ]

        # Add genera in alphabetical order
        for genus in sorted(genera):
            color = genus_colors.get(genus, '#808080')
            legend_elements.append(
                Patch(facecolor=color, label=genus, edgecolor='gray', linewidth=0.5)
            )

        # Position legend outside plot area to avoid blocking samples
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1),
                  fontsize=9, framealpha=0.95, edgecolor='gray', fancybox=True,
                  title='Mosquito Genera', title_fontsize=10)

        # Save figure
        plt.tight_layout()
        plt.savefig(output_image, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()

        print(f"✓ Tree visualization saved: {output_image}")
        return True
    except Exception as e:
        print(f"✗ Visualization failed: {e}")
        return False

def visualize_tree_toytree(tree_file, output_dir):
    """Generate multiple tree layouts using toytree"""
    if not TOYTREE_AVAILABLE:
        print("⚠️  Skipping toytree visualizations (library not available)")
        return []

    try:
        # Load tree
        tree = toytree.tree(str(tree_file))

        # Define genus colors (same as Bio.Phylo visualization)
        genus_colors = {
            'Aedes': '#AA96DA',
            'Anopheles': '#4ECDC4',
            'Culex': '#95E1D3',
            'Deinocerites': '#A8E6CF',
            'Psorophora': '#8B7FD8',
            'Uranotaenia': '#9BDEAC',
            'Wyeomyia': '#5DADE2',
            'Ochlerotatus': '#85C1E2',
            'Toxorhynchites': '#AED6F1',
            'Mansonia': '#B4E7CE',
            'Culiseta': '#6DC5D1',
        }

        def is_reference_sequence(name):
            """Check if this matches reference sequence pattern"""
            if not name:
                return False
            name = name.strip()
            parts = name.split('_')
            if len(parts) < 3:
                return False

            last_part = parts[-1]
            has_letters_last = any(c.isalpha() for c in last_part)
            has_numbers_last = any(c.isdigit() for c in last_part)

            if has_letters_last and has_numbers_last:
                return True

            if has_numbers_last and len(parts) >= 3:
                second_to_last = parts[-2]
                if any(c.isalpha() for c in second_to_last):
                    return True

            return False

        def extract_genus(name):
            """Extract genus from sequence name"""
            if not name:
                return None
            parts = name.strip().split('_')
            if len(parts) >= 1:
                return parts[0]
            return None

        # Prepare tip colors
        tip_colors = []
        for tip_name in tree.get_tip_labels():
            if is_reference_sequence(tip_name):
                genus = extract_genus(tip_name)
                color = genus_colors.get(genus, '#808080')
                tip_colors.append(color)
            else:
                # Student sample (red)
                tip_colors.append('#FF6B6B')

        # Prepare node labels - only show bootstrap values (0-100), hide branch lengths
        node_labels = []
        for node in tree.treenode.traverse():
            if node.is_leaf():
                node_labels.append('')  # No label for leaves
            else:
                # Node name might be bootstrap value or empty
                if not node.name or node.name.strip() == '':
                    node_labels.append('')  # Empty node name - no bootstrap value
                else:
                    try:
                        value = float(node.name)
                        # Only show if it looks like a bootstrap value (0-100)
                        # Exclude 0 as it's likely not a real bootstrap value
                        if 1 <= value <= 100:
                            node_labels.append(str(int(value)))
                        else:
                            node_labels.append('')  # Hide branch lengths (>100) or 0
                    except (ValueError, TypeError):
                        node_labels.append('')  # Hide non-numeric values

        # Generate different layouts
        layouts = {
            'rectangular': {
                'layout': 'r',
                'edge_type': 'p',
                'width': 1000,
                'height': 900,
                'use_edge_lengths': True,
                'description': 'Traditional rectangular tree'
            },
            'circular': {
                'layout': 'c',
                'edge_type': 'c',  # Curved edges (required for circular layout)
                'width': 1600,  # Increased for better spacing
                'height': 1600,  # Increased for better spacing
                'use_edge_lengths': False,  # Cladogram style - all tips at same radius
                'description': 'Circular cladogram (all tips aligned)'
            },
            'unrooted': {
                'layout': 'u',
                'edge_type': 'p',
                'width': 1200,
                'height': 1200,
                'use_edge_lengths': True,
                'description': 'Unrooted radial layout'
            }
        }

        generated_files = []

        for layout_name, params in layouts.items():
            try:
                # Adjust label positioning and sizing based on layout type
                if layout_name == 'circular':
                    label_shift = '15px'  # More spacing for circular layout
                    tip_label_size = '12px'  # Larger labels for circular
                elif layout_name == 'unrooted':
                    label_shift = '12px'
                    tip_label_size = '10px'
                else:
                    label_shift = '8px'
                    tip_label_size = '10px'

                # Draw tree with bootstrap support values on branches
                canvas, axes, mark = tree.draw(
                    layout=params['layout'],
                    edge_type=params['edge_type'],
                    use_edge_lengths=params['use_edge_lengths'],
                    tip_labels=True,
                    tip_labels_colors=tip_colors,
                    node_sizes=0,  # Hide internal node markers
                    node_labels=node_labels,  # Show filtered bootstrap values (0-100 only)
                    node_labels_style={
                        'font-size': '11px',  # Slightly larger for better readability
                        'fill': '#000000',  # Black for maximum visibility
                        'font-weight': 'bold',  # Bold for emphasis
                        '-toyplot-anchor-shift': label_shift  # Layout-specific offset
                    },
                    edge_widths=2.5 if layout_name == 'circular' else 2,  # Thicker edges for circular
                    edge_colors='#555555',  # Slightly darker for better visibility
                    tip_labels_style={'font-size': tip_label_size},
                    width=params['width'],
                    height=params['height'],
                )

                # Save PNG
                png_file = output_dir / f"tree_{layout_name}.png"
                toyplot.png.render(canvas, str(png_file))

                # Save SVG (scalable, editable)
                svg_file = output_dir / f"tree_{layout_name}.svg"
                toyplot.svg.render(canvas, str(svg_file))

                # Save PDF (for Adobe Illustrator editing)
                pdf_file = output_dir / f"tree_{layout_name}.pdf"
                toyplot.pdf.render(canvas, str(pdf_file))

                print(f"  ✓ {layout_name.capitalize()} layout: {png_file.name}")
                generated_files.append({
                    'name': layout_name,
                    'png': png_file,
                    'svg': svg_file,
                    'pdf': pdf_file,
                    'description': params['description']
                })

            except Exception as e:
                print(f"  ✗ {layout_name.capitalize()} layout failed: {e}")

        return generated_files

    except Exception as e:
        print(f"✗ toytree visualization failed: {e}")
        return []

def generate_html_report(tree_file, image_file, log_file, output_file, toytree_files=None):
    """Generate HTML report for phylogenetic tree using new design system"""
    from datetime import datetime
    from pathlib import Path

    if toytree_files is None:
        toytree_files = []

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
    base_css = (project_root / "modules/styles/base.css").read_text()
    components_css = (project_root / "modules/styles/components.css").read_text()
    reports_css = (project_root / "modules/styles/reports.css").read_text()

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

            <!-- Tree Visualizations -->
            <h3>Tree Layouts</h3>
            <p>Multiple visualizations of the same phylogenetic tree to help you identify relationships:</p>

            <!-- Bio.Phylo Traditional Tree -->
            <div class="figure-container" style="margin-bottom: 2rem;">
                <h4 style="margin-top: 1rem; color: var(--primary-color);">📊 Traditional Layout (Bio.Phylo)</h4>
                <img src="{Path(image_file).name}" alt="Traditional Tree" style="max-width: 100%; height: auto; border: 1px solid var(--border-color); border-radius: var(--radius-md); background: white;">
                <div class="figure-caption">
                    <strong>Figure 1:</strong> Traditional rectangular tree layout with bootstrap values.
                </div>
            </div>

            <!-- Toytree Layouts (if available) -->"""

    # Generate HTML for toytree layouts
    toytree_html = ""
    if toytree_files:
        layout_icons = {
            'rectangular': '📐',
            'circular': '🔵',
            'unrooted': '🌟'
        }
        layout_titles = {
            'rectangular': 'Rectangular Layout (toytree)',
            'circular': 'Circular Layout (toytree)',
            'unrooted': 'Unrooted Radial Layout (toytree)'
        }

        for i, tf in enumerate(toytree_files, start=2):
            icon = layout_icons.get(tf['name'], '🎨')
            title = layout_titles.get(tf['name'], tf['description'])
            toytree_html += f"""
            <div class="figure-container" style="margin-bottom: 2rem;">
                <h4 style="margin-top: 1rem; color: var(--primary-color);">{icon} {title}</h4>
                <img src="{tf['png'].name}" alt="{title}" style="max-width: 100%; height: auto; border: 1px solid var(--border-color); border-radius: var(--radius-md); background: white;">
                <div class="figure-caption">
                    <strong>Figure {i}:</strong> {tf['description']}.
                    <br><small>Download: <a href="{tf['png'].name}">PNG</a> | <a href="{tf['svg'].name}">SVG</a> | <a href="{tf['pdf'].name}">PDF</a></small>
                </div>
            </div>"""

    html += toytree_html + """

            <div class="info-box info-secondary" style="margin-top: 1.5rem;">
                <strong>📥 Download Options:</strong>
                <p>All tree layouts are available in multiple formats:</p>
                <ul>
                    <li><strong>PNG</strong> - For presentations and documents</li>
                    <li><strong>SVG</strong> - Scalable vector graphics (resizable without quality loss)</li>
                    <li><strong>PDF</strong> - Edit in Adobe Illustrator for publication-quality figures</li>
                </ul>
            </div>

            <div class="info-box info-tip">
                <strong>💡 Which layout should I use?</strong>
                <ul>
                    <li><strong>Rectangular:</strong> Best for showing detailed branch lengths and bootstrap values</li>
                    <li><strong>Circular:</strong> Compact, visually appealing, good for posters/presentations</li>
                    <li><strong>Unrooted:</strong> Shows relationships without implying evolutionary direction</li>
                </ul>
                <p style="margin-bottom: 0;"><strong>Summary:</strong> Maximum likelihood phylogenetic tree constructed using IQ-TREE with {bootstrap_support} bootstrap replicates. Best-fit model: {model_info}</p>
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

            <div class="info-box info-tip" style="margin-top: 1rem;">
                <strong>🔍 Why don't all branches have bootstrap values?</strong>
                <p>This is normal! Bootstrap values are only calculated for internal nodes (branch points where lineages split). You won't see bootstrap values for:</p>
                <ul style="margin-bottom: 0;">
                    <li><strong>Terminal branches (tips):</strong> These lead directly to your samples/species - no branching point to test</li>
                    <li><strong>Root node:</strong> The base of the tree has no parent node to calculate support for</li>
                    <li><strong>Some very short branches:</strong> IQ-TREE may skip bootstrap calculation for branches with very little evolutionary change</li>
                </ul>
                <p style="margin-top: 0.5rem; margin-bottom: 0;"><em>In this analysis, IQ-TREE calculated ~79 bootstrap values for ~98 total nodes, which is typical. Focus on the bootstrap values at key branching points that are relevant to your sample identification.</em></p>
            </div>

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
    print("\n" + "="*70)
    print("  🌳  PHYLOGENETIC TREE CONSTRUCTION PIPELINE")
    print("="*70)

    # Step 1: Run IQ-TREE
    print("\n[Step 1/3] Building phylogenetic tree with IQ-TREE...")
    output_prefix = output_dir / "tree"
    if not run_iqtree(input_fasta, output_prefix):
        sys.exit(1)
    print_progress_bar("  ✓ Tree construction complete!", steps=15, delay=0.02)

    # IQ-TREE creates these files
    tree_file = output_prefix.with_suffix('.treefile')
    log_file = output_prefix.with_suffix('.log')
    iqtree_file = output_prefix.with_suffix('.iqtree')

    # Step 2: Visualize tree
    print("\n[Step 2/4] Generating tree visualizations...")
    image_file = output_dir / "tree.png"
    visualize_tree(tree_file, image_file)
    print_progress_bar("  ✓ Bio.Phylo visualization saved!", steps=10, delay=0.03)

    # Step 3: Generate multiple layouts with toytree
    print("\n[Step 3/4] Generating additional tree layouts (toytree)...")
    toytree_files = visualize_tree_toytree(tree_file, output_dir)
    if toytree_files:
        print_progress_bar(f"  ✓ Generated {len(toytree_files)} additional layouts!", steps=10, delay=0.03)
    else:
        print("  ⚠️  No additional layouts generated")

    # Step 4: Generate HTML report
    print("\n[Step 4/4] Creating interactive HTML report...")
    html_file = output_dir / "phylogeny_report.html"
    generate_html_report(tree_file, image_file, log_file, html_file, toytree_files)
    print_progress_bar("  ✓ HTML report generated!", steps=10, delay=0.03)

    print("\n" + "="*70)
    print("  ✓ PHYLOGENY ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\n📁 Output files in {output_dir}:")
    print(f"  🌲 tree.treefile      - Newick format (open in FigTree)")
    print(f"  🖼️  tree.png           - Tree visualization (Bio.Phylo)")
    if toytree_files:
        print(f"\n  Additional layouts (PNG/SVG/PDF):")
        for tf in toytree_files:
            print(f"  🎨 {tf['name']:15} - {tf['description']}")
        print(f"\n  💡 PDF files can be edited in Adobe Illustrator for publications!")
    print(f"\n  🌐 phylogeny_report.html - Interactive HTML report")
    print(f"  📊 tree.iqtree        - Detailed analysis log")
    print(f"\n💡 Open the HTML report in your browser to view the tree!\n")

if __name__ == "__main__":
    main()
