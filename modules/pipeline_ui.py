#!/usr/bin/env python3
"""
Pipeline UI - Beautiful terminal output using Rich

Provides consistent, colorful terminal output for all pipeline steps.
Uses Rich library for panels, progress bars, and styled text.
"""

import sys
from typing import Optional

try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
    from rich.text import Text
    from rich.box import DOUBLE, ROUNDED, HEAVY
    from rich import print as rprint
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False

# Create console instance with force_terminal=True
# This ensures colors work even when output is piped through tee
console = Console(force_terminal=True) if RICH_AVAILABLE else None

# Dracula theme colors
COLORS = {
    'purple': '#bd93f9',
    'cyan': '#8be9fd',
    'green': '#50fa7b',
    'orange': '#ffb86c',
    'pink': '#ff79c6',
    'yellow': '#f1fa8c',
    'red': '#ff5555',
    'white': '#f8f8f2',
    'gray': '#6272a4',
    'bg': '#282a36'
}


# ============================================================================
# ASCII Art Banners
# ============================================================================

MAIN_BANNER = """
[bold purple]╔══════════════════════════════════════════════════════════════════════════════╗[/]
[bold purple]║[/]                                                                              [bold purple]║[/]
[bold purple]║[/]      [bold cyan]██████╗ ███╗   ██╗ █████╗     ██████╗  █████╗ ██████╗  ██████╗ ██████╗[/]  [bold purple]║[/]
[bold purple]║[/]      [bold cyan]██╔══██╗████╗  ██║██╔══██╗    ██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔═══██╗[/] [bold purple]║[/]
[bold purple]║[/]      [bold cyan]██║  ██║██╔██╗ ██║███████║    ██████╔╝███████║██████╔╝██║     ██║   ██║[/] [bold purple]║[/]
[bold purple]║[/]      [bold cyan]██║  ██║██║╚██╗██║██╔══██║    ██╔══██╗██╔══██║██╔══██╗██║     ██║   ██║[/] [bold purple]║[/]
[bold purple]║[/]      [bold cyan]██████╔╝██║ ╚████║██║  ██║    ██████╔╝██║  ██║██║  ██║╚██████╗╚██████╔╝[/] [bold purple]║[/]
[bold purple]║[/]      [bold cyan]╚═════╝ ╚═╝  ╚═══╝╚═╝  ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚═════╝[/]  [bold purple]║[/]
[bold purple]║[/]                                                                              [bold purple]║[/]
[bold purple]║[/]                     [bold white]DNA Barcoding Analysis Pipeline[/]                          [bold purple]║[/]
[bold purple]║[/]                                                                              [bold purple]║[/]
[bold purple]║[/]                🦟 [bold magenta]Cosme Lab[/] • [dim]Fall 2025 • UCR Dept. of Entomology[/]           [bold purple]║[/]
[bold purple]║[/]                                                                              [bold purple]║[/]
[bold purple]╚══════════════════════════════════════════════════════════════════════════════╝[/]
"""


WORKFLOW_DIAGRAM = """
[bold cyan]┌──────────────────────────────────────────────────────────────────────────────┐[/]
[bold cyan]│[/]                         [bold white]ANALYSIS WORKFLOW[/]                                    [bold cyan]│[/]
[bold cyan]├──────────────────────────────────────────────────────────────────────────────┤[/]
[bold cyan]│[/]                                                                              [bold cyan]│[/]
[bold cyan]│[/]   [dim].ab1 FILES[/]        [bold green]STEP 1[/]         [bold green]STEP 2[/]          [bold green]STEP 3[/]         [bold green]STEP 4[/]    [bold cyan]│[/]
[bold cyan]│[/]  [dim]┌──────────┐[/]    [cyan]┌──────────┐[/]   [cyan]┌──────────┐[/]   [cyan]┌──────────┐[/]   [cyan]┌──────────┐[/]  [bold cyan]│[/]
[bold cyan]│[/]  [dim]│ ████████ │[/]    [cyan]│    QC    │[/]   [cyan]│ CONSENSUS│[/]   [cyan]│  ALIGN   │[/]   [cyan]│   TREE   │[/]  [bold cyan]│[/]
[bold cyan]│[/]  [dim]│ ████████ │[/][white]───▶[/][cyan]│  FILTER  │[/][white]──▶[/][cyan]│   F + R  │[/][white]──▶[/][cyan]│  MAFFT   │[/][white]──▶[/][cyan]│ IQ-TREE  │[/]  [bold cyan]│[/]
[bold cyan]│[/]  [dim]│ ████████ │[/]    [cyan]│  Phred   │[/]   [cyan]│  merge   │[/]   [cyan]│          │[/]   [cyan]│          │[/]  [bold cyan]│[/]
[bold cyan]│[/]  [dim]└──────────┘[/]    [cyan]└──────────┘[/]   [cyan]└──────────┘[/]   [cyan]└──────────┘[/]   [cyan]└──────────┘[/]  [bold cyan]│[/]
[bold cyan]│[/]  [dim]Chromatograms[/]   [dim]Quality Ctrl[/]   [dim]Combine reads[/]  [dim]Multi-align[/]   [dim]Phylogenetics[/]  [bold cyan]│[/]
[bold cyan]│[/]                                                                    [white]│[/]         [bold cyan]│[/]
[bold cyan]│[/]              [white]┌─────────────────────────────────────────────────────┘[/]         [bold cyan]│[/]
[bold cyan]│[/]              [white]│[/]                                                               [bold cyan]│[/]
[bold cyan]│[/]              [white]▼[/]                                                               [bold cyan]│[/]
[bold cyan]│[/]        [cyan]┌──────────┐[/]   [cyan]┌──────────┐[/]                                          [bold cyan]│[/]
[bold cyan]│[/]        [cyan]│  BLAST   │[/]   [cyan]│   LAB    │[/]                                          [bold cyan]│[/]
[bold cyan]│[/]        [cyan]│   NCBI   │[/][white]──▶[/][cyan]│  REPORT  │[/]                                          [bold cyan]│[/]
[bold cyan]│[/]        [cyan]│ database │[/]   [cyan]│          │[/]                                          [bold cyan]│[/]
[bold cyan]│[/]        [cyan]└──────────┘[/]   [cyan]└──────────┘[/]                                          [bold cyan]│[/]
[bold cyan]│[/]          [bold green]STEP 5[/]         [bold green]STEP 6[/]                                              [bold cyan]│[/]
[bold cyan]│[/]        [dim]Species ID[/]    [dim]Visualizations[/]                                         [bold cyan]│[/]
[bold cyan]│[/]                                                                              [bold cyan]│[/]
[bold cyan]└──────────────────────────────────────────────────────────────────────────────┘[/]
"""


TEAM_BATTLE = """
[bold orange]╔══════════════════════════════════════════════════════════════════════════════╗[/]
[bold orange]║[/]                                                                              [bold orange]║[/]
[bold orange]║[/]                    [bold white]⚔️  TEAM CHALLENGE: SPIN vs MAGNET  ⚔️[/]                     [bold orange]║[/]
[bold orange]║[/]                                                                              [bold orange]║[/]
[bold orange]╠═══════════════════════════════════╦══════════════════════════════════════════╣[/]
[bold orange]║[/]                                   [bold orange]║[/]                                          [bold orange]║[/]
[bold orange]║[/]         [bold purple]TEAM SPIN[/]                 [bold orange]║[/]           [bold yellow]TEAM MAGNET[/]                    [bold orange]║[/]
[bold orange]║[/]       [dim](Column-based)[/]              [bold orange]║[/]         [dim](Bead-based)[/]                     [bold orange]║[/]
[bold orange]║[/]                                   [bold orange]║[/]                                          [bold orange]║[/]
[bold orange]║[/]    [purple]┌───────────────────┐[/]          [bold orange]║[/]      [yellow]┌───────────────────┐[/]               [bold orange]║[/]
[bold orange]║[/]    [purple]│   ◯◯◯◯◯◯◯◯◯◯     │[/]          [bold orange]║[/]      [yellow]│   ●●●●●●●●●●     │[/]               [bold orange]║[/]
[bold orange]║[/]    [purple]│   Spin columns    │[/]          [bold orange]║[/]      [yellow]│   Magnetic beads  │[/]               [bold orange]║[/]
[bold orange]║[/]    [purple]└───────────────────┘[/]          [bold orange]║[/]      [yellow]└───────────────────┘[/]               [bold orange]║[/]
[bold orange]║[/]                                   [bold orange]║[/]                                          [bold orange]║[/]
[bold orange]║[/]    [bold purple]JR, HV, TW, JM, WL[/]             [bold orange]║[/]      [bold yellow]MA, BR, WA, KG, JA[/]                  [bold orange]║[/]
[bold orange]║[/]                                   [bold orange]║[/]                                          [bold orange]║[/]
[bold orange]╚═══════════════════════════════════╩══════════════════════════════════════════╝[/]
"""


COMPLETION_BANNER = """
[bold green]╔══════════════════════════════════════════════════════════════════════════════╗[/]
[bold green]║[/]                                                                              [bold green]║[/]
[bold green]║[/]                        [bold white]✓ ANALYSIS COMPLETE! ✓[/]                                [bold green]║[/]
[bold green]║[/]                                                                              [bold green]║[/]
[bold green]╠══════════════════════════════════════════════════════════════════════════════╣[/]
[bold green]║[/]                                                                              [bold green]║[/]
[bold green]║[/]   Your results are ready in the [cyan]results/[/] directory:                          [bold green]║[/]
[bold green]║[/]                                                                              [bold green]║[/]
[bold green]║[/]   [bold cyan]📊 01_qc/[/]          Quality control report                                  [bold green]║[/]
[bold green]║[/]   [bold cyan]🧬 02_consensus/[/]   Consensus sequences                                     [bold green]║[/]
[bold green]║[/]   [bold cyan]📐 03_alignment/[/]   Sequence alignment                                      [bold green]║[/]
[bold green]║[/]   [bold cyan]🌳 04_phylogeny/[/]   Phylogenetic tree                                       [bold green]║[/]
[bold green]║[/]   [bold cyan]🦟 05_blast/[/]       Species identification                                  [bold green]║[/]
[bold green]║[/]   [bold cyan]📈 06_lab/[/]         Lab data visualizations                                 [bold green]║[/]
[bold green]║[/]                                                                              [bold green]║[/]
[bold green]║[/]   [bold white]Open the HTML reports in your browser to explore![/]                          [bold green]║[/]
[bold green]║[/]                                                                              [bold green]║[/]
[bold green]╠══════════════════════════════════════════════════════════════════════════════╣[/]
[bold green]║[/]                                                                              [bold green]║[/]
[bold green]║[/]                      🦟 [bold magenta]Cosme Lab[/] • [dim]DNA Barcoding Pipeline[/]                     [bold green]║[/]
[bold green]║[/]                           [dim]Fall 2025 • UCR Dept. of Entomology[/]                    [bold green]║[/]
[bold green]║[/]                                                                              [bold green]║[/]
[bold green]╚══════════════════════════════════════════════════════════════════════════════╝[/]
"""


ENVIRONMENT_CHOICE = """
[bold yellow]╭─────────────────────────────────────────────────────────────────────╮[/]
[bold yellow]│[/]                    [bold white]CHOOSE YOUR ENVIRONMENT[/]                          [bold yellow]│[/]
[bold yellow]├─────────────────────────────────────────────────────────────────────┤[/]
[bold yellow]│[/]                                                                     [bold yellow]│[/]
[bold yellow]│[/]   [bold cyan]OPTION A: GitHub Codespaces[/] [green](No Installation!)[/]                    [bold yellow]│[/]
[bold yellow]│[/]   [dim]──────────────────────────────────────────────[/]                    [bold yellow]│[/]
[bold yellow]│[/]   • Click "Code" → "Codespaces" → "Create codespace"                [bold yellow]│[/]
[bold yellow]│[/]   • Everything runs in the cloud                                    [bold yellow]│[/]
[bold yellow]│[/]   • Use: [cyan]./tutorial-cs.sh[/] and [cyan]./run-analysis-cs.sh[/]                  [bold yellow]│[/]
[bold yellow]│[/]                                                                     [bold yellow]│[/]
[bold yellow]│[/]   [bold purple]OPTION B: Local Docker[/]                                            [bold yellow]│[/]
[bold yellow]│[/]   [dim]──────────────────────────────────────────────[/]                    [bold yellow]│[/]
[bold yellow]│[/]   • Install Docker Desktop on your computer                         [bold yellow]│[/]
[bold yellow]│[/]   • Works on Mac (Intel/M1), Windows, Linux                         [bold yellow]│[/]
[bold yellow]│[/]   • Use: [purple]./tutorial.sh[/] and [purple]./run-analysis.sh[/]                        [bold yellow]│[/]
[bold yellow]│[/]                                                                     [bold yellow]│[/]
[bold yellow]╰─────────────────────────────────────────────────────────────────────╯[/]
"""


# ============================================================================
# Step-specific banners
# ============================================================================

STEP_BANNERS = {
    1: """
[bold green]┌─────────────────────────────────────────────────────────────┐[/]
[bold green]│[/]  [bold white]STEP 1: QUALITY CONTROL[/] 🔬                                 [bold green]│[/]
[bold green]├─────────────────────────────────────────────────────────────┤[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]   [cyan]Input:[/]  .ab1 chromatogram files                            [bold green]│[/]
[bold green]│[/]   [cyan]Check:[/]  Phred quality scores (Q20+ required)               [bold green]│[/]
[bold green]│[/]   [cyan]Check:[/]  Sequence length (>200bp required)                  [bold green]│[/]
[bold green]│[/]   [cyan]Output:[/] passed_sequences.fasta (trimmed)                   [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]   [green]▓▓▓▓▓▓▓▓▓▓[/][yellow]▓▓▓[/][red]░░░[/]  [dim]Quality Score Profile[/]                [bold green]│[/]
[bold green]│[/]   [green]ACGTACGT[/][yellow]NNN[/][red]...[/]   [dim]Low-quality ends trimmed[/]              [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]└─────────────────────────────────────────────────────────────┘[/]
""",
    2: """
[bold green]┌─────────────────────────────────────────────────────────────┐[/]
[bold green]│[/]  [bold white]STEP 2: CONSENSUS SEQUENCES[/] 🧬                              [bold green]│[/]
[bold green]├─────────────────────────────────────────────────────────────┤[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]   [cyan]Forward:[/]  5'──ACGT...──────▶ 3'                            [bold green]│[/]
[bold green]│[/]                    [green]║║║║[/]                                       [bold green]│[/]
[bold green]│[/]   [cyan]Reverse:[/]  3'◀──────...TGCA──5'                             [bold green]│[/]
[bold green]│[/]                    [green]║║║║[/]                                       [bold green]│[/]
[bold green]│[/]   [cyan]Consensus:[/] ════ACGT...════                                 [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]   [dim]Combines F+R for 2x accuracy![/]                              [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]└─────────────────────────────────────────────────────────────┘[/]
""",
    3: """
[bold green]┌─────────────────────────────────────────────────────────────┐[/]
[bold green]│[/]  [bold white]STEP 3: SEQUENCE ALIGNMENT[/] 📐                              [bold green]│[/]
[bold green]├─────────────────────────────────────────────────────────────┤[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]   Sample1:  [cyan]ACGT--ACGTACGT[/]                                   [bold green]│[/]
[bold green]│[/]   Sample2:  [cyan]ACGTACACGTACGT[/]                                   [bold green]│[/]
[bold green]│[/]   Sample3:  [cyan]ACGT--ACGT-CGT[/]                                   [bold green]│[/]
[bold green]│[/]   Ref:      [cyan]ACGTACACGTACGT[/]                                   [bold green]│[/]
[bold green]│[/]             [green]****  **** ***[/]                                   [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]   [dim]MAFFT aligns all sequences[/]                                 [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]└─────────────────────────────────────────────────────────────┘[/]
""",
    4: """
[bold green]┌─────────────────────────────────────────────────────────────┐[/]
[bold green]│[/]  [bold white]STEP 4: PHYLOGENETIC TREE[/] 🌳                               [bold green]│[/]
[bold green]├─────────────────────────────────────────────────────────────┤[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]          [cyan]┌────[/] [purple]Aedes aegypti[/]                              [bold green]│[/]
[bold green]│[/]       [cyan]┌──┤[/]                                                  [bold green]│[/]
[bold green]│[/]       [cyan]│[/]  [cyan]└────[/] [purple]Aedes albopictus[/]                          [bold green]│[/]
[bold green]│[/]    [cyan]───┤[/]                                                     [bold green]│[/]
[bold green]│[/]       [cyan]│[/]  [cyan]┌────[/] [orange3]Culex pipiens[/]                             [bold green]│[/]
[bold green]│[/]       [cyan]└──┤[/]                                                  [bold green]│[/]
[bold green]│[/]          [cyan]└────[/] [bold yellow]YOUR SAMPLE ◄[/]                            [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]   [dim]IQ-TREE builds ML phylogeny with bootstrap support[/]        [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]└─────────────────────────────────────────────────────────────┘[/]
""",
    5: """
[bold green]┌─────────────────────────────────────────────────────────────┐[/]
[bold green]│[/]  [bold white]STEP 5: SPECIES IDENTIFICATION (BLAST)[/] 🦟                  [bold green]│[/]
[bold green]├─────────────────────────────────────────────────────────────┤[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]   Your sequence ──▶ [bold cyan]NCBI GenBank[/]                            [bold green]│[/]
[bold green]│[/]                         [cyan]│[/]                                    [bold green]│[/]
[bold green]│[/]                         [cyan]▼[/]                                    [bold green]│[/]
[bold green]│[/]   [cyan]┌─────────────────────────────┐[/]                            [bold green]│[/]
[bold green]│[/]   [cyan]│[/] [bold white]Top Hit:[/] [purple]Culex pipiens[/]      [cyan]│[/]                            [bold green]│[/]
[bold green]│[/]   [cyan]│[/] [bold white]Identity:[/] [green]99.2%[/]             [cyan]│[/]                            [bold green]│[/]
[bold green]│[/]   [cyan]│[/] [bold white]Accession:[/] KY123456         [cyan]│[/]                            [bold green]│[/]
[bold green]│[/]   [cyan]└─────────────────────────────┘[/]                            [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]   [dim]BLAST searches global database for species match[/]          [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]└─────────────────────────────────────────────────────────────┘[/]
""",
    6: """
[bold green]┌─────────────────────────────────────────────────────────────┐[/]
[bold green]│[/]  [bold white]STEP 6: LAB DATA ANALYSIS[/] 📊                               [bold green]│[/]
[bold green]├─────────────────────────────────────────────────────────────┤[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]   [cyan]DNA Extraction[/] ──▶ [cyan]PCR[/] ──▶ [cyan]Sequencing[/]                     [bold green]│[/]
[bold green]│[/]        [dim]│[/]              [dim]│[/]        [dim]│[/]                             [bold green]│[/]
[bold green]│[/]        [dim]▼[/]              [dim]▼[/]        [dim]▼[/]                             [bold green]│[/]
[bold green]│[/]   [cyan]┌─────────┐[/]   [cyan]┌────────┐[/]  [cyan]┌────────┐[/]                     [bold green]│[/]
[bold green]│[/]   [cyan]│[/] [white]Yield[/]   [cyan]│[/]   [cyan]│[/][white]Success[/] [cyan]│[/]  [cyan]│[/]  [white]QC[/]   [cyan]│[/]                     [bold green]│[/]
[bold green]│[/]   [cyan]│[/] [purple]▓▓▓▓[/][dim]░░░[/] [cyan]│[/]   [cyan]│[/]  [green]60%[/]   [cyan]│[/]  [cyan]│[/] [yellow]40%[/]  [cyan]│[/]                     [bold green]│[/]
[bold green]│[/]   [cyan]└─────────┘[/]   [cyan]└────────┘[/]  [cyan]└────────┘[/]                     [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]│[/]   [dim]Interactive Plotly visualizations of class results[/]        [bold green]│[/]
[bold green]│[/]                                                             [bold green]│[/]
[bold green]└─────────────────────────────────────────────────────────────┘[/]
"""
}


# ============================================================================
# Public API Functions
# ============================================================================

def print_banner():
    """Print the main DNA Barcoding banner."""
    if RICH_AVAILABLE:
        console.print(MAIN_BANNER)
    else:
        print("=" * 70)
        print("   DNA BARCODING ANALYSIS PIPELINE")
        print("   ENTM 201L - UC Riverside")
        print("=" * 70)


def print_workflow():
    """Print the workflow diagram."""
    if RICH_AVAILABLE:
        console.print(WORKFLOW_DIAGRAM)
    else:
        print("\n[.ab1] → [QC] → [CONSENSUS] → [ALIGN] → [TREE] → [BLAST] → [REPORT]\n")


def print_step_banner(step_num: int):
    """Print banner for a specific step."""
    if RICH_AVAILABLE and step_num in STEP_BANNERS:
        console.print(STEP_BANNERS[step_num])
    else:
        step_names = {
            1: "Quality Control",
            2: "Consensus Sequences",
            3: "Sequence Alignment",
            4: "Phylogenetic Tree",
            5: "Species Identification",
            6: "Lab Data Analysis"
        }
        name = step_names.get(step_num, f"Step {step_num}")
        print(f"\n{'='*60}")
        print(f"  STEP {step_num}: {name}")
        print(f"{'='*60}\n")


def print_step_header(step_num: int, total_steps: int, title: str):
    """Print a step header with progress indication."""
    if RICH_AVAILABLE:
        console.print(f"\n[bold cyan]━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[/]")
        console.print(f"[bold white]STEP {step_num} of {total_steps}:[/] [bold green]{title}[/]")
        console.print(f"[bold cyan]━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[/]\n")
    else:
        print(f"\n{'━' * 68}")
        print(f"STEP {step_num} of {total_steps}: {title}")
        print(f"{'━' * 68}\n")


def print_completion():
    """Print completion banner."""
    if RICH_AVAILABLE:
        console.print(COMPLETION_BANNER)
    else:
        print("\n" + "=" * 70)
        print("   ✓ ANALYSIS COMPLETE!")
        print("=" * 70)
        print("\nYour results are in the results/ directory.\n")


def print_team_battle():
    """Print team challenge banner."""
    if RICH_AVAILABLE:
        console.print(TEAM_BATTLE)
    else:
        print("\n" + "=" * 70)
        print("   TEAM CHALLENGE: SPIN vs MAGNET")
        print("=" * 70 + "\n")


def print_environment_choice():
    """Print environment choice banner."""
    if RICH_AVAILABLE:
        console.print(ENVIRONMENT_CHOICE)
    else:
        print("\nCHOOSE YOUR ENVIRONMENT:")
        print("  A: GitHub Codespaces (use -cs.sh scripts)")
        print("  B: Local Docker (use .sh scripts)\n")


def print_success(message: str):
    """Print a success message."""
    if RICH_AVAILABLE:
        console.print(f"[bold green]✓[/] {message}")
    else:
        print(f"✓ {message}")


def print_error(message: str):
    """Print an error message."""
    if RICH_AVAILABLE:
        console.print(f"[bold red]✗[/] {message}")
    else:
        print(f"✗ {message}")


def print_warning(message: str):
    """Print a warning message."""
    if RICH_AVAILABLE:
        console.print(f"[bold yellow]⚠[/] {message}")
    else:
        print(f"⚠ {message}")


def print_info(message: str):
    """Print an info message."""
    if RICH_AVAILABLE:
        console.print(f"[bold cyan]ℹ[/] {message}")
    else:
        print(f"ℹ {message}")


def create_progress_bar(description: str = "Processing"):
    """Create a progress bar context manager."""
    if RICH_AVAILABLE:
        return Progress(
            SpinnerColumn(),
            TextColumn("[bold blue]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            console=console
        )
    else:
        # Return a dummy context manager for non-Rich environments
        class DummyProgress:
            def __enter__(self):
                return self
            def __exit__(self, *args):
                pass
            def add_task(self, description, total):
                return 0
            def update(self, task_id, advance=1):
                pass
        return DummyProgress()


def create_status_table(title: str, data: dict) -> None:
    """Create and print a status table."""
    if RICH_AVAILABLE:
        table = Table(title=title, box=ROUNDED)
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="green")

        for key, value in data.items():
            table.add_row(key, str(value))

        console.print(table)
    else:
        print(f"\n{title}")
        print("-" * 40)
        for key, value in data.items():
            print(f"  {key}: {value}")
        print()


# ============================================================================
# Main - Demo all UI elements
# ============================================================================

if __name__ == "__main__":
    print_banner()
    print_workflow()
    print_environment_choice()

    for i in range(1, 7):
        print_step_banner(i)

    print_team_battle()

    # Demo messages
    print_success("This is a success message")
    print_error("This is an error message")
    print_warning("This is a warning message")
    print_info("This is an info message")

    # Demo table
    create_status_table("Analysis Summary", {
        "Total Sequences": 30,
        "Passed QC": 20,
        "Species Identified": 4
    })

    print_completion()
